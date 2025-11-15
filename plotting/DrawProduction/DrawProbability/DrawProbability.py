import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from matplotlib import rc
from matplotlib.ticker import MaxNLocator
from matplotlib.pyplot import figure
from matplotlib.colors import LinearSegmentedColormap, LogNorm
import pandas as pd
import os, glob

# geometry + detection
NPOT = 1e20
A_DET = 3
N_GAMMA = 2.5e5
R_DET = 0.5
L_BASE = 40.0
THETA_CUT = np.arctan(R_DET / L_BASE)

# --- Proton bremsstrahlung (V2) helpers ---
def _effective_lumi_pb(N_POT=1.0e20, rho=7.87, L=16.8, A=55.845):
    """
    L_eff [pb^-1] = N_POT * (rho*L/A) * N_A  [1/cm^2]  / 1e36
    """
    NA = 6.02214076e23  # /mol
    nuclei_per_cm2 = (rho * L / A) * NA
    L_eff_cm2 = N_POT * nuclei_per_cm2
    return L_eff_cm2 / 1e36

def brem_v2_to_grid(m_grid_masses, m1d, y1d_baseline, charges, a=3, N_gamma=2.5e5):
    """
    Interpolate 'baseline' (no ε, no Poisson) onto mass grid, then broadcast with
    detection factor: (1 - exp(-Nγ ε^2))^a * ε^2
    """
    if m1d.size == 0:
        return np.zeros((charges.shape[0], m_grid_masses.size), dtype=float)
    out = np.zeros_like(m_grid_masses, dtype=float)
    pos = (y1d_baseline > 0) & (m1d > 0)
    if np.any(pos):
        xm = np.log10(m1d[pos]); ym = np.log10(y1d_baseline[pos])
        xq = np.log10(np.clip(m_grid_masses, m1d[pos].min()*0.9999, m1d[pos].max()*1.0001))
        yq = np.interp(xq, xm, ym, left=-np.inf, right=-np.inf)
        valid = np.isfinite(yq)
        out[valid] = 10.0**yq[valid]
    det = (1.0 - np.exp(-N_gamma * charges**2))**a * (charges**2)
    return det[:, None] * out[None, :]# shape: [n_eps, n_m]

def _brem_v2_folder_baseline_single_or_shared(
        directory,
        pattern="Brem_120GeV_*.txt",
        det_angle=np.arctan(0.5/40.0),
        lambda_idx=1,                 # 0→Λp=1.0, 1→1.5, 2→2.0
        N_POT=1.0e20, rho=7.87, L=16.8, A=55.845,
):
    """
    Read SINGLE-OR-SHARED files:
      log10(theta), log10(p/GeV), sigma_L1, sigma_L1p5, sigma_L2   [pb/bin]

    Baseline returned is *events* (no ε^2, no Poisson), and *no extra ×2*.
    Saeid’s note: two χ in same bar already accounted (no additional factor).
    """
    files = sorted(glob.glob(os.path.join(directory, pattern)))
    if not files:
        return np.array([]), np.array([])

    L_eff_pb = _effective_lumi_pb(N_POT=N_POT, rho=rho, L=L, A=A)
    log10_theta_cut = np.log10(det_angle)

    masses, y_events = [], []
    for f in files:
        try:
            mchi = float(os.path.basename(f).split('_')[-1].replace('.txt',''))
        except Exception:
            continue

        arr = np.loadtxt(f, comments="#")
        if arr.ndim == 1: arr = arr[None, :]
        if arr.shape[1] < 5:  # expect [log10θ, log10p, sigma_L1, sigma_L1p5, sigma_L2]
            continue

        log10theta = arr[:,0]
        sigmas = arr[:,2:5]
        mask = (log10theta < log10_theta_cut)
        sigma_sum_pb = sigmas[mask, lambda_idx].sum()

        # SINGLE-OR-SHARED: baseline is *events* = L_eff * Σσ (no ×2, no ε^2)
        y_events.append(L_eff_pb * sigma_sum_pb)
        masses.append(mchi)

    masses = np.asarray(masses); y_events = np.asarray(y_events)
    idx = np.argsort(masses)
    return masses[idx], y_events[idx]

# --- use SINGLE-OR-SHARED baseline only ---
BREM_SHARED_DIR = "/Users/leobailloeul/Documents/coding/decaysimulation/single_hit_or_shared_bar"

def _load_ageo_aligned(path, fileMass):
    """
    Load an acceptance file and interpolate its acceptance vs mass onto fileMass.
    Supports 2+ column files (col0=mass, others=acceptances).
    """
    arr = np.loadtxt(path)
    if arr.ndim == 1:
        raise ValueError(f"{path} has 1 column but length != fileMass; cannot align without mass column.")
    x = arr[:, 0].astype(float)

    # Choose which acceptance column to use
    # Special case for Eta DARK files that store acceptance in col 2
    if ("Eta" in os.path.basename(path)) and (arr.shape[1] >= 3) and ("dark" in path.lower()):
        y = arr[:, 2].astype(float)
    else:
        y = arr[:, 1].astype(float)

    # Interpolate to fileMass (fill=0 outside range)
    y_on_grid = _interp_to(x, y, fileMass, fill=0.0)
    return y_on_grid  # shape matches fileMass

def mesonProduction_2d(
        fileNamePi, fileNameEta, fileNameRho, fileNamePhi, fileNameOmega, fileNameJpsi,
        fileMass, NPOT, charges, a, N_gamma=2.5e5, SHIP=False
):
    # constants
    alpha = 1.0/137.0
    m_e   = 0.00051
    m_pi, m_eta, m_rho, m_omega, m_phi, m_jpsi, m_ups = 0.135, 0.548, 0.775, 0.782, 1.019, 3.1, 9.46
    # meson / NPOT (PYTHIA)
    c_pi, c_eta, c_rho, c_omega, c_phi, c_jpsi, c_ups = 3.98, 0.51, 0.24, 0.24, 4.9e-3, 3.81e-5, 2.5e-9
    # BR to e+e-
    br_pi, br_eta = 0.98, 0.39
    br_rho, br_omega, br_phi = 4.72e-5, 7.28e-5, 2.95e-4
    br_jpsi, br_ups = 0.05971, 0.0238

    # load + align acceptances to fileMass (→ shape (nx,))
    ageo_pi  = _load_ageo_aligned(fileNamePi,   fileMass)
    ageo_eta = _load_ageo_aligned(fileNameEta,  fileMass)
    ageo_rho = _load_ageo_aligned(fileNameRho,  fileMass)
    ageo_phi = _load_ageo_aligned(fileNamePhi,  fileMass)
    ageo_om  = _load_ageo_aligned(fileNameOmega,fileMass)
    ageo_j   = _load_ageo_aligned(fileNameJpsi, fileMass)
    if SHIP is True:
        ageo_ups = np.full(fileMass.size, 0.011338)
    else:
        ageo_ups = np.full(fileMass.size, 0.010276)  # your fixed Υ acceptance


    # helpers
    def I2(x, y):
        return ((1 + 2*x) * np.sqrt(1 - 4*x)) / ((1 + 2*y) * np.sqrt(1 - 4*y))
    def _I3_integrand(z, x):
        return 2/(3*3.14) * np.sqrt(1 - 4*x/z) * (1 - z)**3 * (2*x + z) / (z**2)
    def I3(x):
        from scipy.integrate import quad
        return quad(_I3_integrand, 4*x, 1, args=(x))[0]
    vI3 = np.vectorize(I3)

    m = fileMass  # (nx,)

    # masks
    ms_pi   = m < m_pi/2
    ms_eta  = m < m_eta/2
    ms_rho  = m < m_rho/2
    ms_om   = m < m_omega/2
    ms_phi  = m < m_phi/2
    ms_jpsi = m < m_jpsi/2
    ms_ups  = m < m_ups/2

    # per-channel baselines (nx,)
    base_pi   = np.zeros_like(m)
    base_eta  = np.zeros_like(m)
    base_rho  = np.zeros_like(m)
    base_om   = np.zeros_like(m)
    base_phi  = np.zeros_like(m)
    base_jpsi = np.zeros_like(m)
    base_ups  = np.zeros_like(m)

    if np.any(ms_pi):
        x = (m[ms_pi]**2) / (m_pi**2)
        base_pi[ms_pi] = ageo_pi[ms_pi] * (2 * c_pi * br_pi * alpha * vI3(x))
    if np.any(ms_eta):
        x = (m[ms_eta]**2) / (m_eta**2)
        base_eta[ms_eta] = ageo_eta[ms_eta] * (2 * c_eta * br_eta * alpha * vI3(x))
    if np.any(ms_rho):
        x = (m[ms_rho]**2) / (m_rho**2); y = (m_e**2) / (m_rho**2)
        base_rho[ms_rho] = ageo_rho[ms_rho] * (2 * c_rho * br_rho * I2(x, y))
    if np.any(ms_om):
        x = (m[ms_om]**2) / (m_omega**2); y = (m_e**2) / (m_omega**2)
        base_om[ms_om] = ageo_om[ms_om] * (2 * c_omega * br_omega * I2(x, y))
    if np.any(ms_phi):
        x = (m[ms_phi]**2) / (m_phi**2); y = (m_e**2) / (m_phi**2)
        base_phi[ms_phi] = ageo_phi[ms_phi] * (2 * c_phi * br_phi * I2(x, y))
    if np.any(ms_jpsi):
        x = (m[ms_jpsi]**2) / (m_jpsi**2); y = (m_e**2) / (m_jpsi**2)
        base_jpsi[ms_jpsi] = ageo_j[ms_jpsi] * (2 * c_jpsi * br_jpsi * I2(x, y))
    if np.any(ms_ups):
        x = (m[ms_ups]**2) / (m_ups**2); y = (m_e**2) / (m_ups**2)
        base_ups[ms_ups] = ageo_ups[ms_ups] * (2 * c_ups * br_ups * I2(x, y))

    baseline_1d = (base_pi + base_eta + base_rho + base_om + base_phi + base_jpsi + base_ups)  # (nx,)

    # detection factor and broadcast to 2D (ny, nx)
    det = (1.0 - np.exp(-N_gamma * charges**2))**a * (charges**2)
    return NPOT * det[:, None] * baseline_1d[None, :]

def _load_xy(path):
    arr = np.loadtxt(path)
    if arr.ndim == 1:                 # single column → y only
        return None, arr.astype(float)
    return arr[:,0].astype(float), arr[:,1].astype(float)  # x, y

def _interp_to(x, y, xgrid, fill=0.0):
    if x is None:
        if y.size != xgrid.size:
            raise ValueError("No mass column and length != fileMass.")
        return y
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    o = np.argsort(x)
    return np.interp(xgrid, x[o], y[o], left=fill, right=fill)
def dyProduction_2d(fileNamecross, filenamedy, fileMass, NPOT, charges, a,
                    N_gamma=2.5e5, totalCrossSection=300e-3):
    """
    Produces a 2D array of DY yields vs (epsilon, mass)
    using a shared fileMass grid and matching lengths.
    """
    # Load cross section and acceptance arrays
    sigma = np.loadtxt(fileNamecross)   # shape (n_mass,) or (n_mass, 1)
    ageo  = np.loadtxt(filenamedy)[:, 1]  # take only acceptance column

    # Sanity check
    if sigma.ndim > 1:
        sigma = sigma[:, 1]  # assume column 1 is σ(m)
    if sigma.size != fileMass.size or ageo.size != fileMass.size:
        raise ValueError("Cross-section or acceptance arrays don't match fileMass size.")

    # Build 1D baseline (yields vs m)
    baseline_1d = NPOT * (sigma * 1e-12) / totalCrossSection * ageo

    # Detection factor vs epsilon
    det = (1.0 - np.exp(-N_gamma * charges**2))**a * (charges**2)

    # Broadcast to 2D
    return det[:, None] * baseline_1d[None, :]

# mass SHIP vs Dark
NPOTSHiP = 2e20
NPOT = 1e20

# masses SHIP vs Dark
# mass = np.loadtxt('/Users/leobailloeul/Documents/coding/decaysimulation/decay/sensitivity-plot/mass_ship.txt')
mass = np.loadtxt('/Users/leobailloeul/Documents/coding/decaysimulation/plotting/sensitivity-plot/mship_values.txt')
massSHiP = np.loadtxt('/Users/leobailloeul/Documents/coding/decaysimulation/plotting/sensitivity-plot/mchi_values.txt')
base = '/Users/leobailloeul/Documents/coding/decaysimulation/plotting/sensitivity-plot/'
charge = np.logspace(-5, 1, 500)
# masses, charges = np.meshgrid(mass, charge)

mesonDecay_dune = mesonProduction_2d(
    base+'total_efficiency_outputdecayPionDARK.txt',
    base+'total_efficiency_outputdecayEtaDARK.txt',
    base+'total_efficiency_output2bodydecay-rhoDARK.txt',
    base+'total_efficiency_output2bodydecay-phiDARK.txt',
    base+'total_efficiency_output2bodydecay-omegaDARK.txt',
    base+'total_efficiency_output2bodydecay-jsiDARK.txt',
    fileMass=mass, NPOT=NPOT, charges=charge, a=3,
    )
mesonDecay_SHiP = mesonProduction_2d(base + 'total_efficiency_outputdecayPion_SHIP-1m.txt',
    base + 'total_efficiency_outputdecayEta_SHIP-1m.txt',
    base + 'total_efficiency_output2bodydecay-rho-SHIP-1m.txt',
    base + 'total_efficiency_output2bodydecay-phi-SHIP-1m.txt',
    base + 'total_efficiency_output2bodydecay-omega-SHIP-1m.txt',
    base + 'total_efficiency_output2bodydecay-jsi-SHIP-1m.txt',
    fileMass=massSHiP, NPOT=NPOTSHiP, charges=charge, a=3, SHIP=True
    )

# get DY production
dy_dune = dyProduction_2d(
    fileNamecross=base+'dy_cross_091525.txt',
    filenamedy=base+'ageo_dy_dark_0.5m.txt',
    fileMass=mass, NPOT=1e20, charges=charge, a=3,
    totalCrossSection=300e-3,
)

dy_dune_new = dyProduction_2d(
    fileNamecross=base+'dy_cross_091525.txt',
    filenamedy=base+'ageo_dy_dark_0.5m.txt',
    fileMass=mass, NPOT=1e20, charges=charge, a=3,
    totalCrossSection=40e-3,
)

dy_SHiP = dyProduction_2d(
    fileNamecross=base+'dy_cross_ship_final.txt',
    filenamedy=base+'ageo_dy_ship_0.5m_final.txt',
    fileMass=massSHiP, NPOT=2e20, charges=charge, a=3,
    totalCrossSection=300e-3,
)

dy_SHiP_new = dyProduction_2d(
    fileNamecross=base+'dy_cross_ship_final.txt',
    filenamedy=base+'ageo_dy_ship_0.5m_final.txt',
    fileMass=massSHiP, NPOT=2e20, charges=charge, a=3,
    totalCrossSection=40e-3,
)

# dy_dune_2layers = dyProduction(base+'dy_cross_091525.txt', base+'ageo_dy_dark_0.5m.txt', NPOT, charges, 2)
# dy_dune = dyProduction(base+'efficiency_dy_dark_0.5m.txt', base+'dy_cross_dark.txt', NPOT, charges, 2)

total_dune = mesonDecay_dune.copy()
total_dune += dy_dune_new

total_dune_new = mesonDecay_dune.copy()
total_dune_new += dy_dune_new

total_SHiP = mesonDecay_SHiP.copy()
total_SHiP += dy_SHiP_new

total_SHiP_new = mesonDecay_SHiP.copy()
total_SHiP_new += dy_SHiP_new

# total_dune_2layers = mesonDecay_dune_2layers.copy()
# total_dune_2layers += dy_dune_2layers


m_shared, y_shared0 = _brem_v2_folder_baseline_single_or_shared(
    directory=BREM_SHARED_DIR,
    det_angle=THETA_CUT, lambda_idx=1,  # Λp = 1.5 GeV
    N_POT=NPOT, rho=7.87, L=16.8, A=55.845,
)

# Map brems baseline → 2D (ε grid) with detection factor
brem_2d = brem_v2_to_grid(
    m_grid_masses=mass,          # your sensitivity mass grid (trimmed to L above)
    m1d=m_shared,                # single-or-shared masses
    y1d_baseline=y_shared0,      # single-or-shared baseline events (no ε, no Poisson)
    charges=charge,
    a=A_DET,
    N_gamma=N_GAMMA,
)

# brem_2d_2layers = brem_v2_to_grid(
#     m_grid_masses=mass,          # your sensitivity mass grid (trimmed to L above)
#     m1d=m_shared,                # single-or-shared masses
#     y1d_baseline=y_shared0,      # single-or-shared baseline events (no ε, no Poisson)
#     charges=charge,
#     a=2,
#     N_gamma=N_GAMMA,
# )

# total_dune = total_dune + brem_2d
total_dune_new = total_dune_new + brem_2d
# total_dune_2layers = total_dune_2layers + brem_2d_2layers

print("Generation completed")

figure(figsize=(22, 16), dpi=500)
plt.tick_params(axis='both', which='both', labelsize=25)

plt.contour(mass, charge, total_dune, levels=[3.09], colors = 'skyblue', linewidths = 3)
plt.plot(10, 5, color='skyblue', label='1 year at DarkQuest 3 layers, 315 bars/layer, bkg = 0, meson + Drell-Yan production only')

plt.contour(mass, charge, total_dune_new, levels=[3.09], colors = 'salmon', linewidths = 3)
plt.plot(10, 5, color='salmon', label='1 year at DarkQuest 3 layers, 315 bars/layer, bkg = 0')

# plt.contour(mass, charge, dy_dune_new, levels=[3.09], colors = 'greenyellow', linewidths = 3)
# plt.plot(10, 5, color='greenyellow', label='3 layers, 315 bars/layer, bkg = 0, Updated Drell-Yan with 10x production')

plt.contour(massSHiP, charge, total_SHiP_new, levels=[3.09], colors = 'greenyellow', linewidths = 3)
plt.plot(10, 5, color='greenyellow', label='5 years at SHiP 3 layers, 315 bars/layer, bkg = 0, meson + Drell-Yan production only')

# slac = pd.read_csv('../../experiment-contours/slac.csv',header=None)
# colliders = pd.read_csv('../../experiment-contours/colliders.csv',header=None)
# bebc = pd.read_csv('../../experiment-contours/bebc.csv',header=None)
# charmii = pd.read_csv('../../experiment-contours/charmii.csv',header=None)
# mq_demo = pd.read_csv('../../experiment-contours/mq_demonstrator_sort.csv',header=None)
# argoneut = pd.read_csv('../../experiment-contours/argoneut_sort.csv',header=None)
# lsnd = pd.read_csv('../../experiment-contours/LSND.csv', header=None)
# sensei = pd.read_csv('../../experiment-contours/SENSEI_POLISHED.csv', header=None)
# mq_run3 = pd.read_csv('../../experiment-contours/MilliQanRun3.csv', header=None)

slac_small = pd.read_csv('../../experiment-contours-small/SLACmQ.csv',header=None)
# colliders_small = pd.read_csv('../../experiment-contours-small/.csv',header=None)
MiniBooNE_small = pd.read_csv('../../experiment-contours-small/MiniBooNE.csv',header=None)
bebc_small = pd.read_csv('../../experiment-contours-small/BEBC.csv',header=None)
charmii_small = pd.read_csv('../../experiment-contours-small/CHARMII.csv',header=None)
# mq_demo_small = pd.read_csv('../../experiment-contours-small/mq_demonstrator_sort.csv',header=None)
argoneut_small = pd.read_csv('../../experiment-contours-small/ArgoNeuT.csv',header=None)
lsnd_small = pd.read_csv('../../experiment-contours-small/LSND.csv', header=None)
sensei_small = pd.read_csv('../../experiment-contours-small/SENSEI_POLISHED.csv', header=None)
mq_run3_small = pd.read_csv('../../experiment-contours-small/milliQan.csv', header=None)
lep_small = pd.read_csv('../../experiment-contours-small/LEP.csv', header=None)

#
# colliders_plot = plt.fill_between(colliders[0], colliders[1], 2, label = 'Colliders', alpha=0.5, color = 'mediumseagreen') #'green')
# slac_plot = plt.fill_between(slac[0], slac[1], 2, label = 'SLAC', alpha=0.5, color = 'gold') #'orange')
# bebc_plot = plt.fill_between(bebc[0], bebc[1], 2, label = 'BEBC', alpha=0.5, color = 'gray')
# lsnd_plot = plt.fill_between(lsnd[0]/1000, lsnd[1], 2, label='LSND', alpha=0.5, color = 'skyblue')
# charmii_plot = plt.fill_between(charmii[0], charmii[1], 2, label = 'Charm II', alpha=0.5, color = 'lightgray')
# argonuet_plot = plt.fill_between(argoneut[0], argoneut[1], 2, label = 'ArgoNeuT', alpha=0.5, color = 'cornflowerblue') #color = 'mediumorchid')
# mq_demo_plot = plt.fill_between(mq_demo[0], mq_demo[1], 2, label = 'MilliQan demonstrator', alpha=0.5, color = 'lightcoral')
# sensei_plot = plt.fill_between(sensei[0]/1000, sensei[1], 2, label = 'SENSEI', alpha=0.5, color = '#FFE6CC')
# mq_run3_plot = plt.fill_between(mq_run3[0], mq_run3[1], 2, label = 'MilliQan Run3', alpha=0.5, color = '#E6E6FA')

MiniBooNE_plot = plt.fill_between(MiniBooNE_small[0], MiniBooNE_small[1], 2, label = 'MiniBooNE', alpha=0.5, color = 'mediumseagreen') #'green')
slac_plot = plt.fill_between(slac_small[0], slac_small[1], 2, label = 'SLAC mQ', alpha=0.5, color = 'gold') #'orange')
bebc_plot = plt.fill_between(bebc_small[0], bebc_small[1], 2, label = 'BEBC', alpha=0.5, color = 'gray')
lsnd_plot = plt.fill_between(lsnd_small[0], lsnd_small[1], 2, label='LSND', alpha=0.5, color = 'skyblue')
charmii_plot = plt.fill_between(charmii_small[0], charmii_small[1], 2, label = 'Charm II', alpha=0.5, color = 'lightgray')
argonuet_plot = plt.fill_between(argoneut_small[0], argoneut_small[1], 2, label = 'ArgoNeuT', alpha=0.5, color = 'cornflowerblue') #color = 'mediumorchid')
lep_plot = plt.fill_between(lep_small[0], lep_small[1], 2, label = 'LEP', alpha=0.5, color = 'lightcoral')
sensei_plot = plt.fill_between(sensei_small[0]/1000, sensei_small[1], 2, label = 'SENSEI', alpha=0.5, color = '#FFE6CC')
mq_run3_plot = plt.fill_between(mq_run3_small[0], mq_run3_small[1], 2, label = 'milliQan', alpha=0.5, color = '#E6E6FA')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('$m_{\chi}$ [$\mathrm{GeV}/\mathrm{c}^2$]', fontsize = 35)
plt.xticks(fontsize = 35)
plt.ylabel('$\epsilon=Q/e$', fontsize = 35)
plt.yticks(fontsize = 35)
plt.xlim(0.001, 10)
plt.ylim(0.000005, 1)
plt.legend(loc='upper left', fontsize=18.5)
plt.savefig('limit-plot-DarkQuest-SHiP.pdf')
print("limit plot drawn")

#

# # Corrected second plot for heatmap
# figure(figsize=(22, 16), dpi=250)
# plt.tick_params(axis='both', which='both', labelsize=25)
#
# # Display the heatmap
# cmap = LinearSegmentedColormap.from_list('white_to_blue', ['white', 'darkblue'])
# plt.pcolormesh(mass, charge, total_dune, norm=LogNorm(vmin=1e0, vmax=1e16), shading='auto', cmap=cmap)
# plt.colorbar()
#
# # Adding the contour lines
# plt.contour(mass, charge, total_dune, levels=[10], colors='greenyellow')
# plt.plot(10, 5, color='greenyellow', label='2 layers, 80 bars/layer, bkg = 20')
#
# plt.contour(mass, charge, total_dune, levels=[21], colors='tomato')
# plt.plot(10, 5, color='tomato', label='2 layers, 80 bars/layer, bkg = 100')
#
# plt.contour(mass, charge, total_dune, levels=[29], colors='skyblue')
# plt.plot(10, 5, color='skyblue', label='2 layers, 80 bars/layer, bkg = 200')
#
# plt.contour(mass, charge, total_dune, levels=[45], colors='cyan')
# plt.plot(10, 5, color='cyan', label='2 layers, 80 bars/layer, bkg = 500')
#
# plt.contour(mass, charge, total_dune, levels=[63], colors='orange')
# plt.plot(10, 5, color='orange', label='2 layers, 80 bars/layer, bkg = 1000')
#
# # Setting logarithmic scale for both axes
# plt.xscale('log')
# plt.yscale('log')
#
# # Setting labels and limits
# plt.xlabel('$m_{\chi}$ [$\mathrm{GeV}/\mathrm{c}^2$]', fontsize=35)
# plt.xticks(fontsize=35)
# plt.ylabel('$\epsilon=Q/e$', fontsize=35)
# plt.yticks(fontsize=35)
# plt.xlim(0.01, 8)
# plt.ylim(0.000005, 1)
#
# # Adding the legend
# plt.legend(loc='lower right', fontsize=19.5)
#
# # Saving the heatmap figure
# plt.savefig('heatmap_longquest.png')
# print("heatmap_longquest.png drawn")


# def mesonDecayProduction(fileNamePi, fileNameEta, fileNameRho, fileNamePhi, fileNameOmega, fileNameJpsi, NPOT, mass, charge, a, fileMass, N_gamma=2.5e5): # num of layers
#     # EM constant
#     alpha = 1.0 / 137.0
#     # mass in GeV
#     m_e = 0.00051
#     m_pi = 0.135
#     m_eta = 0.548
#     m_rho = 0.775
#     m_omega = 0.782
#     m_phi = 1.019
#     m_jpsi = 3.1
#     m_upsilon = 9.46
#     # meson / NPOT obtained from PYTHIA
#     c_pi = 3.98
#     c_eta = 0.51
#     c_rho = 0.24
#     c_omega = 0.24
#     c_phi = 4.9e-03
#     c_jpsi = 3.81e-5
#     c_upsilon = 2.5e-9
#     # m -> e+e- branching ratio
#     branch_pi = 0.98
#     branch_eta = 0.39
#     branch_rho = 4.72e-5
#     branch_omega = 7.28e-5
#     branch_phi = 2.95e-4
#     branch_jpsi = 0.05971
#     branch_upsilon = 0.0238
#     # Import data from efficiency files
#     ageo_datapi = np.loadtxt(fileNamePi)
#     ageo_dataeta = np.loadtxt(fileNameEta)
#     ageo_datarho = np.loadtxt(fileNameRho)
#     ageo_dataphi = np.loadtxt(fileNamePhi)
#     ageo_dataomega = np.loadtxt(fileNameOmega)
#     ageo_datajpsi = np.loadtxt(fileNameJpsi)
#     mass_upsilon = fileMass[:]
#     # Separate efficiency and mass
#     ageo_rho = ageo_datarho[:,1]
#     ageo_omega = ageo_dataomega[:,1]
#     ageo_phi = ageo_dataphi[:,1]
#     ageo_pi = ageo_datapi[:,1]
#     if fileNameEta == '/Users/leobailloeul/Documents/coding/decaysimulation/plotting/sensitivity-plot/ageo_darkquest_0.5m.txt':
#         ageo_eta = ageo_dataeta[:,2]
#     else:
#         ageo_eta = ageo_dataeta[:,1]
#     ageo_jpsi = ageo_datajpsi[:,1]
#     # SHIP vs LongQuest
#     # ageo_upsilon = np.full(mass_upsilon.size, 0.010127)
#     ageo_upsilon = np.full(mass_upsilon.size, 0.010276)
#     # 0.000474
#
#     list_ageo = [ageo_rho, ageo_omega, ageo_phi, ageo_pi, ageo_eta, ageo_jpsi, ageo_upsilon]
#
#     # print(fileMass.size)
#     # print(charge.size)
#
#     #padding
#     for i in range(len(list_ageo)):
#         difference = fileMass.size - np.size(list_ageo[i])
#         list_ageo[i] = np.pad(list_ageo[i], (0, difference), 'constant',  constant_values=0)
#         list_ageo[i] = np.tile(list_ageo[i], (int(charge.size / fileMass.size), 1))
#
#     # Re-assign back to the original variables
#     ageo_rho, ageo_omega, ageo_phi, ageo_pi, ageo_eta, ageo_jpsi, ageo_upsilon = list_ageo
#
#     # define phase space integrals
#     def I2(x, y):
#         return ((1 + 2 * x) * np.sqrt(1 - 4 * x)) / ((1 + 2 * y) * np.sqrt(1 - 4 * y))
#
#     def I3_integrand(z, x):
#         return 2 / (3 * 3.14) * np.sqrt(1 - 4 * x / z)  * ((1 - z) ** 3) * (2 * x + z) / (z ** 2)
#
#     def I3(x):
#         return quad(I3_integrand, 4 * x, 1, args=(x))[0]  # integrate
#     v_I3 = np.vectorize(I3) # vectorize
#
#     pi = np.zeros(np.shape(masses))
#     eta = np.zeros(np.shape(masses))
#     phi = np.zeros(np.shape(masses))
#     omega = np.zeros(np.shape(masses))
#     jpsi = np.zeros(np.shape(masses))
#     rho = np.zeros(np.shape(masses))
#     upsilon = np.zeros(np.shape(masses))
#
#     # Dalitz decay
#     pi[mass < m_pi/2] = ageo_pi[mass < m_pi/2] * 2 * c_pi * branch_pi * alpha * v_I3(mass[mass < m_pi/2] ** 2 / m_pi ** 2) * (1 - np.exp(-N_gamma * charge[mass < m_pi / 2] ** 2)) ** a * charge[mass < m_pi / 2] ** 2
#
#     eta[mass < m_eta/2] = ageo_eta[mass < m_eta/2] * 2 * c_eta * branch_eta * alpha * v_I3(mass[mass < m_eta/2] ** 2 / m_eta ** 2) * (1 - np.exp(-N_gamma * charge[mass < m_eta / 2] ** 2)) ** a * charge[mass < m_eta / 2] ** 2
#
#     # Direct decay
#     jpsi[mass < m_jpsi/2] = ageo_jpsi[mass < m_jpsi/2] * 2 * c_jpsi * branch_jpsi * I2( mass[mass < m_jpsi/2] ** 2 / m_jpsi ** 2, m_e ** 2 / m_jpsi ** 2  ) * (1 - np.exp(-N_gamma * charge[mass < m_jpsi / 2] ** 2)) ** a * charge[mass < m_jpsi / 2] ** 2
#
#     rho[mass < m_rho/2] = ageo_rho[mass < m_rho/2] * 2 * c_rho * branch_rho * I2( mass[mass < m_rho/2] ** 2 / m_rho ** 2, m_e ** 2 / m_rho ** 2  ) * (1 - np.exp(-N_gamma * charge[mass < m_rho / 2] ** 2)) ** a * charge[mass < m_rho / 2] ** 2
#
#     omega[mass < m_omega/2] = ageo_omega[mass < m_omega/2] * 2 * c_omega * branch_omega * I2( mass[mass < m_omega/2] ** 2 / m_omega ** 2, m_e ** 2 / m_omega ** 2  ) * (1 - np.exp(-N_gamma * charge[mass < m_omega / 2] ** 2)) ** a * charge[mass < m_omega / 2] ** 2
#
#     phi[mass < m_phi/2] = ageo_phi[mass < m_phi/2] * 2 * c_phi * branch_phi * I2( mass[mass < m_phi/2] ** 2 / m_phi ** 2, m_e ** 2 / m_phi ** 2  ) * (1 - np.exp(-N_gamma * charge[mass < m_phi / 2] ** 2)) ** a * charge[mass < m_phi / 2] ** 2
#
#     upsilon[mass < m_upsilon/2] = ageo_upsilon[mass < m_upsilon/2] * 2 * c_upsilon * branch_upsilon * I2( mass[mass < m_upsilon/2] ** 2 / m_upsilon ** 2, m_e ** 2 / m_upsilon ** 2  )* (1 - np.exp(-N_gamma * charge[mass < m_upsilon / 2] ** 2)) ** a * charge[mass < m_upsilon / 2] ** 2
#
#     sensitivity = upsilon
#
#     return sensitivity * NPOT

# def dyProduction(fileNamecross, filenamedy, fileMass, NPOT, charges, a, N_gamma=2.5e5, totalCrossSection = 300e-3):
#
#     # Import Data from efficiency files
#     print(fileMass.size)
#
#     ageo_datacross = np.loadtxt(fileNamecross)
#     ageo_datady = np.loadtxt(filenamedy)
#
#
#
#     cross_dy = ageo_datacross[:] * 1e-12 # pico barn
#     ageo_dy = ageo_datady[:,1]
#
#     difference = fileMass.size - ageo_dy.size
#     ageo_dy = np.pad(ageo_dy, (0, difference), 'constant', constant_values=0)
#
#     return NPOT * cross_dy / totalCrossSection * ageo_dy * (1 - np.exp(-N_gamma * charges ** 2)) ** a * charges ** 2
#
# print(sensei[0].max(), sensei[0].min(), sensei[1].max(), sensei[1].min())

# mesonDecay_dune = list(mesonDecayProduction(base + 'DQ-total_efficiency_outputdecayPion1m.txt',
#                               base + 'DQ-total_efficiency_outputdecayEta1m.txt',
#                               base + 'total_efficiency_output2bodydecay-rho1m.txt',
#                               base + 'total_efficiency_output2bodydecay-phi1m.txt',
#                               base + 'total_efficiency_output2bodydecay-omega1m.txt',
#                               base + 'total_efficiency_output2bodydecay-jsi1m.txt',
#                               NPOT, masses, charges, 2, mass))
# mesonDecay_dune = mesonDecayProduction(base + 'total_efficiency_outputdecayPionDARK.txt',
#                                        base + 'total_efficiency_outputdecayEtaDARK.txt',
#                                        base + 'total_efficiency_output2bodydecay-rhoDARK.txt',
#                                        base + 'total_efficiency_output2bodydecay-phiDARK.txt',
#                                        base + 'total_efficiency_output2bodydecay-omegaDARK.txt',
#                                        base + 'total_efficiency_output2bodydecay-jsiDARK.txt',
#                                        NPOT, masses, charges, 3, mass)

# mesonDecay_dune_2layers = mesonDecayProduction(base + 'total_efficiency_outputdecayPionDARK.txt',
#                                        base + 'total_efficiency_outputdecayEtaDARK.txt',
#                                        base + 'total_efficiency_output2bodydecay-rhoDARK.txt',
#                                        base + 'total_efficiency_output2bodydecay-phiDARK.txt',
#                                        base + 'total_efficiency_output2bodydecay-omegaDARK.txt',
#                                        base + 'total_efficiency_output2bodydecay-jsiDARK.txt',
#                                        NPOT, masses, charges, 2, mass)

# mesonDecay_dune = list(mesonDecayProduction(base + 'total_efficiency_outputdecayPion1m.txt',
#                                             base + 'total_efficiency_outputdecayEta1m.txt',
#                                             base + 'total_efficiency_output2bodydecay-rho1m.txt',
#                                             base + 'total_efficiency_output2bodydecay-phi1m.txt',
#                                             base + 'total_efficiency_output2bodydecay-omega1m.txt',
#                                             base + 'total_efficiency_output2bodydecay-jsi1m.txt',
#                                             NPOT, masses, charges, 2, mass))

# plt.contour(mass, charge, total_dune_2layers, levels=[456.41], colors = 'lightgreen', linewidths = 3)
# plt.plot(10, 5, color='lightgreen', label='2 layers, 315 bars/layer, bkg = 50,000')
#
# # plt.contour(mass, charge, total_dune2, levels=[21], colors = 'skyblue')
# # plt.plot(10, 5, color='skyblue', label='2 layers, 0.5 m radius (~ 315 bars/layer), bkg = 100')
#
# plt.contour(mass, charge, total_dune2, levels=[29], colors = 'lightgreen')
# plt.plot(10, 5, color='lightgreen', label='2 layers, 0.5 m radius (~ 315 bars/layer), bkg = 200')
#
# plt.contour(mass, charge, total_dune2, levels=[45], colors = 'gray')
# plt.plot(10, 5, color='gray', label='2 layers, 0.5 m radius (~ 315 bars/layer), bkg = 500')
# #
# plt.contour(mass, charge, total_dune, levels=[10], colors = 'greenyellow')
# plt.plot(10, 5, color='greenyellow', label='2 layers,0.5 m radius (~ 315 bars/layer), bkg = 20, without proton brem')
#
# # plt.contour(mass, charge, total_dune, levels=[21], colors = 'tomato')
# # plt.plot(10, 5, color='tomato', label='2 layers, 0.5 m radius (~ 315 bars/layer), bkg = 100, without proton brem')
#
# plt.contour(mass, charge, total_dune, levels=[29], colors = 'skyblue')
# plt.plot(10, 5, color='skyblue', label='2 layers, 0.5 m radius (~ 315 bars/layer), bkg = 200, without proton brem')
#
# plt.contour(mass, charge, total_dune, levels=[45], colors = 'cyan')
# plt.plot(10, 5, color='cyan', label='2 layers, 0.5 m radius (~ 315 bars/layer), bkg = 500, without proton brem')


# plt.contour(mass, charge, total_dune2, levels=[10], colors = 'greenyellow')
# plt.plot(10, 5, color='greenyellow', label='2 layers, 1 m radius (~ 315 bars/layer), bkg = 20')
#
# plt.contour(mass, charge, total_dune2, levels=[21], colors = 'tomato')
# plt.plot(10, 5, color='tomato', label='2 layers, 0.5 m radius (~ 315 bars/layer), bkg = 100')
#
# plt.contour(mass, charge, total_dune2, levels=[29], colors = 'skyblue')
# plt.plot(10, 5, color='skyblue', label='2 layers, 0.5 m radius (~ 315 bars/layer), bkg = 200')
#
# plt.contour(mass, charge, total_dune2, levels=[45], colors = 'cyan')
# plt.plot(10, 5, color='cyan', label='2 layers, 0.5 m radius (~ 315 bars/layer), bkg = 500')
#
# plt.contour(mass, charge, total_dune2, levels=[63], colors='orange')
# plt.plot(10, 5, color='orange', label='2 layers, 0.5 m radius (~ 315 bars/layer), bkg = 1000')
#
# plt.contour(mass, charge, total_dune2, levels=[201], colors='purple')
# plt.plot(10, 5, color='purple', label='2 layers, 0.5 m radius (~ 315 bars/layer), bkg = 10 000')
#
# plt.contour(mass, charge, total_dune2, levels=[448], colors='darkgreen')
# plt.plot(10, 5, color='darkgreen', label='2 layers, 0.5 m radius (~ 315 bars/layer), bkg = 50 000')

# # ------------------ one-hue combined-constraints envelope (with overlays) ------------------
# def df_to_xy(df, name, x_scale=1.0):
#     if not isinstance(df, pd.DataFrame) or df.shape[1] < 2:
#         print(f"[{name}] not a 2-col DataFrame")
#         return None
#     x = (df[0].to_numpy(dtype=float) * x_scale)
#     y = df[1].to_numpy(dtype=float)
#     m = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
#     x, y = x[m], y[m]
#     if x.size < 2:
#         print(f"[{name}] empty/too small after cleaning")
#         return None
#     o = np.argsort(x)
#     x, y = x[o], y[o]
#     print(f"[{name}] x-range = [{x.min():.3e}, {x.max():.3e}] (n={x.size})")
#     return name, x, y
#
# # pick ONE color for everything (tweak if you like)
# envelope_color = "#E41A1C"   # red similar to your reference figure
# fill_alpha      = 0.15
# parts_alpha     = 0.85
# parts_lw        = 1.6
# envelope_lw     = 2.8
#
# # pretty labels (optional)
# label_map = {
#     "colliders": "LEP/CMS/Colliders",
#     "slac":      "SLAC",
#     "bebc":      "BEBC",
#     "lsnd":      "LSND",
#     "charmii":   "CHARM-II",
#     "argoneut":  "ArgoNeuT",
#     "mq_demo":   "milliQan (demo)",
#     "sensei":    "SENSEI",
#     "mq_run3":   "milliQan (Run 3)",
# }
#
# # Only keep: milliQan (mq_run3), BEBC, SENSEI, LSND, ArgoNeuT, colliders
# curves = []
# for name, df, scale in [
#     ("colliders",  colliders, 1.0),
#     ("bebc",       bebc,      1.0),
#     ("slac",       slac,      1.0),       # <-- added
#     ("lsnd",       lsnd,      1/1000.0),   # MeV -> GeV
#     ("argoneut",   argoneut,  1.0),
#     ("sensei",     sensei,    1/1000.0),   # MeV -> GeV
#     ("mq_run3",    mq_run3,   1.0),
# ]:
#     out = df_to_xy(df, name, x_scale=scale)
#     if out is not None:
#         curves.append(out)
#
# if not curves:
#     print("No curves loaded; check CSV paths.")
# else:
#     # grid built from union of curve ranges (clamped to your axes)
#     global_min = min(x.min() for _, x, _ in curves)
#     global_max = max(x.max() for _, x, _ in curves)
#     m_lo = max(global_min, 1e-3)
#     m_hi = min(global_max, 1e1)
#     print(f"[grid] using [{m_lo:.3e}, {m_hi:.3e}]")
#     mgrid = np.logspace(np.log10(m_lo), np.log10(m_hi), 2000)
#
#     envelope     = np.full_like(mgrid, np.inf, dtype=float)
#     contributed  = np.zeros_like(mgrid, dtype=bool)
#     parts_on_grid = {}  # name -> y(mgrid) with NaNs outside its range
#
#     for name, x, y in curves:
#         lx = np.log10(x)
#         ux, idx = np.unique(lx, return_index=True)
#         if ux.size < 2:
#             print(f"[{name}] skipped (not enough unique x for interp)")
#             continue
#         ly = np.log10(y)[idx]
#
#         yg = np.full_like(mgrid, np.nan, dtype=float)
#         in_rng = (mgrid >= x.min()) & (mgrid <= x.max())
#         if np.any(in_rng):
#             g_lx = np.log10(mgrid[in_rng])
#             g_ly = np.interp(g_lx, ux, ly)
#             yg[in_rng] = 10.0**g_ly
#             envelope[in_rng] = np.minimum(envelope[in_rng], yg[in_rng])
#             contributed[in_rng] = True
#         else:
#             print(f"[{name}] skipped (no overlap with grid)")
#
#         parts_on_grid[name] = yg
#
#     if not np.any(contributed):
#         print("Combined envelope is empty after interpolation (check mass ranges / columns).")
#     else:
#         envelope[~contributed] = np.nan
#         print(f"Envelope stats: min ε = {np.nanmin(envelope):.3e} "
#               f"max ε = {np.nanmax(envelope):.3e} #pts = {np.isfinite(envelope).sum()}")
#
#         # -------------------- plotting (single hue + inline labels) --------------------
#         ax = plt.gca()
#         y_top = ax.get_ylim()[1]
#         mask  = np.isfinite(envelope)
#
#         # combined fill + line
#         ax.fill_between(mgrid[mask], envelope[mask], y_top,
#                         alpha=fill_alpha, color=envelope_color,
#                         zorder=1, label="Existing constraints (combined)  ArgoNeuT + BEBC + LSND + milliQan + SENSEI")
#         # ax.plot(mgrid[mask], envelope[mask],
#         #         color=envelope_color, lw=envelope_lw,
#         #         solid_capstyle="round", zorder=5,
#         #         label="Combined envelope")
#
#         # dashed overlays (same hue) + name printed at right edge
#         dash_bank = [(4,2), (2,2), (6,3,2,3), (1,2), (8,2), (5,3,1,3)]
#         pretty = {
#             "mq_run3":  "milliQan",
#             "bebc":     "BEBC",
#             "sensei":   "SENSEI",
#             "lsnd":     "LSND",
#             "argoneut": "ArgoNeuT",
#             "colliders":"Colliders",
#             "slac":     "SLAC",     # <-- added
#         }
#         color_map = {
#             "colliders": "#1f77b4",  # blue
#             "slac":      "#ff7f0e",  # orange
#             "bebc":      "#2ca02c",  # green
#             "lsnd":      "#9467bd",  # purple
#             "argoneut":  "#8c564b",  # brown
#             "sensei":    "#e377c2",  # pink
#             "mq_run3":   "#17becf",  # teal (milliQan)
#         }
#
#         for i, (name, yg) in enumerate(parts_on_grid.items()):
#             msk = np.isfinite(yg)
#             if not np.any(msk):
#                 continue
#             ax.plot(mgrid[msk], yg[msk],
#                     color=envelope_color, lw=parts_lw,
#                     dashes=dash_bank[i % len(dash_bank)],
#                     alpha=parts_alpha, zorder=3)
#
#             # Label near the rightmost valid point (slight right offset)
#             # --- label positions ---
#             if name == "argoneut":
#                 # beginning of curve, shifted right
#                 x_pos = mgrid[msk][0] * 1.35   # more rightward
#                 y_pos = yg[msk][0] * 0.85
#                 ha = "left"
#
#             elif name == "slac":
#                 # beginning of curve, shifted right
#                 x_pos = mgrid[msk][0] * 0.80
#                 y_pos = yg[msk][0]
#                 ha = "left"
#
#             elif name == "lsnd":
#                 # beginning of curve, shifted right
#                 x_pos = mgrid[msk][0] * 1.9
#                 y_pos = yg[msk][0]* 1.2
#                 ha = "left"
#
#             elif name == "bebc":
#                 # beginning of curve, slight shift
#                 x_pos = mgrid[msk][0] * 1.1
#                 y_pos = yg[msk][0]
#                 ha = "left"
#
#             elif name == "sensei":
#                 # middle of curve, shifted downward
#                 mid_idx = len(mgrid[msk]) // 2
#                 x_pos = mgrid[msk][mid_idx]
#                 y_pos = yg[msk][mid_idx] * 0.85   # push down
#                 ha = "center"
#
#             else:
#                 # default = right end
#                 x_pos = mgrid[msk][-1] * 1.05
#                 y_pos = yg[msk][-1]
#                 ha = "left"
#
#             ax.text(x_pos, y_pos,
#                     pretty.get(name, name),
#                     color=envelope_color,
#                     fontsize=16,            # larger font
#                     va="center", ha=ha,
#                     fontweight="bold")
#
#
#             ax.text(x_pos, y_pos,
#                     pretty.get(name, name),
#                     color=envelope_color,
#                     fontsize=16,            # larger font
#                     va="center", ha=ha,
#                     fontweight="bold")
#
#
#     # Optional legend (you can comment this out if inline labels suffice)
#         handles, labels = ax.get_legend_handles_labels()
#         if "Combined envelope" in labels:
#             i_comb = labels.index("Combined envelope")
#             order = [i_comb] + [i for i in range(len(labels)) if i != i_comb]
#             ax.legend([handles[i] for i in order], [labels[i] for i in order],
#                       frameon=False, loc="best")
#         else:
#             ax.legend(frameon=False, loc="best")
# # ---------------------------------------------------------------------------

# # =========================
# # SECOND PLOT (FORESEE-style)
# # =========================
# from matplotlib import patheffects as pe
#
# def styled_fill(ax, x, y, color, label, ymax=1.0, alpha=0.45, lw=2.2, z=2):
#     """Fill constraint to ymax, then draw a crisp outline with white halo."""
#     ax.fill_between(x, y, ymax, facecolor=color, edgecolor="none", alpha=alpha, zorder=z)
#     ax.plot(x, y, color=color, lw=lw, zorder=z+1,
#             path_effects=[pe.Stroke(linewidth=lw+1.2, foreground="white"), pe.Normal()],
#             label=label)
#
# # Make a new figure with SAME geometry as your first
# fig2 = plt.figure(figsize=(22, 16), dpi=500)
# ax2  = fig2.add_subplot(111)
#
# # Axes config identical to the first plot
# ax2.set_xscale('log')
# ax2.set_yscale('log')
# ax2.set_xlim(0.001, 10)
# ax2.set_ylim(0.000005, 1)
# ax2.set_xlabel('$m_{\\chi}$ [$\\mathrm{GeV}/\\mathrm{c}^2$]', fontsize=35)
# ax2.set_ylabel('$\\epsilon=Q/e$', fontsize=35)
# ax2.tick_params(axis='both', which='both', labelsize=25)
# ax2.set_rasterization_zorder(1)  # heavy fills rasterized, lines stay vector
#
# # ---- Constraints (same sources you already loaded) ----
# # Use the *_small dataframes you already created above.
# styled_fill(ax2, MiniBooNE_small[0].to_numpy(), MiniBooNE_small[1].to_numpy(),
#             'mediumseagreen', 'MiniBooNE', ymax=1.0, alpha=0.45)
# styled_fill(ax2, slac_small[0].to_numpy(), slac_small[1].to_numpy(),
#             'gold', 'SLAC mQ', ymax=1.0, alpha=0.45)
# styled_fill(ax2, bebc_small[0].to_numpy(), bebc_small[1].to_numpy(),
#             'gray', 'BEBC', ymax=1.0, alpha=0.45)
# styled_fill(ax2, lsnd_small[0].to_numpy(), lsnd_small[1].to_numpy(),
#             'skyblue', 'LSND', ymax=1.0, alpha=0.45)
# styled_fill(ax2, charmii_small[0].to_numpy(), charmii_small[1].to_numpy(),
#             'lightgray', 'CHARM II', ymax=1.0, alpha=0.45)
# styled_fill(ax2, argoneut_small[0].to_numpy(), argoneut_small[1].to_numpy(),
#             'cornflowerblue', 'ArgoNeuT', ymax=1.0, alpha=0.45)
# styled_fill(ax2, lep_small[0].to_numpy(), lep_small[1].to_numpy(),
#             'lightcoral', 'LEP', ymax=1.0, alpha=0.45)
# styled_fill(ax2, mq_run3_small[0].to_numpy(), mq_run3_small[1].to_numpy(),
#             '#E6E6FA', 'milliQan', ymax=1.0, alpha=0.45)
#
# # ---- Your sensitivity contour (same level & dataset as first plot) ----
# cs = ax2.contour(mass, charge, total_dune_new, levels=[3.09],
#                  colors=['salmon'], linewidths=3.0, zorder=5)
# # add a halo to the contour lines
# # ---- Your sensitivity contour (same level & dataset as first plot) ----
# cs = ax2.contour(
#     mass, charge, total_dune_new,
#     levels=[3.09],
#     colors=['#e63946'],      # vivid salmon-red tone
#     linewidths=6.5,          # thicker line
#     zorder=10
# )
#
# # optional subtle shadow to make it stand out, no halo
# if hasattr(cs, "collections"):
#     for c in cs.collections:
#         c.set_path_effects([
#             pe.SimpleLineShadow(offset=(2, -2), alpha=0.35),
#             pe.Normal()
#         ])
# else:
#     print("Warning: contour has no .collections attribute; skipping shadow effect.")
#
# # Legend entry (matches the contour color)
# ax2.plot([], [], color='#e63946', lw=6.5, label='DarkQuest Sensitivity (3 events)')
# ax2.legend(loc='upper left', fontsize=18.5, frameon=False)
#
# fig2.savefig('limit-plot-DarkQuest-FORESEEstyle-' + str(NPOT) + '.pdf')
# print("FORESEE-style plot drawn")
