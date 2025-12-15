import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from matplotlib.pyplot import figure
import pandas as pd
import os, glob

plt.rcParams['mathtext.fontset'] = 'cm'   # use Computer Modern
plt.rcParams['font.family'] = 'serif'

# geometry + detection
NPOT = 1e20
NPOT_SHiP = 2e20
A_DET = 3
N_GAMMA = 2.5e5
R_DET = 0.5
L_BASE = 40.0
L_BASE_SHiP = 100.0
THETA_CUT = np.arctan(R_DET / L_BASE)
THETA_CUT_SHiP = np.arctan(R_DET / L_BASE_SHiP)

def _effective_lumi_pb(N_POT=1.0e20, rho=7.87, L=16.8, A=55.845):
    """
    Compute effective luminosity in pb^-1: L_eff = N_POT * (rho * L / A) * N_A
    """
    NA = 6.02214076e23  # /mol
    nuclei_per_cm2 = (rho * L / A) * NA
    L_eff_cm2 = N_POT * nuclei_per_cm2
    return L_eff_cm2 / 1e36

def brem_v2_to_grid(m_grid_masses, m1d, y1d_baseline, charges, a=3, N_gamma=2.5e5):
    """
    Interpolate 'baseline' onto mass grid, then broadcast with detection factor:
    (1 - exp(-Nγ ε^2))^a * ε^2
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
        pattern="Brem_*.txt",
        det_angle=np.arctan(0.5/40.0),
        lambda_idx=1,                 # 0→Λp=1.0, 1→1.5, 2→2.0
        N_POT=1.0e20, rho=7.87, L=16.8, A=55.845,
):
    """
    Reads single_hit_or_shared files with the following this structure:
      log10(theta), log10(p/GeV), sigma_L1, sigma_L1p5, sigma_L2   [pb/bin]
    Returns: masses, # of track
    Additionally two chi in same bar already accounted
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

def mesonProduction_2d(fileNamePi, fileNameEta, fileNameRho, fileNamePhi, fileNameOmega, fileNameJpsi, fileMass, NPOT, charges, a, N_gamma=2.5e5, SHIP=False):
    # constants
    alpha = 1.0/137.0
    m_e   = 0.00051
    m_pi, m_eta, m_rho, m_omega, m_phi, m_jpsi, m_ups = 0.135, 0.548, 0.775, 0.782, 1.019, 3.1, 9.46
    # meson / NPOT (PYTHIA)
    if SHIP is True:
        c_pi, c_eta, c_rho, c_omega, c_phi, c_jpsi, c_ups = 7.5, 0.85, 1.0, 1.0, 4.1e-2, 8.3e-5, 5.5e-9
    else:
        c_pi, c_eta, c_rho, c_omega, c_phi, c_jpsi, c_ups = 4.7, 0.53, 0.61, 0.61, 2.2e-2, 4.0e-5, 2.5e-9
    # BR to e+e-
    br_pi, br_eta = 0.98, 0.39
    br_rho, br_omega, br_phi = 4.72e-5, 7.28e-5, 2.95e-4
    br_jpsi, br_ups = 0.05971, 0.0238

    # load + align acceptances to fileMass (shape (nx,))
    ageo_pi  = _load_ageo_aligned(fileNamePi,   fileMass)
    ageo_eta = _load_ageo_aligned(fileNameEta,  fileMass)
    ageo_rho = _load_ageo_aligned(fileNameRho,  fileMass)
    ageo_phi = _load_ageo_aligned(fileNamePhi,  fileMass)
    ageo_om  = _load_ageo_aligned(fileNameOmega,fileMass)
    ageo_j   = _load_ageo_aligned(fileNameJpsi, fileMass)
    if SHIP is True:
        ageo_ups = np.full(fileMass.size, 0.011338)
    else:
        ageo_ups = np.full(fileMass.size, 0.011357)

    # define phase space integrals
    def I2(x, y):
        return ((1 + 2*x) * np.sqrt(1 - 4*x)) / ((1 + 2*y) * np.sqrt(1 - 4*y))
    def _I3_integrand(z, x):
        return 2/(3*3.14) * np.sqrt(1 - 4*x/z) * (1 - z)**3 * (2*x + z) / (z**2)
    def I3(x):
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
                    N_gamma=2.5e5, totalCrossSection=13e-3):
    """
    Produces a 2D array of DY yields vs (epsilon, mass) using a shared fileMass grid
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


# masses SHIP vs Dark
mass = np.loadtxt('/Users/leobailloeul/Documents/coding/Archive-dec-9-2025/data-backup/mship_values.txt')
massSHiP = np.loadtxt('/Users/leobailloeul/Documents/coding/Archive-dec-9-2025/data-backup/mship_values_copy.txt')
mass = mass[:-1]
charge = np.logspace(-5, 1, 500)
base = '/Users/leobailloeul/Documents/coding/Archive-dec-9-2025/data-backup/'

BREM_SHARED_DIR = "/Users/leobailloeul/Documents/coding/Archive-dec-9-2025/DarkQuest-brem-backup"
BREM_SHARED_DIR_SHiP = "/Users/leobailloeul/Documents/coding/Archive-dec-9-2025/SHiP-brem-backup"

mesonDecay_dune = mesonProduction_2d(
    base + 'total_efficiency_outputdecayPionDARK.txt',
    base + 'total_efficiency_outputdecayEtaDARK.txt',
    base + 'total_efficiency_output2bodydecay-rho.txt',
    base + 'total_efficiency_output2bodydecay-phi.txt',
    base + 'total_efficiency_output2bodydecay-omega.txt',
    base + 'total_efficiency_output2body_decay-jsi.txt',
    fileMass=mass, NPOT=NPOT, charges=charge, a=3, SHIP=False
)
mesonDecay_SHiP = mesonProduction_2d(
    base + 'total_efficiency_outputdecayPion_SHIP-1m-small.txt',
    base + 'total_efficiency_outputdecayEta_SHIP-1m-small.txt',
    base + 'total_efficiency_output2bodydecay-rho-SHIP.txt',
    base + 'total_efficiency_output2bodydecay-phi-SHIP.txt',
    base + 'total_efficiency_output2bodydecay-omega-SHIP.txt',
    base + 'total_efficiency_output2body_decay-jsi-SHiP.txt',
    fileMass=massSHiP, NPOT=NPOT_SHiP, charges=charge, a=3, SHIP=True
)

dy_dune_new = dyProduction_2d(
    fileNamecross=base+'dy_cross_darkquest_nCTEQ15_iron.txt',
    filenamedy=base+'ageo_darkquest_nCTEQ15_iron.txt',
    fileMass=mass, NPOT=NPOT, charges=charge, a=3,
    totalCrossSection=13e-3,
)

dy_SHiP_new = dyProduction_2d(
    fileNamecross=base+'dy_cross_ship_nCTEQ15_iron.txt',
    filenamedy=base+'ageo_ship_nCTEQ_iron.txt',
    fileMass=massSHiP, NPOT=NPOT_SHiP, charges=charge, a=3,
    totalCrossSection=13e-3,
)

# prepping sensitivity for proton brem comparison
total_dune = mesonDecay_dune.copy()
total_dune += dy_dune_new

total_dune_new = mesonDecay_dune.copy()
total_dune_new += dy_dune_new

total_SHiP = mesonDecay_SHiP.copy()
total_SHiP += dy_SHiP_new

total_SHiP_new = mesonDecay_SHiP.copy()
total_SHiP_new += dy_SHiP_new

m_shared, y_shared0 = _brem_v2_folder_baseline_single_or_shared(
    directory=BREM_SHARED_DIR,
    det_angle=THETA_CUT, lambda_idx=1,  # Λp = 1.5 GeV
    N_POT=NPOT, rho=7.87, L=16.8, A=55.845,
)

# Maps pbrem with detection factor
brem_2d = brem_v2_to_grid(
    m_grid_masses=mass,
    m1d=m_shared,
    y1d_baseline=y_shared0,
    charges=charge,
    a=A_DET,
    N_gamma=N_GAMMA,
)

m_shared_SHiP, y_shared0_SHiP = _brem_v2_folder_baseline_single_or_shared(
    directory=BREM_SHARED_DIR_SHiP,
    det_angle=THETA_CUT_SHiP, lambda_idx=1,  # Λp = 1.5 GeV
    N_POT=NPOT_SHiP, rho=10.2, L=15.27, A=95.95,
)
brem_2d_SHiP = brem_v2_to_grid(
    m_grid_masses=massSHiP,          # your sensitivity mass grid (trimmed to L above)
    m1d=m_shared_SHiP,                # single-or-shared masses
    y1d_baseline=y_shared0_SHiP,      # single-or-shared baseline events (no ε, no Poisson)
    charges=charge,
    a=A_DET,
    N_gamma=N_GAMMA,
)

total_dune_new = total_dune_new + brem_2d
total_SHiP_new = total_SHiP_new + brem_2d_SHiP

print("Generation completed")

figure(figsize=(22, 16), dpi=500)
plt.tick_params(axis='both', which='both', labelsize=25)

cs1 = plt.contour(mass, charge, total_dune, levels=[3.09], colors = '#C94F65', linewidths = 3)
plt.plot(10, 5, color='#C94F65', linewidth=3.1, linestyle='--', label='$10^{20}$ POT at SpinQuest')

plt.contour(mass, charge, total_dune_new,levels=[3.09], colors='#C94F65', linewidths=3)
plt.plot(10, 5, color='#C94F65',  linewidth=3.1, label='$10^{20}$ POT at SpinQuest, including proton brem')

for c in cs1.collections:
    c.set_linestyle('--')

cs2 = plt.contour(massSHiP, charge, total_SHiP, levels=[3.09], colors = '#0487FF', linewidths = 3)
plt.plot(10, 5, color='#0487FF', linewidth=3.1, linestyle='--', label='5 years at SHiP ($2\\times10^{20}$ POT)')

plt.contour(massSHiP, charge, total_SHiP_new,levels=[3.09], colors='#0487FF', linewidths=3)
for c in cs2.collections:
    c.set_linestyle('--')

plt.plot(10, 5, color='#0487FF',  linewidth=3.1, label='5 years at SHiP ($2\\times10^{20}$ POT), including proton brem')

def make_envelope(x, y, m_grid, y_top=2.0):
    x = np.asarray(x)
    y = np.asarray(y)

    order = np.argsort(x)
    x_sorted = x[order]
    y_sorted = y[order]

    y_interp = np.interp(m_grid, x_sorted, y_sorted, left=y_top, right=y_top)
    return y_interp

slac_small = pd.read_csv('../experiment-contours-small/SLACmQ.csv',header=None)
MiniBooNE_small = pd.read_csv('../experiment-contours-small/MiniBooNE.csv',header=None)
bebc_small = pd.read_csv('../experiment-contours-small/BEBC.csv',header=None)
charmii_small = pd.read_csv('../experiment-contours-small/CHARMII.csv',header=None)
argoneut_small = pd.read_csv('../experiment-contours-small/ArgoNeuT.csv',header=None)
lsnd_small = pd.read_csv('../experiment-contours-small/LSND.csv', header=None)
sensei_small = pd.read_csv('../experiment-contours-small/SENSEI_POLISHED.csv', header=None)
mq_run3_small = pd.read_csv('../experiment-contours-small/milliQanRun3Fix.csv', header=None)
lep_small = pd.read_csv('../experiment-contours-small/LEP_edited.csv', header=None)

# Collect all mass / epsilon arrays (make sure units are consistent!)
all_mass = [
    MiniBooNE_small[0],
    slac_small[0],
    bebc_small[0],
    lsnd_small[0],
    charmii_small[0],
    argoneut_small[0],
    # lep_small[0],
    sensei_small[0] / 1000.,  # convert SENSEI to GeV if in MeV
    mq_run3_small[0],
]

all_eps = [
    MiniBooNE_small[1],
    slac_small[1],
    bebc_small[1],
    lsnd_small[1],
    charmii_small[1],
    argoneut_small[1],
    # lep_small[1],
    sensei_small[1],
    mq_run3_small[1],
]

m_grid = np.logspace(np.log10(0.0009916136673313), np.log10(11.3), 800)
y_top = 2.0

# curves_on_grid will store each experiment's smoothed curve on m_grid
curves_on_grid = []
env = np.full_like(m_grid, y_top, dtype=float)

for mx, eps in zip(all_mass, all_eps):
    mx = np.asarray(mx)
    eps = np.asarray(eps)

    # sort in x (mass)
    order = np.argsort(mx)
    mx_sorted = mx[order]
    eps_sorted = eps[order]
    valid = (mx_sorted > 0) & (eps_sorted > 0)
    mx_log = np.log10(mx_sorted[valid])
    eps_log = np.log10(eps_sorted[valid])

    eps_interp_log = np.interp(
        np.log10(m_grid),
        mx_log,
        eps_log,
        left=np.log10(y_top),
        right=np.log10(y_top)
    )

    eps_interp = 10 ** eps_interp_log

    curves_on_grid.append(eps_interp)
    env = np.minimum(env, eps_interp)

plt.fill_between(
    m_grid, env, y_top,
    color='lightgray', alpha=0.7,
    label='Existing constraints'
)
line_kwargs = dict(color='dimgray', lw=2.0, alpha=0.9)  # thicker lines

constraint_labels = [
    'MiniBooNE',
    'SLAC mQ',
    'BEBC',
    'LSND',
    'CHARM II',
    'ArgoNeuT',
    # 'LEP',
    'SENSEI',
    'milliQan',
]

ax = plt.gca()

for i, y_curve in enumerate(curves_on_grid):
    y_plot = np.where(y_curve < y_top, y_curve, np.nan)
    # Plot the line
    ax.plot(
        m_grid, y_plot,
        **line_kwargs,
    )

    label = constraint_labels[i]
    mask = np.isfinite(y_plot)
    if not np.any(mask):
        continue

    x_last = m_grid[mask][-1]
    y_last = y_plot[mask][-1]
    rot = 0

    x_text = x_last
    y_text = y_last
    if label == 'SENSEI':
        x_text = x_last * 0.0465
        y_text = y_last * 0.0145
        rot = 12.5
    if label == 'SLAC mQ':
        x_text = x_last * 0.038
        y_text = y_last * 0.000595
        rot = 17
    if label == 'LSND':
        x_text = x_last * 0.016
        y_text = y_last * 0.00049
        rot = 9
    if label == 'LEP':
        x_text = x_last * 0.61
        y_text = y_last * 1.15
    if label == 'BEBC':
        x_text = x_last * 0.210
        y_text = y_last * 0.0243
        rot = 13
    if label == 'CHARM II':
        x_text = x_last * 0.002
        y_text = y_last * 0.0036
        rot = 6
    if label == 'ArgoNeuT':
        x_text = x_last * 0.003
        y_text = y_last * 0.0105
        rot = 2
    if label == 'MiniBooNE':
        x_text = x_last * 0.009
        y_text = y_last * 0.00175
        rot = 7
    if label == 'milliQan':
        x_text = x_last * 0.270
        y_text = y_last * 0.3
        rot = 45

    ax.text(
        x_text, y_text,
        label,
        fontsize=22,
        color='dimgray',
        ha='left',
        va='center',
        rotation=rot,
        rotation_mode='anchor',
    )

plt.xscale('log')
plt.yscale('log')
plt.xlabel('$m_{\chi}$ [$\mathrm{GeV}/\mathrm{c}^2$]', fontsize = 35)
plt.xticks(fontsize = 35)
plt.ylabel('$\epsilon=Q/e$', fontsize = 35)
plt.yticks(fontsize = 35)
plt.xlim(0.001, 11.2)
plt.ylim(0.000005, 1)
plt.legend(loc='upper left', fontsize=24)

plt.savefig('limit_plot_DarkQuest_SHiP.pdf')
print("limit plot drawn")
