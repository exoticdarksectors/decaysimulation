import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import glob, os

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.family'] = 'serif'
plt.rcParams['xtick.labelsize'] = 16
plt.rcParams['ytick.labelsize'] = 16

def mesonDecayProduction(fileNamePi, fileNameEta, fileNameRho, fileNamePhi, fileNameOmega, fileNameJpsi, NPOT, fileMass):
    # Constants
    # EM constant
    alpha = 1.0 / 137
    # mass in GeV
    m_e = 0.00051
    m_pi = 0.135
    m_eta = 0.548
    m_rho = 0.775
    m_omega = 0.782
    m_phi = 1.019
    m_jpsi = 3.1
    m_upsilon = 9.46
    # meson / NPOT obtained from PYTHIA
    c_pi = 7.5
    c_eta = 0.85
    c_rho = 1.0
    c_omega = 1.0
    c_phi = 4.1e-02
    c_jpsi = 8.3e-5
    c_upsilon = 5.5e-9
    # m -> e+e- branching ratio
    branch_pi = 0.98
    branch_eta = 0.39
    branch_rho = 4.72e-5
    branch_omega = 7.28e-5
    branch_phi = 2.95e-4
    branch_jpsi = 0.05971
    branch_upsilon = 0.0238

    # Import data from efficiency files
    ageo_datapi = np.loadtxt(fileNamePi)
    ageo_dataeta = np.loadtxt(fileNameEta)
    ageo_datarho = np.loadtxt(fileNameRho)
    ageo_dataphi = np.loadtxt(fileNamePhi)
    ageo_dataomega = np.loadtxt(fileNameOmega)
    ageo_datajpsi = np.loadtxt(fileNameJpsi)

    ageo_rho = ageo_datarho[:, 1]
    ageo_omega = ageo_dataomega[:, 1]
    ageo_phi = ageo_dataphi[:, 1]
    ageo_pi = ageo_datapi[:, 1]
    ageo_jpsi = ageo_datajpsi[:, 1]
    ageo_eta = ageo_dataeta[:, 1]

    # Separate efficiency and mass
    mass_rho = ageo_datarho[:,0]
    mass_omega = ageo_dataomega[:,0]
    mass_phi = ageo_dataphi[:,0]
    mass_pi = ageo_datapi[:,0]
    mass_eta = ageo_dataeta[:,0]
    mass_jpsi = ageo_datajpsi[:,0]
    mass_upsilon = fileMass[:]

    ageo_upsilon = np.full(mass_upsilon.size, 0.011338) #taken from J/psi efficiency at 1e-3 GeV

    # define phase space integrals
    def I2(x, y):
        return ((1 + 2 * x) * np.sqrt(1 - 4 * x)) / ((1 + 2 * y) * np.sqrt(1 - 4 * y))
    def I3_integrand(z, x):
        return 2 / (3 * 3.14) * np.sqrt(1 - 4 * x / z) * ((1 - z) ** 3) * (2 * x + z) / (z ** 2)
    def I3(x):
        return quad(I3_integrand, 4 * x, 1, args=x)[0]  # integrate
    v_I3 = np.vectorize(I3) # vectorize

    # define productions
    pi = np.zeros(np.shape(mass_pi))
    eta = np.zeros(np.shape(mass_eta))
    phi = np.zeros(np.shape(mass_phi))
    omega = np.zeros(np.shape(mass_omega))
    jpsi = np.zeros(np.shape(mass_jpsi))
    rho = np.zeros(np.shape(mass_rho))
    upsilon = np.zeros(np.shape(mass_upsilon))

    # Dalitz decay
    pi[mass_pi < m_pi/2] = NPOT * ageo_pi[mass_pi < m_pi/2] * 2 * c_pi * branch_pi * alpha * v_I3( mass_pi[mass_pi < m_pi/2] ** 2 / m_pi ** 2)
    eta[mass_eta < m_eta/2] = NPOT * ageo_eta[mass_eta < m_eta/2] * 2 * c_eta * branch_eta * alpha * v_I3( mass_eta[mass_eta < m_eta/2] ** 2 / m_eta ** 2)

    # Direct decay
    jpsi[mass_jpsi < m_jpsi/2] = NPOT * ageo_jpsi[mass_jpsi < m_jpsi/2] * 2  * c_jpsi * branch_jpsi * I2( mass_jpsi[mass_jpsi < m_jpsi/2] ** 2 / m_jpsi ** 2, m_e ** 2 / m_jpsi ** 2  )
    rho[mass_rho < m_rho/2] = NPOT * ageo_rho[mass_rho < m_rho/2] * 2  * c_rho * branch_rho * I2( mass_rho[mass_rho < m_rho/2] ** 2 / m_rho ** 2, m_e ** 2 / m_rho ** 2  )
    omega[mass_omega < m_omega/2] = NPOT * ageo_omega[mass_omega < m_omega/2] * 2  * c_omega * branch_omega * I2( mass_omega[mass_omega < m_omega/2] ** 2 / m_omega ** 2, m_e ** 2 / m_omega ** 2  )
    phi[mass_phi < m_phi/2] = NPOT * ageo_phi[mass_phi < m_phi/2] * 2  * c_phi * branch_phi * I2( mass_phi[mass_phi < m_phi/2] ** 2 / m_phi ** 2, m_e ** 2 / m_phi ** 2  )
    upsilon[mass_upsilon < m_upsilon/2] = NPOT * ageo_upsilon[mass_upsilon < m_upsilon/2] * 2  * c_upsilon * branch_upsilon * I2( mass_upsilon[mass_upsilon < m_upsilon/2] ** 2 / m_upsilon ** 2, m_e ** 2 / m_upsilon ** 2  )

    return pi, eta, jpsi, upsilon, rho, omega, phi

def _effective_lumi_pb(N_POT=0, rho=0, L=0, A=0):
    """
    Compute effective luminosity in pb^-1: L_eff = N_POT * (rho * L / A) * N_A
    """
    NA = 6.02214076e23  # /mol
    nuclei_per_cm2 = (rho * L / A) * NA
    L_eff_cm2 = N_POT * nuclei_per_cm2
    return L_eff_cm2 / 1e36 # convert 1/cm^2 units to pb^-1

def brem_single_hit_or_shared_v2(
        directory,
        pattern="Brem_*.txt",
        det_angle=0,   # radians
        epsilon=1.0,
        lambda_idx=1, # 0:1.0 GeV, 1:1.5 GeV, 2:2.0 GeV
        N_POT=0, rho=0, L=0, A=0
):
    """
    Reads single_hit_or_shared files with the following this structure:
      log10(theta), log10(p/GeV), sigma_L1, sigma_L1p5, sigma_L2   [pb/bin]
    Returns: masses, # of track
    """
    files = sorted(glob.glob(os.path.join(directory, pattern)))
    if not files:
        return np.array([]), np.array([])
    L_eff = _effective_lumi_pb(N_POT, rho, L, A)
    log10_theta_cut = np.log10(det_angle)

    masses, tracks = [], []
    for f in files:
        try:
            mchi = float(os.path.basename(f).split('_')[-1].replace('.txt', ''))
        except Exception:
            continue

        arr = np.loadtxt(f, comments="#")
        if arr.ndim == 1: arr = arr[None, :]
        if arr.shape[1] < 5:  # need 3 sigma columns
            continue

        log10theta = arr[:, 0]
        sigmas = arr[:, 2:5]  # [nbins, 3]
        if lambda_idx not in (0, 1, 2):
            raise ValueError("lambda_idx must be 0, 1, or 2")

        # acceptance on the listed particle
        mask = (log10theta < log10_theta_cut)
        sigma_sum_pb = sigmas[mask, lambda_idx].sum()

        N_tracks = (epsilon**2) * L_eff * sigma_sum_pb

        masses.append(mchi)
        tracks.append(N_tracks)

    masses = np.asarray(masses); tracks = np.asarray(tracks)
    idx = np.argsort(masses)
    return masses[idx], tracks[idx]

def dyProductionTOT(fileNamecross, filenamedy, NPOT, totalCrossSection):
    # import data from efficiency files
    ageo_datacross = np.loadtxt(fileNamecross)
    ageo_datady = np.loadtxt(filenamedy)

    # mass = ageo_datady[:, 0]
    cross_dy = ageo_datacross[:] * 1e-12 # pico barn
    ageo_dy = ageo_datady[:,1]

    return NPOT * cross_dy / totalCrossSection * ageo_dy

def dyProduction(fileNamecross, filenamedy, NPOT, totalCrossSection, cutoff=1.0):
    # Load data
    ageo_datacross = np.loadtxt(fileNamecross)
    ageo_datady    = np.loadtxt(filenamedy)

    # Extract columns
    mass    = ageo_datady[:, 0]
    ageo_dy = ageo_datady[:, 1]

    # cross_dy
    if ageo_datacross.ndim == 1:
        cross_dy = ageo_datacross
    else:
        cross_dy = ageo_datacross[:, 0]

    cross_dy = cross_dy * 1e-12

    n = min(mass.shape[0], cross_dy.shape[0], ageo_dy.shape[0])
    mass    = mass[:n]
    ageo_dy = ageo_dy[:n]
    cross_dy = cross_dy[:n]

    # Compute production
    dy_prod = NPOT * cross_dy / totalCrossSection * ageo_dy
    dy_prod[mass < cutoff] = np.nan
    return dy_prod

if __name__ == '__main__':
    base = '/Users/leobailloeul/Documents/coding/Archive-dec-9-2025/data-backup/'
    mass_ship = np.loadtxt(base + 'mship_values_copy.txt')
    NPOT_SHiP = 2e20

    mesonDecay_ship = list(mesonDecayProduction(base + 'total_efficiency_outputdecayPion_SHIP-1m-small.txt',
                                                base + 'total_efficiency_outputdecayEta_SHIP-1m-small.txt',
                                                base + 'total_efficiency_output2bodydecay-rho-SHIP.txt',
                                                base + 'total_efficiency_output2bodydecay-phi-SHIP.txt',
                                                base + 'total_efficiency_output2bodydecay-omega-SHIP.txt',
                                                base + 'total_efficiency_output2body_decay-jsi-SHiP.txt',
                                                NPOT_SHiP, mass_ship))

    # DY production
    dy_ship = dyProduction(base+'dy_cross_ship_nCTEQ15_iron.txt', base+'ageo_ship_nCTEQ_iron.txt', NPOT_SHiP, 13e-3)
    dy_shiptot = dyProductionTOT(base+'dy_cross_ship_nCTEQ15_iron.txt', base+'ageo_ship_nCTEQ_iron.txt', NPOT_SHiP, 13e-3)

    # add paddings
    for i in range(len(mesonDecay_ship)):
        difference = mass_ship.size - mesonDecay_ship[i].size
        mesonDecay_ship[i] = np.pad(mesonDecay_ship[i], (0, difference), 'constant',  constant_values=0)

    difference1 = mass_ship.size - dy_shiptot.size
    dy_shiptot = np.pad(dy_shiptot, (0, difference1), 'constant', constant_values=0)

    total_ship = dy_shiptot.copy()
    for i in range(len(mesonDecay_ship)):
        total_ship += mesonDecay_ship[i]

    brem_dir_shared = "/Users/leobailloeul/Documents/coding/Archive-dec-9-2025/SHiP-brem-backup"
    det_angle = np.arctan(0.5/100.0)
    epsilon = 1.0

    # Sigma p = 1.0, 1.5, 2.0 GeV (lower, central, upper)
    m_low, y_low = brem_single_hit_or_shared_v2(
        brem_dir_shared, det_angle=det_angle, epsilon=epsilon,
        lambda_idx=0, N_POT=NPOT_SHiP, rho=10.2, L=15.27, A=95.95
    )
    m_cent, y_cent = brem_single_hit_or_shared_v2(
        brem_dir_shared, det_angle=det_angle, epsilon=epsilon,
        lambda_idx=1, N_POT=NPOT_SHiP, rho=10.2, L=15.27, A=95.95
    )
    m_high, y_high = brem_single_hit_or_shared_v2(
        brem_dir_shared, det_angle=det_angle, epsilon=epsilon,
        lambda_idx=2, N_POT=NPOT_SHiP, rho=10.2, L=15.27, A=95.95
    )

    # put brem on a common internal mass grid
    m_brem_grid = np.unique(np.concatenate([m_low, m_cent, m_high]))

    y_low_g = np.interp(m_brem_grid, m_low, y_low, left=0.0, right=0.0)
    y_cent_g = np.interp(m_brem_grid, m_cent, y_cent, left=0.0, right=0.0)
    y_high_g = np.interp(m_brem_grid, m_high, y_high, left=0.0, right=0.0)

    # interpolate all three onto the main mass grid used in the plot (mass_ship)
    y_low_ship = np.interp(mass_ship, m_brem_grid, y_low_g, left=0.0, right=0.0)
    y_cent_ship = np.interp(mass_ship, m_brem_grid, y_cent_g, left=0.0, right=0.0)
    y_high_ship = np.interp(mass_ship, m_brem_grid, y_high_g, left=0.0, right=0.0)

    # add ONLY the central Λp = 1.5 GeV brem to the total yield
    total_ship += y_cent_ship
    y_cent_plot = y_cent_ship.copy()

    # plotting
    fig, ax = plt.subplots(figsize=(8.09, 5))
    # axes + scales
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.margins(x=0, y=0)
    ax.set_xlabel(r'$m_{\chi}\, [\mathrm{GeV}/c^{2}]$', fontsize=18)
    ax.set_ylabel(r'$N_{\chi}/\epsilon^{2}$', fontsize=18)

    # ticks with powers of 10
    x_ticks = [10**i for i in range(-2, 2)]
    x_tick_labels = [r'$10^{-2}$', r'$10^{-1}$', r'$1$', r'$10^{1}$']
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_tick_labels)

    y_ticks = [10**i for i in range(0, 20)]
    y_tick_labels = [r'$1$'] + [ (rf'$10^{{{i}}}$' if i % 2 == 0 else ' ') for i in range(1, 20) ]
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_tick_labels)
    ax.set_ylim(1e6, 1e20)
    ax.set_xlim(1e-3, 1e1)

    # lines (use your data; order/labels corrected)
    lw = 3.0
    ax.plot(mass_ship, mesonDecay_ship[0], lw=lw, label=r'$\pi^{0}\!\to\!\gamma\chi\bar{\chi}$')
    ax.plot(mass_ship, mesonDecay_ship[1], lw=lw, label=r'$\eta\!\to\!\gamma\chi\bar{\chi}$')
    ax.plot(mass_ship, mesonDecay_ship[2], lw=lw, label=r'$J/\psi\!\to\!\chi\bar{\chi}$')
    ax.plot(mass_ship, mesonDecay_ship[3], lw=lw, label=r'$\Upsilon\!\to\!\chi\bar{\chi}$')
    ax.plot(mass_ship, mesonDecay_ship[4], lw=lw, label=r'$\rho\!\to\!\chi\bar{\chi}$')
    ax.plot(mass_ship, mesonDecay_ship[5], lw=lw, label=r'$\omega\!\to\!\chi\bar{\chi}$')
    ax.plot(mass_ship, mesonDecay_ship[6], lw=lw, label=r'$\phi\!\to\!\chi\bar{\chi}$')
    ax.plot(mass_ship, dy_ship, lw=3, label = r'$q\overline{q}\to\gamma^{*}\to\chi\overline{\chi}$')
    ax.plot(mass_ship, y_cent_plot, color='#bcbd22', lw=3, ls='--',
            label=r'$p$ brem ($\Lambda_p = 1.5~\mathrm{GeV}$)', zorder=2)
    ax.plot(mass_ship, total_ship, color='black', lw=3, label='Total')

    y_low_plot  = y_low_ship.copy()
    y_high_plot = y_high_ship.copy()
    y_low_plot[y_low_plot <= 0]   = np.nan
    y_high_plot[y_high_plot <= 0] = np.nan

    # light yellow band between Λp = 1.0 and 2.0 GeV
    ax.fill_between(mass_ship, y_low_plot, y_high_plot, color='#bcbd22', alpha=0.5, linewidth=0.0, zorder=1)

    # legend placement and dimensions
    leg = ax.legend(loc='lower left', ncol=2, frameon=True, fancybox=True, framealpha=0.85,
                    fontsize=12.5, handlelength=2.2, markerscale=1.1)
    # for line in leg.get_lines():
    #     line.set_linewidth(2.2)

    textbox_props = dict(boxstyle='round', facecolor='white', edgecolor='none', alpha=0.8)
    ax.text(0.47, 0.93,r'$5$ Years at SHiP ($2\times10^{20}$ POT)' '\n'
                       r'$100$ m, $0.5$ m radius cylindrical detector',
            transform=ax.transAxes, ha='left', va='center', bbox=textbox_props, fontsize=12.5)
    
    # save plot
    fig.tight_layout()
    fig.savefig('MCP_production_SHiP.pdf', bbox_inches='tight', pad_inches=0.1)
