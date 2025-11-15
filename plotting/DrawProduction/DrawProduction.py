import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import ScalarFormatter
from scipy.integrate import quad
# define matplotlib pgf output settings
# matplotlib.use("pgf")
import os, glob

# Update Matplotlib parameters to not use LaTeX for rendering
matplotlib.rcParams.update({
    'font.family': 'serif',  # Use a serif font
    'text.usetex': False,    # Disable LaTeX rendering
})

# define matplotlib parameters
# SMALL_SIZE = 20
# MEDIUM_SIZE = 20
# BIGGER_SIZE = 20
#
# plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
# plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
# plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
# plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
# plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
# plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
# plt.rc('figure', titlesize=BIGGER_SIZE + 5)  # fontsize of the figure title


def mesonDecayProduction(fileNamePi, fileNameEta, fileNameRho, fileNamePhi, fileNameOmega, fileNameJpsi, NPOT, fileMass):
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
    c_pi = 4.678401
    c_eta = 0.527576
    c_rho = 0.24
    c_omega = 0.24
    c_phi = 4.9e-03
    c_jpsi = 4.087e-5
    c_upsilon = 2.5e-9

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

    # Separate efficiency and mass
    mass_rho = ageo_datarho[:,0]
    mass_omega = ageo_dataomega[:,0]
    mass_phi = ageo_dataphi[:,0]
    mass_pi = ageo_datapi[:,0]
    mass_eta = ageo_dataeta[:,0]
    mass_jpsi = ageo_datajpsi[:,0]
    mass_upsilon = fileMass[:]

    ageo_rho = ageo_datarho[:,1]
    ageo_omega = ageo_dataomega[:,1]
    ageo_phi = ageo_dataphi[:,1]
    ageo_pi = ageo_datapi[:,1]
    ageo_jpsi = ageo_datajpsi[:,1]
    # special case for abnormal file structures
    if fileNameEta == '/Users/leobailloeul/Documents/coding/decaysimulation/plotting/sensitivity-plot/ageo_darkquest_0.5m.txt':
        ageo_eta = ageo_dataeta[:,2]
    else:
        ageo_eta = ageo_dataeta[:,1]
    if fileNameJpsi == '/Users/leobailloeul/Documents/coding/decaysimulation/plotting/sensitivity-plot/total_efficiency_output2bodydecay-jsi.txt':
        ageo_upsilon = np.full(mass_upsilon.size, 0.010276)
    else:
        ageo_upsilon = np.full(mass_upsilon.size, 0.010276)

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

    # masking to achieve kinematical validity
    # Dalitz decay
    #    pi[mass] = NPOT * ageo_pi[mass] * 2 * c_pi * branch_pi * alpha * v_I3( mass[mass] ** 2 / mass[mass] ** 2)
    #    eta[mass] = NPOT * ageo_eta[mass] * 2 * c_eta * branch_eta * alpha * v_I3( mass[mass] ** 2 / mass[mass] ** 2)
    #    # Direct decay
    #    jpsi[mass] = NPOT * ageo_jpsi[mass] * 2  * c_jpsi * branch_jpsi * I2( mass[mass] ** 2 / m_jpsi ** 2, m_e ** 2 / m_jpsi ** 2  )
    #    upsilon[mass] = NPOT * ageo_upsilon[mass] * 2 * c_upsilon * branch_upsilon * I2( mass[mass] ** 2 / m_upsilon ** 2, m_e ** 2 / m_upsilon ** 2  )

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

def _effective_lumi_pb(N_POT=1.0e20, rho=7.87, L=16.8, A=55.845):
    """
    Compute effective luminosity in pb^-1:
      L_eff = N_POT * (rho * L / A) * N_A   [1/cm^2]  -> convert to pb^-1
    """
    NA = 6.02214076e23  # /mol
    nuclei_per_cm2 = (rho * L / A) * NA           # 1/cm^2
    L_eff_cm2 = N_POT * nuclei_per_cm2            # 1/cm^2
    return L_eff_cm2 / 1e36                        # pb^-1


def protonBremSeries_v2(
        directory,
        pattern="Brem_120GeV_*.txt",
        det_angle=0.5/40.0,      # radians
        epsilon=1.0,             # σ ∝ ε^2
        lambda_idx=1,            # 0 → Λp=1.0 GeV, 1 → 1.5 GeV, 2 → 2.0 GeV
        N_POT=1.0e20,
        rho=7.87,
        L=16.8,
        A=55.845,
):
    """
    Read v2 single-track proton-brem files with columns:
      log10(theta), log10(p/GeV), sigma_L1[pb/bin], sigma_L1p5[pb/bin], sigma_L2[pb/bin]
    Returns:
      masses (sorted), tracks (same order), where 'tracks' already includes:
        - θ < det_angle acceptance on the *listed* particle
        - ε^2 scaling
        - ×2 (χ1/χ2 single-track symmetry)
        - multiplication by effective luminosity
    """
    files = sorted(glob.glob(os.path.join(directory, pattern)))
    if not files:
        print(f"⚠️ No v2 files found in {directory} matching {pattern}")
        return np.array([]), np.array([])

    L_eff_pb = _effective_lumi_pb(N_POT=N_POT, rho=rho, L=L, A=A)
    masses, tracks = [], []

    log10_theta_cut = np.log10(det_angle)
    for f in files:
        base = os.path.basename(f)
        # parse mass from tail, e.g. "..._0.0104.txt"
        try:
            mchi = float(base.split('_')[-1].replace('.txt', ''))
        except Exception:
            # ignore files that don't match naming
            continue

        arr = np.loadtxt(f, comments="#")
        if arr.ndim == 1:
            arr = arr[None, :]
        if arr.shape[1] < 5:
            print(f"⚠️ Unexpected column count in {f}; expected ≥5 (got {arr.shape[1]}). Skipping.")
            continue

        log10theta = arr[:, 0]
        sigmas = arr[:, 2:5]  # three Λp columns

        if lambda_idx not in (0, 1, 2):
            raise ValueError("lambda_idx must be 0 (1.0 GeV), 1 (1.5 GeV), or 2 (2.0 GeV).")

        # acceptance on the *listed* particle (single-track files already enforce partner out)
        mask = (log10theta < log10_theta_cut)

        # sum accepted σ for the chosen Λp (already per-bin, so simple sum)
        sigma_sum_pb = sigmas[mask, lambda_idx].sum()

        # single-track yields (events = tracks for singles),
        # then ×2 for χ1/χ2 symmetry; also σ × ε^2; and multiply by L_eff
        N_tracks = 2.0 * (epsilon**2) * L_eff_pb * sigma_sum_pb

        masses.append(mchi)
        tracks.append(N_tracks)

    if not masses:
        return np.array([]), np.array([])

    masses = np.asarray(masses)
    tracks = np.asarray(tracks)

    order = np.argsort(masses)
    return masses[order], tracks[order]


def _L_eff_pb(N_POT=1e20, rho=7.87, L=16.8, A=55.845):
    NA = 6.02214076e23
    return (N_POT * (rho * L / A) * NA) / 1e36  # pb^-1

def brem_singletrack_v2(
        directory,
        pattern="Brem_120GeV_*.txt",
        det_angle=0.5/40.0,   # radians
        epsilon=1.0,
        lambda_idx=1,         # 0→1.0 GeV, 1→1.5 GeV, 2→2.0 GeV
        N_POT=1e20, rho=7.87, L=16.8, A=55.845
):
    """
    Reads v2 SINGLE-TRACK files (second χ outside acceptance), with columns:
      log10(theta), log10(p/GeV), sigma_L1, sigma_L1p5, sigma_L2   [pb/bin]

    Returns: masses, N_tracks (already includes ε^2, L_eff, and ×2 χ1/χ2 symmetry)
    """
    import glob, os
    import numpy as np

    files = sorted(glob.glob(os.path.join(directory, pattern)))
    if not files:
        return np.array([]), np.array([])

    L_eff = _L_eff_pb(N_POT, rho, L, A)
    log10_theta_cut = np.log10(det_angle)

    masses, tracks = [], []
    for f in files:
        # parse mχ from filename suffix "..._<mass>.txt"
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

        mask = (log10theta < log10_theta_cut)
        sigma_sum_pb = sigmas[mask, lambda_idx].sum()

        # SINGLE-TRACK: events=tracks, then ×2 (χ1/χ2 symmetry)
        N_tracks = 2.0 * (epsilon**2) * L_eff * sigma_sum_pb

        masses.append(mchi)
        tracks.append(N_tracks)

    masses = np.asarray(masses); tracks = np.asarray(tracks)
    idx = np.argsort(masses)
    return masses[idx], tracks[idx]

def brem_doubletrack_v2(
        directory,
        pattern="Brem_120GeV_*.txt",
        det_angle=0.5/40.0,   # radians
        epsilon=1.0,
        lambda_idx=1,         # 0→1.0 GeV, 1→1.5 GeV, 2→2.0 GeV
        N_POT=1e20, rho=7.87, L=16.8, A=55.845
):
    """
    Reads v2 DOUBLE-TRACK files (both χ inside acceptance), with columns:
      log10(theta), log10(p/GeV), sigma_L1, sigma_L1p5, sigma_L2   [pb/bin]

    Returns: masses, N_tracks (already includes ε^2, L_eff, and ×2 tracks/event)
    """
    import glob, os
    import numpy as np

    files = sorted(glob.glob(os.path.join(directory, pattern)))
    if not files:
        return np.array([]), np.array([])

    L_eff = _L_eff_pb(N_POT, rho, L, A)
    log10_theta_cut = np.log10(det_angle)

    masses, tracks = [], []
    for f in files:
        # parse mχ from filename suffix "..._<mass>.txt"
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

        # DOUBLE-TRACK: two tracks per event (×2), no symmetry factor needed
        N_tracks = 2.0 * (epsilon**2) * L_eff * sigma_sum_pb

        masses.append(mchi)
        tracks.append(N_tracks)

    masses = np.asarray(masses); tracks = np.asarray(tracks)
    idx = np.argsort(masses)
    return masses[idx], tracks[idx]

def brem_single_hit_or_shared_v2(
        directory,
        pattern="Brem_120GeV_*.txt",
        det_angle=0.5/40.0,   # radians
        epsilon=1.0,
        lambda_idx=1,         # 0→1.0 GeV, 1→1.5 GeV, 2→2.0 GeV
        N_POT=1e20, rho=7.87, L=16.8, A=55.845
):
    """
    Reads v2 DOUBLE-TRACK files (both χ inside acceptance), with columns:
      log10(theta), log10(p/GeV), sigma_L1, sigma_L1p5, sigma_L2   [pb/bin]

    Returns: masses, N_tracks (already includes ε^2, L_eff, and ×2 tracks/event)
    """
    import glob, os
    import numpy as np

    files = sorted(glob.glob(os.path.join(directory, pattern)))
    if not files:
        return np.array([]), np.array([])

    L_eff = _L_eff_pb(N_POT, rho, L, A)
    log10_theta_cut = np.log10(det_angle)

    masses, tracks = [], []
    for f in files:
        # parse mχ from filename suffix "..._<mass>.txt"
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

        # DOUBLE-TRACK: two tracks per event (×2), no symmetry factor needed
        N_tracks = (epsilon**2) * L_eff * sigma_sum_pb

        masses.append(mchi)
        tracks.append(N_tracks)

    masses = np.asarray(masses); tracks = np.asarray(tracks)
    idx = np.argsort(masses)
    return masses[idx], tracks[idx]


def dyProductionTOT(fileNamecross, filenamedy, NPOT, totalCrossSection):

    # totalCrossSection = 300e-3

    # Import Data from efficiency files

    ageo_datacross = np.loadtxt(fileNamecross)
    ageo_datady = np.loadtxt(filenamedy)

    mass = ageo_datady[:, 0]
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

    # cross_dy may be 1- or 2-D; take first column if needed
    if ageo_datacross.ndim == 1:
        cross_dy = ageo_datacross
    else:
        cross_dy = ageo_datacross[:, 0]  # change index if your cross section is in another column

    # Convert pb → cm^2 (if that's your target unit)
    cross_dy = cross_dy * 1e-12

    # Ensure aligned lengths (use shortest)
    n = min(mass.shape[0], cross_dy.shape[0], ageo_dy.shape[0])
    mass    = mass[:n]
    ageo_dy = ageo_dy[:n]
    cross_dy = cross_dy[:n]

    # Compute production
    dy_prod = NPOT * cross_dy / totalCrossSection * ageo_dy

    # Zero out values below cutoff WITHOUT shifting indices
    dy_prod[mass < cutoff] = np.nan

    return dy_prod

if __name__ == '__main__':
    # read the original mass file
    #    mass = np.array([])
    #    with open('mass.txt', 'r') as f:
    #        data = f.readline()
    #        for line in data:
    # values = line.strip().split()
    #            print(line)
    #            np.append(mass, float(line) )
    mass_dark = np.loadtxt('/Users/leobailloeul/Documents/coding/decaysimulation/plotting/sensitivity-plot/mchi_values.txt')
    mass_ship = np.loadtxt('/Users/leobailloeul/Documents/coding/decaysimulation/plotting/sensitivity-plot/mship_values.txt')
    # mass_ship = np.loadtxt('/Users/leobailloeul/Documents/coding/decaysimulation/decay/sensitivity-plot/mass_ship.txt')
    # get decay production
    NPOT_dark = 1e20
    NPOT2_dark = 1e17
    NPOT_ship = 4*1e19
    #   pi_numi, eta_numi, jpsi_numi, upsilon_numi = mesonDecayProduction(numi)
    #   pi_dune, eta_dune, jpsi_dune, upsilon_dune = mesonDecayProduction(dune)
    base = '/Users/leobailloeul/Documents/coding/decaysimulation/plotting/sensitivity-plot/'
    # mesonDecay_dark = list(mesonDecayProduction(base + 'total_efficiency_outputdecayPion574mSpin.txt',
    #                                             base + 'total_efficiency_outputdecayEta574mSpin.txt',
    #                                             base + 'total_efficiency_output2bodydecay-rho574mSpin.txt',
    #                                             base + 'total_efficiency_output2bodydecay-phi574mSpin.txt',
    #                                             base + 'total_efficiency_output2bodydecay-omega574mSpin.txt',
    #                                             base + 'total_efficiency_output2bodydecay-jsi574mSpin.txt',
    #                                             NPOT_dark, mass_dark))
    # mesonDecay_dark = list(mesonDecayProduction(base + 'total_efficiency_outputdecayPion.txt',
    #                                             base + 'total_efficiency_outputdecayEta.txt',
    #                                             base + 'total_efficiency_output2bodydecay-rho.txt',
    #                                             base + 'total_efficiency_output2bodydecay-phi.txt',
    #                                             base + 'total_efficiency_output2bodydecay-omega.txt',
    #                                             base + 'total_efficiency_output2bodydecay-jsi.txt',
    #                                             NPOT_dark, mass_dark))
    # mesonDecay_dark = list(mesonDecayProduction(base + 'total_efficiency_outputdecayPion.txt',
    #                                             base + 'ageo_darkquest_0.5m.txt',
    #                                             base + 'total_efficiency_output2bodydecay-rho.txt',
    #                                             base + 'total_efficiency_output2bodydecay-phi.txt',
    #                                             base + 'total_efficiency_output2bodydecay-omega.txt',
    #                                             base + 'total_efficiency_output2bodydecay-jsi.txt',
    #                                             NPOT_dark, mass_dark))
    # mesonDecay_ship = list(mesonDecayProduction(base + 'total_efficiency_outputdecayPion.txt',
    #                                             base + 'ageo_darkquest_0.5m.txt',
    #                                             base + 'total_efficiency_output2bodydecay-rho.txt',
    #                                             base + 'total_efficiency_output2bodydecay-phi.txt',
    #                                             base + 'total_efficiency_output2bodydecay-omega.txt',
    #                                             base + 'total_efficiency_output2bodydecay-jsi.txt',
    #                                             NPOT2_dark, mass_dark))
    mesonDecay_ship = list(mesonDecayProduction(base + 'total_efficiency_outputdecayPionDARK.txt',
                                                base + 'total_efficiency_outputdecayEtaDARK.txt',
                                                base + 'total_efficiency_output2bodydecay-rhoDARK.txt',
                                                base + 'total_efficiency_output2bodydecay-phiDARK.txt',
                                                base + 'total_efficiency_output2bodydecay-omegaDARK.txt',
                                                base + 'total_efficiency_output2bodydecay-jsiDARK.txt',
                                                NPOT_dark, mass_ship))

    # mesonDecay_ship = list(mesonDecayProduction(base + 'total_efficiency_outputdecayPion1m.txt',
    #                                             base + 'total_efficiency_outputdecayEta1m.txt',
    #                                             base + 'total_efficiency_output2bodydecay-rho1m.txt',
    #                                             base + 'total_efficiency_output2bodydecay-phi1m.txt',
    #                                             base + 'total_efficiency_output2bodydecay-omega1m.txt',
    #                                             base + 'total_efficiency_output2bodydecay-jsi1m.txt',
    #                                             NPOT2_dark, mass_dark))


    # mesonDecay_ship = list(mesonDecayProduction(base + 'total_efficiency_outputdecayPion1mSpin.txt',
    #                                             base + 'total_efficiency_outputdecayEta1mSpin.txt',
    #                                             base + 'total_efficiency_output2bodydecay-rho1mSpin.txt',
    #                                             base + 'total_efficiency_output2bodydecay-phi1mSpin.txt',
    #                                             base + 'total_efficiency_output2bodydecay-omega1mSpin.txt',
    #                                             base + 'total_efficiency_output2bodydecay-jsi1mSpin.txt',
    #                                             NPOT_ship, mass_ship))

    # print(mesonDecay_numi[6])

    # get DY production
    # dy_dark = dyProduction(base+'ageo_dy.txt', base+'dy_cross.txt', NPOT_dark, 300e-3)
    # dy_numi = dyProduction(base+'efficiency_dy_dark.txt', base+'dy_cross_dark.txt', NPOT)
    dy_ship = dyProduction(base+'dy_cross_091525.txt', base+'ageo_dy_dark_0.5m.txt', NPOT_dark, 40e-3)
    dy_shiptot = dyProductionTOT(base+'dy_cross_091525.txt', base+'ageo_dy_dark_0.5m.txt', NPOT_dark, 40e-3)
    # dy_ship = dyProduction(base+'efficiency_dy_ship_1m.txt', base+'dy_cross_ship.txt', NPOT_ship, 260e-3)

    # add paddings
    for i in range(len(mesonDecay_ship)):
        difference = mass_ship.size - mesonDecay_ship[i].size
        mesonDecay_ship[i] = np.pad(mesonDecay_ship[i], (0, difference), 'constant',  constant_values=0)
    # for i in range(len(mesonDecay_dark)):
    #     difference = mass_dark.size - mesonDecay_dark[i].size
    #     mesonDecay_dark[i] = np.pad(mesonDecay_dark[i], (0, difference), 'constant',  constant_values=0)

    difference1 = mass_ship.size - dy_shiptot.size
    # difference2 = mass_dark.size - dy_dark.size
    #
    # dy_ship = np.pad(dy_ship, (0, difference1), 'constant', constant_values=0)
    dy_shiptot = np.pad(dy_shiptot, (0, difference1), 'constant', constant_values=0)
    # dy_dark = np.pad(dy_dark, (0, difference2), 'constant', constant_values=0)
    #
    total_ship = dy_shiptot.copy()

    for i in range(len(mesonDecay_ship)):
        total_ship += mesonDecay_ship[i]

    # --- v2 proton-brem: SINGLE vs DOUBLE and SUM ---
    brem_dir_single = "/Users/leobailloeul/Documents/coding/decaysimulation/single_track_updated"   # update to your single-track folder
    brem_dir_double = "/Users/leobailloeul/Documents/coding/decaysimulation/double_track_updated"   # update to your double-track folder
    brem_dir_shared = "/Users/leobailloeul/Documents/coding/decaysimulation/single_hit_or_shared_bar"
    det_angle = 0.5/40.0
    epsilon = 1.0
    lambda_idx = 1  # 0→1.0 GeV, 1→1.5 GeV, 2→2.0 GeV

    # read single-track (tracks already)
    m_single, y_single = brem_singletrack_v2(
        brem_dir_single, det_angle=det_angle, epsilon=epsilon,
        lambda_idx=lambda_idx, N_POT=NPOT_dark, rho=7.87, L=16.8, A=55.845
    )

    # read double-track (tracks already)
    m_double, y_double = brem_doubletrack_v2(
        brem_dir_double, det_angle=det_angle, epsilon=epsilon,
        lambda_idx=lambda_idx, N_POT=NPOT_dark, rho=7.87, L=16.8, A=55.845
    )

    m_shared, y_shared = brem_single_hit_or_shared_v2(
        brem_dir_double, det_angle=det_angle, epsilon=epsilon,
        lambda_idx=lambda_idx, N_POT=NPOT_dark, rho=7.87, L=16.8, A=55.845
    )

    # put them on a common mass grid for a clean Sum (choose union or a reference grid)
    m_grid = np.unique(np.concatenate([m_single, m_double]))
    y_single_g = np.interp(m_grid, m_single, y_single, left=0.0, right=0.0)
    y_double_g = np.interp(m_grid, m_double, y_double, left=0.0, right=0.0)
    y_shared_g = np.interp(m_grid, m_shared, y_shared, left=0.0, right=0.0)

    # SUM = tracks_single + tracks_double   (both already “tracks”)
    y_sum = y_single_g + y_double_g
    total_ship += np.interp(mass_ship, m_shared, y_shared, left=0.0, right=0.0) #edited

# Add v2 brem to total (fast path if grids identical; else safe fallback)
#     if (mass_ship.shape == m_brem_v2.shape) and np.allclose(mass_ship, m_brem_v2, rtol=1e-12, atol=0.0):
#         total_ship += y_brem_v2
#     else:
#         total_ship += np.interp(mass_ship, m_brem_v2, y_brem_v2, left=0.0, right=0.0)

    # (remove any other np.allclose(mass_ship, m_brem_v2) checks below)




# for i in range(len(mesonDecay_dark)):
    #     total_dark = mesonDecay_dark[i]

    # print(mass.size)
    # print(dy_numi.size)
    # print(mesonDecay_numi[1].size)


    # # draw plots
    # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))
    # plt.subplots_adjust(left=0.08, right=0.92, wspace=0.3)
    # x_ticks = [10**i for i in range(-2, 2)]
    # x_tick_labels = ['$10^{-2}$', '$10^{-1}$', '1', '$10^{1}$']
    # y_ticks = [10**i for i in range(0, 19)]
    # y_tick_labels = ['$1$'] + [(lambda x : f'$10^{{{x}}}$' if x % 2 == 0 else ' ')(i) for i in range(1, 19)]
    # textbox_props = dict(boxstyle='round', facecolor='white', edgecolor='none', alpha=0.7)
    #
    #
    # ax1.margins(x=0, y=0)
    # ax1.set_xscale('log')
    # ax1.set_yscale('log')
    # ax1.set_xlabel(r'$m_{\chi}$ [$\mathrm{GeV}/\mathrm{c}^2$]')
    # ax1.set_xticks(x_ticks)
    # ax1.set_xticklabels(x_tick_labels)
    # ax1.set_ylabel(r'$N_{\chi} / \epsilon^{2}$')
    # ax1.set_ylim(1, 1e18)
    # ax1.set_yticks(y_ticks)
    # ax1.set_yticklabels(y_tick_labels)
    # ax1.plot(mass_ship, mesonDecay_ship[0], label = r'$\pi^{0}\to\gamma\chi\overline{\chi}$')
    # ax1.plot(mass_ship, mesonDecay_ship[1], label = r'$\eta\to\gamma\chi\overline{\chi}$')
    # ax1.plot(mass_ship, mesonDecay_ship[2], label = r'$J/\psi\to\chi\overline{\chi}$')
    # ax1.plot(mass_ship, mesonDecay_ship[3], label = r'$\Upsilon\to\chi\overline{\chi}$')
    # ax1.plot(mass_ship, mesonDecay_ship[2], label = r'$\rho\to\chi\overline{\chi}$')
    # ax1.plot(mass_ship, mesonDecay_ship[3], label = r'$\Phi\to\chi\overline{\chi}$')
    # ax1.plot(mass_ship, mesonDecay_ship[4], label = r'$\Omega\to\chi\overline{\chi}$')
    #
    #
    # # ax1.plot(mass_ship, dy_ship,            label = r'$q\overline{q}\to\gamma^{*}\to\chi\overline{\chi}$')
    # # ax1.plot(mass_ship, total_ship, color='black', linestyle='solid', linewidth=2, label='Total')
    # ax1.legend(loc = 'lower left', ncol = 2)
    # # ax1.text(0.55, 0.93, '$1$ Year at SHIP ($4 x 10^{19}$ POT) \n $100$ m, $1$ m radius cylindrical detector', transform=ax1.transAxes, fontsize=10, verticalalignment='center', multialignment = 'left',  bbox=textbox_props)
    # ax1.text(0.55, 0.93, '$1$ Year at DarkQuest ($10^{18}$ POT) \n $40$ m, $0.5$ m radius cylindrical detector', transform=ax1.transAxes, fontsize=10, verticalalignment='center', multialignment = 'left', bbox=textbox_props)
    #
    #
    # ax2.margins(x=0, y=0)
    # ax2.set_xscale('log')
    # ax2.set_yscale('log')
    # ax2.set_xlabel(r'$m_{\chi}$ [$\mathrm{GeV}/\mathrm{c}^2$]')
    # ax2.set_xticks(x_ticks)
    # ax2.set_xticklabels(x_tick_labels)
    # ax2.set_ylabel(r'$N_{\chi} / \epsilon^{2}$')
    # ax2.set_xlim(1e-2, 1e-1)
    # ax2.set_ylim(1, 1e16)
    # ax2.set_yticks(y_ticks)
    # ax2.set_yticklabels(y_tick_labels)
    # ax2.plot(mass_dark, mesonDecay_dark[0], label = r'$\pi^{0}\to\gamma\chi\overline{\chi}$')
    # # ax2.plot(mass_dark, mesonDecay_dark[1], label = r'$\eta\to\gamma\chi\overline{\chi}$')
    # # ax2.plot(mass_dark, mesonDecay_dark[2], label = r'$J/\psi\to\chi\overline{\chi}$')
    # # ax2.plot(mass_dark, mesonDecay_dark[3], label = r'$\Upsilon\to\chi\overline{\chi}$')
    #
    # # ax2.plot(mass_dark, mesonDecay_dark[2], label = r'$J/\psi\to\chi\overline{\chi}$')
    # # ax2.plot(mass_dark, mesonDecay_dark[3], label = r'$\Upsilon\to\chi\overline{\chi}$')
    # # ax2.plot(mass_dark, mesonDecay_dark[4], label = r'$\Omega\to\chi\overline{\chi}$')
    #
    # # ax2.plot(mass_dark, dy_dark,            label = r'$q\overline{q}\to\gamma^{*}\to\chi\overline{\chi}$')
    # # ax2.plot(mass_dark, total_dark, color='black', linestyle='solid', linewidth=2, label='Total')
    # ax2.legend(loc = 'lower left', ncol = 2)
    # ax2.text(0.55, 0.93, '$1$ Year at ArgoNeut ($10^{20}$ POT) \n $975$ m, $0.23$ m radius cylindrical detector', transform=ax2.transAxes, fontsize=10, verticalalignment='center', multialignment = 'left', bbox=textbox_props)
    #
    #
    #
    # # plt.show()
    # # save plot as pgf format
    # fig.savefig('MCP-production.png')
    # draw plot
# --- ONE-PANEL, WIDE PLOT LIKE YOUR EXAMPLE ---
fig, ax = plt.subplots(figsize=(8.09, 5)) # wide like the screenshot

# axes + scales
ax.set_xscale('log')
ax.set_yscale('log')
ax.margins(x=0, y=0)
ax.set_xlabel(r'$m_{\chi}\, [\mathrm{GeV}/c^{2}]$', fontsize=13)
ax.set_ylabel(r'$N_{\chi}/\epsilon^{2}$', fontsize=13)

# nice ticks (powers of 10 with sparse labels)
x_ticks = [10**i for i in range(-2, 2)]
x_tick_labels = [r'$10^{-2}$', r'$10^{-1}$', r'$1$', r'$10^{1}$']
ax.set_xticks(x_ticks)
ax.set_xticklabels(x_tick_labels)

y_ticks = [10**i for i in range(0, 20)]
y_tick_labels = [r'$1$'] + [ (rf'$10^{{{i}}}$' if i % 2 == 0 else ' ') for i in range(1, 20) ]
ax.set_yticks(y_ticks)
ax.set_yticklabels(y_tick_labels)
ax.set_ylim(1e6, 1e20)



# thin grid (optional, helps readability)
# ax.grid(True, which='both', linewidth=0.5, alpha=0.25)

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

print(dy_ship)
print(dy_shiptot)

lambda_labels = {0: "Λp=1.0 GeV", 1: "Λp=1.5 GeV", 2: "Λp=2.0 GeV"}

# ax.plot(m_single, y_single, lw=3, label=rf'$p$ brem (single, {lambda_labels[lambda_idx]})')
# ax.plot(m_double, y_double, lw=3, ls='--', label=rf'$p$ brem (double, {lambda_labels[lambda_idx]})')
ax.plot(m_shared, y_shared, lw=3, ls='--', label=rf'$p$ brem')
# ax.plot(m_shared, y_shared, lw=3, ls='--', label=rf'$p$ brem (shared, {lambda_labels[lambda_idx]})')
# ax.plot(m_grid,  y_sum,    lw=3, color='k', label='p brem (sum)')

# paths/inputs
# brem_dir_st = "/Users/leobailloeul/Documents/coding/decaysimulation/single_track"   # your v2 single-track folder
# det_angle = 0.5/40.0
# epsilon = 1.0
#
# lambda_labels = {0: "Λp=1.0 GeV", 1: "Λp=1.5 GeV", 2: "Λp=2.0 GeV"}
#
# # central (Λp=1.5 GeV) as solid
# m_st_1p5, y_st_1p5 = brem_singletrack_v2(brem_dir_st, det_angle=det_angle,
#                                          epsilon=epsilon, lambda_idx=1,
#                                          N_POT=NPOT_dark)
# ax.plot(m_st_1p5, y_st_1p5, lw=3,
#         label=rf'$p$ bremsstrahlung (single, {lambda_labels[1]})')
#
# # Λp=1.0 and 2.0 as dashed
# m_st_1p0, y_st_1p0 = brem_singletrack_v2(brem_dir_st, det_angle=det_angle,
#                                          epsilon=epsilon, lambda_idx=0,
#                                          N_POT=NPOT_dark)
# ax.plot(m_st_1p0, y_st_1p0, lw=3, ls='--',
#         label=rf'$p$ bremsstrahlung (single, {lambda_labels[0]})')
#
# m_st_2p0, y_st_2p0 = brem_singletrack_v2(brem_dir_st, det_angle=det_angle,
#                                          epsilon=epsilon, lambda_idx=2,
#                                          N_POT=NPOT_dark)
# ax.plot(m_st_2p0, y_st_2p0, lw=3, ls='--',
#         label=rf'$p$ bremsstrahlung (single, {lambda_labels[2]})')
#
# print(y_st_1p0.shape, y_st_2p0.shape)


ax.plot(mass_ship, total_ship, color='black', lw=3, label='Total')

# legend styled like your example
# leg = ax.legend(loc='lower left', ncol=2, frameon=True, fancybox=True, framealpha=0.85)
# for line in leg.get_lines():
#     line.set_linewidth(2.0)

leg = ax.legend(
    loc='lower left',
    ncol=2,
    frameon=True,
    fancybox=True,
    framealpha=0.85,

    fontsize=12.5,       # was default ~12 → ~10% bump
    handlelength=2.2,  # slightly longer lines
    markerscale=1.1,   # 10% thicker line appearance
)

for line in leg.get_lines():
    line.set_linewidth(2.2)   # tiny increase from 2.0 → 2.2


# top-right annotation (two lines)
textbox_props = dict(boxstyle='round', facecolor='white', edgecolor='none', alpha=0.8)
ax.text(0.48, 0.93,
        r'$1$ Year at DarkQuest ($10^{20}$ POT)' '\n'
        r'$40$ m, $0.5$ m radius cylindrical detector',
        transform=ax.transAxes, ha='left', va='center', bbox=textbox_props, fontsize=12.5)

# save full-page-ish
fig.tight_layout()
fig.savefig('MCP-production.pdf', dpi=300, bbox_inches='tight')
# fig.savefig('MCP-production.png', dpi=300, bbox_inches='tight')