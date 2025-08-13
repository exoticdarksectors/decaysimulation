import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.interpolate import interp1d
import os

# Use Agg backend if running headless
matplotlib.use('Agg')
matplotlib.rcParams.update({
    'font.family': 'serif',
    'text.usetex': False,
})

def dalitz_yields(eff_file, mass_cut, c_const, branch_frac, alpha, NPOT, N_in):
    """
    Compute Dalitz yields and errors for a given parent spectrum.
    """
    data = np.loadtxt(eff_file, skiprows=0)
    mass   = data[:,0]
    eff    = data[:,1]
    counts = data[:,2]

    def I3_integrand(z, x):
        return (2/(3*np.pi)
                * np.sqrt(1 - 4*x/z)
                * (1 - z)**3
                * (2*x + z)/z**2)
    def I3(x):
        return quad(I3_integrand, 4*x, 1, args=(x,))[0]
    vI3 = np.vectorize(I3)

    mask = mass < mass_cut
    x = mass[mask]**2 / (2*mass_cut)**2
    pref = 2 * c_const * branch_frac * alpha

    y    = np.zeros_like(mass)
    yerr = np.zeros_like(mass)
    y[mask]    = NPOT * eff[mask] * pref * vI3(x)
    yerr[mask] = NPOT * (np.sqrt(counts[mask]) / N_in) * pref * vI3(x)

    return mass, y, yerr


def read_spectrum(filename):
    data = np.genfromtxt(
        filename, delimiter=',', comments='#', invalid_raise=False
    )
    mass   = data[:,0]
    counts = data[:,1]
    errs   = np.sqrt(counts)
    return mass, counts, errs


def interpolate_spectrum(x_src, y_src, yerr_src, x_target):
    f_y    = interp1d(x_src,    y_src,    kind='linear', bounds_error=False, fill_value=np.nan)
    f_yerr = interp1d(x_src, yerr_src, kind='linear', bounds_error=False, fill_value=np.nan)
    return f_y(x_target), f_yerr(x_target)


def compute_ratio(y_num, yerr_num, y_den, yerr_den):
    ratio    = y_num / y_den
    frac_err = np.sqrt((yerr_num/y_num)**2 + (yerr_den/y_den)**2)
    return ratio, ratio * frac_err


def main():
    # Constants
    alpha    = 1/137
    NPOT     = 1e20
    N_in_pi  = 4.6773361e7
    N_in_eta = 5.2799435e7

    # File paths (customize as needed)
    pi_eff_file  = '/Users/leobailloeul/Documents/coding/decaysimulation/plotting/sensitivity-plot/ArgoNeutTwobyTwototal_efficiency_outputdecayPion.txt'
    eta_eff_file = '/Users/leobailloeul/Documents/coding/decaysimulation/plotting/sensitivity-plot/TwoByTwototal_efficiency_outputdecayEtaDBG.txt'
    pi_data_file = '/Users/leobailloeul/Desktop/ArgoNeutPi0Data.txt'
    eta_data_file= '/Users/leobailloeul/Desktop/ArgoNeutEtaDataSet.txt'
    out_png      = 'combined_comparison.png'

    # Dalitz yields
    mass_pi,  y_pi,  yerr_pi  = dalitz_yields(pi_eff_file, 0.135/2, 4.6, 0.98, alpha, NPOT, N_in_pi)
    mass_eta, y_eta, yerr_eta = dalitz_yields(eta_eff_file,0.548/2, 0.51,0.39, alpha, NPOT, N_in_eta)

    # Experimental spectra
    m_pi_data, pi_counts, pi_errs = read_spectrum(pi_data_file)
    m_eta_data, eta_counts, eta_errs = read_spectrum(eta_data_file)

    pi_counts = 1e4 * pi_counts
    eta_counts = 1e4 * eta_counts

    # Interpolations for ratios
    # Pythia π⁰ onto homemade π grid
    pythia_pi_on_pi, pythia_pi_err_on_pi = interpolate_spectrum(
        m_pi_data, pi_counts, pi_errs, mass_pi)
    # Pythia η onto homemade η grid
    pythia_eta_on_eta, pythia_eta_err_on_eta = interpolate_spectrum(
        m_eta_data, eta_counts, eta_errs, mass_eta)

    # Compute ratios
    # π⁰: homemade / Pythia
    pi_ratio, pi_ratio_err = compute_ratio(
        y_pi, yerr_pi, pythia_pi_on_pi, pythia_pi_err_on_pi)
    # η: homemade / Pythia
    eta_ratio, eta_ratio_err = compute_ratio(
        y_eta, yerr_eta, pythia_eta_on_eta, pythia_eta_err_on_eta)

    # Plotting
    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, sharex=True,
        gridspec_kw={'height_ratios':[3,1]},
        figsize=(8,6)
    )

    # Top: yields and spectra
    ax_top.set_yscale('log')
    ax_top.errorbar(mass_pi, y_pi,  yerr=yerr_pi,  fmt='o-', label='Homemade π⁰', capsize=3)
    ax_top.errorbar(mass_eta,y_eta, yerr=yerr_eta, fmt='s--', label='Homemade η', capsize=3)
    ax_top.errorbar(m_pi_data, pi_counts, yerr=pi_errs, fmt='^', label='Pythia π⁰', capsize=3)
    ax_top.errorbar(m_eta_data, eta_counts, yerr=eta_errs, fmt='v', label='Pythia η', capsize=3)
    ax_top.set_ylabel(r'$N_{\chi} (\epsilon = 1)$')
    ax_top.legend(loc='upper right')
    ax_top.grid(which='both', ls='--', lw=0.5)

    # Bottom: homemade / Pythia ratios
    ax_bot.set_xscale('log')
    ax_bot.set_yscale('log')
    ax_bot.errorbar(mass_pi, pi_ratio**(-1), yerr=pi_ratio_err, fmt='o-', label='π⁰ Pythia/Homemade', capsize=3)
    ax_bot.errorbar(mass_eta, eta_ratio**(-1), yerr=eta_ratio_err, fmt='s--', label='η Pythia/Homemade', capsize=3)
    ax_bot.axhline(1.0, ls='--', color='gray')
    ax_bot.set_ylabel('Pythia/Homemade')
    ax_bot.set_xlabel(r'$m_{\chi}$ [GeV/$c^2$]')
    ax_bot.set_ylim(100, 100000)
    ax_bot.legend(loc='upper right')
    ax_bot.grid(which='both', ls='--', lw=0.5)

    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    print(f"Saved combined plot to: {out_png}")

if __name__ == '__main__':
    main()
