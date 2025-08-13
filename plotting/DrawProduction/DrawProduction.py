import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import ScalarFormatter
from scipy.integrate import quad
# define matplotlib pgf output settings
# matplotlib.use("pgf")

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
    # ageo_datarho = np.loadtxt(fileNameRho)
    # ageo_dataphi = np.loadtxt(fileNamePhi)
    # ageo_dataomega = np.loadtxt(fileNameOmega)
    ageo_datajpsi = np.loadtxt(fileNameJpsi)

    # Separate efficiency and mass
    # mass_rho = ageo_datarho[:,0]
    # mass_omega = ageo_dataomega[:,0]
    # mass_phi = ageo_dataphi[:,0]
    mass_pi = ageo_datapi[:,0]
    mass_eta = ageo_dataeta[:,0]
    mass_jpsi = ageo_datajpsi[:,0]
    mass_upsilon = fileMass[:]

    # ageo_rho = ageo_datarho[:,1]
    # ageo_omega = ageo_dataomega[:,1]
    # ageo_phi = ageo_dataphi[:,1]
    ageo_pi = ageo_datapi[:,1]
    ageo_jpsi = ageo_datajpsi[:,1]
    # special case for abnormal file structures
    if fileNameEta == '/Users/leobailloeul/Documents/coding/decaysimulation/plotting/sensitivity-plot/ageo_darkquest_0.5m.txt':
        ageo_eta = ageo_dataeta[:,2]
    else:
        ageo_eta = ageo_dataeta[:,1]
    if fileNameJpsi == '/Users/leobailloeul/Documents/coding/decaysimulation/plotting/sensitivity-plot/total_efficiency_output2bodydecay-jsi.txt':
        ageo_upsilon = np.full(mass_upsilon.size, 0.011462)
    else:
        ageo_upsilon = np.full(mass_upsilon.size, 0.000474)

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
    # phi = np.zeros(np.shape(mass_phi))
    # omega = np.zeros(np.shape(mass_omega))
    jpsi = np.zeros(np.shape(mass_jpsi))
    # rho = np.zeros(np.shape(mass_rho))
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
    # rho[mass_rho < m_rho/2] = NPOT * ageo_rho[mass_rho < m_rho/2] * 2  * c_rho * branch_rho * I2( mass_rho[mass_rho < m_rho/2] ** 2 / m_rho ** 2, m_e ** 2 / m_rho ** 2  )
    # omega[mass_omega < m_omega/2] = NPOT * ageo_omega[mass_omega < m_omega/2] * 2  * c_omega * branch_omega * I2( mass_omega[mass_omega < m_omega/2] ** 2 / m_omega ** 2, m_e ** 2 / m_omega ** 2  )
    # phi[mass_phi < m_phi/2] = NPOT * ageo_phi[mass_phi < m_phi/2] * 2  * c_phi * branch_phi * I2( mass_phi[mass_phi < m_phi/2] ** 2 / m_phi ** 2, m_e ** 2 / m_phi ** 2  )
    upsilon[mass_upsilon < m_upsilon/2] = NPOT * ageo_upsilon[mass_upsilon < m_upsilon/2] * 2  * c_upsilon * branch_upsilon * I2( mass_upsilon[mass_upsilon < m_upsilon/2] ** 2 / m_upsilon ** 2, m_e ** 2 / m_upsilon ** 2  )

    return pi, eta, jpsi, upsilon #, rho, omega, phi

def dyProduction(fileNamecross, filenamedy, NPOT, totalCrossSection):

    # totalCrossSection = 300e-3

    # Import Data from efficiency files
    ageo_datacross = np.loadtxt(fileNamecross)
    ageo_datady = np.loadtxt(filenamedy)

    cross_dy = ageo_datacross[:,1] * 1e-12 # pico barn
    ageo_dy = ageo_datady[:]

    return NPOT * cross_dy / totalCrossSection * ageo_dy

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
    NPOT_dark = 1e21
    NPOT2_dark = 1e17
    NPOT_ship = 4*1e19
    #   pi_numi, eta_numi, jpsi_numi, upsilon_numi = mesonDecayProduction(numi)
    #   pi_dune, eta_dune, jpsi_dune, upsilon_dune = mesonDecayProduction(dune)
    base = '/Users/leobailloeul/Documents/coding/decaysimulation/plotting/sensitivity-plot/'
    mesonDecay_dark = list(mesonDecayProduction(base + 'total_efficiency_outputdecayPion574mSpin.txt',
                                                base + 'total_efficiency_outputdecayEta574mSpin.txt',
                                                base + 'total_efficiency_output2bodydecay-rho574mSpin.txt',
                                                base + 'total_efficiency_output2bodydecay-phi574mSpin.txt',
                                                base + 'total_efficiency_output2bodydecay-omega574mSpin.txt',
                                                base + 'total_efficiency_output2bodydecay-jsi574mSpin.txt',
                                                NPOT_dark, mass_dark))
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
    mesonDecay_ship = list(mesonDecayProduction(base + 'DQ-total_efficiency_outputdecayPion1m.txt',
                                                base + 'DQ-total_efficiency_outputdecayEta1m.txt',
                                                base + 'total_efficiency_output2bodydecay-rho1m.txt',
                                                base + 'total_efficiency_output2bodydecay-phi1m.txt',
                                                base + 'total_efficiency_output2bodydecay-omega1m.txt',
                                                base + 'total_efficiency_output2bodydecay-jsi1m.txt',
                                                NPOT2_dark, mass_ship))

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
    dy_dark = dyProduction(base+'ageo_dy.txt', base+'dy_cross.txt', NPOT_dark, 300e-3)
    # dy_numi = dyProduction(base+'efficiency_dy_dark.txt', base+'dy_cross_dark.txt', NPOT)
    dy_ship = dyProduction(base+'ageo_dy.txt', base+'dy_cross.txt', NPOT2_dark, 300e-3)

    # dy_ship = dyProduction(base+'efficiency_dy_ship_1m.txt', base+'dy_cross_ship.txt', NPOT_ship, 260e-3)

    # add paddings
    for i in range(len(mesonDecay_ship)):
        difference = mass_ship.size - mesonDecay_ship[i].size
        mesonDecay_ship[i] = np.pad(mesonDecay_ship[i], (0, difference), 'constant',  constant_values=0)
    for i in range(len(mesonDecay_dark)):
        difference = mass_dark.size - mesonDecay_dark[i].size
        mesonDecay_dark[i] = np.pad(mesonDecay_dark[i], (0, difference), 'constant',  constant_values=0)

    # difference1 = mass_ship.size - dy_ship.size
    # difference2 = mass_dark.size - dy_dark.size
    #
    # dy_ship = np.pad(dy_ship, (0, difference1), 'constant', constant_values=0)
    # dy_dark = np.pad(dy_dark, (0, difference2), 'constant', constant_values=0)
    #
    # total_ship = dy_ship
    # total_dark = dy_dark

    for i in range(len(mesonDecay_ship)):
        total_ship = mesonDecay_ship[i]
    for i in range(len(mesonDecay_dark)):
        total_dark = mesonDecay_dark[i]

    # print(mass.size)
    # print(dy_numi.size)
    # print(mesonDecay_numi[1].size)


    # draw plots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))
    plt.subplots_adjust(left=0.08, right=0.92, wspace=0.3)
    x_ticks = [10**i for i in range(-2, 2)]
    x_tick_labels = ['$10^{-2}$', '$10^{-1}$', '1', '$10^{1}$']
    y_ticks = [10**i for i in range(0, 19)]
    y_tick_labels = ['$1$'] + [(lambda x : f'$10^{{{x}}}$' if x % 2 == 0 else ' ')(i) for i in range(1, 19)]
    textbox_props = dict(boxstyle='round', facecolor='white', edgecolor='none', alpha=0.7)


    ax1.margins(x=0, y=0)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel(r'$m_{\chi}$ [$\mathrm{GeV}/\mathrm{c}^2$]')
    ax1.set_xticks(x_ticks)
    ax1.set_xticklabels(x_tick_labels)
    ax1.set_ylabel(r'$N_{\chi} / \epsilon^{2}$')
    ax1.set_ylim(1, 1e18)
    ax1.set_yticks(y_ticks)
    ax1.set_yticklabels(y_tick_labels)
    ax1.plot(mass_ship, mesonDecay_ship[0], label = r'$\pi^{0}\to\gamma\chi\overline{\chi}$')
    ax1.plot(mass_ship, mesonDecay_ship[1], label = r'$\eta\to\gamma\chi\overline{\chi}$')
    ax1.plot(mass_ship, mesonDecay_ship[2], label = r'$J/\psi\to\chi\overline{\chi}$')
    ax1.plot(mass_ship, mesonDecay_ship[3], label = r'$\Upsilon\to\chi\overline{\chi}$')
    # ax1.plot(mass_ship, mesonDecay_ship[2], label = r'$\rho\to\chi\overline{\chi}$')
    # ax1.plot(mass_ship, mesonDecay_ship[3], label = r'$\Phi\to\chi\overline{\chi}$')
    # ax1.plot(mass_ship, mesonDecay_ship[4], label = r'$\Omega\to\chi\overline{\chi}$')


    # ax1.plot(mass_ship, dy_ship,            label = r'$q\overline{q}\to\gamma^{*}\to\chi\overline{\chi}$')
    # ax1.plot(mass_ship, total_ship, color='black', linestyle='solid', linewidth=2, label='Total')
    ax1.legend(loc = 'lower left', ncol = 2)
    # ax1.text(0.55, 0.93, '$1$ Year at SHIP ($4 x 10^{19}$ POT) \n $100$ m, $1$ m radius cylindrical detector', transform=ax1.transAxes, fontsize=10, verticalalignment='center', multialignment = 'left',  bbox=textbox_props)
    ax1.text(0.55, 0.93, '$1$ Year at LongQuest ($10^{17}$ POT) \n $40$ m, $0.11$ m radius cylindrical detector', transform=ax1.transAxes, fontsize=10, verticalalignment='center', multialignment = 'left', bbox=textbox_props)


    ax2.margins(x=0, y=0)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlabel(r'$m_{\chi}$ [$\mathrm{GeV}/\mathrm{c}^2$]')
    ax2.set_xticks(x_ticks)
    ax2.set_xticklabels(x_tick_labels)
    ax2.set_ylabel(r'$N_{\chi} / \epsilon^{2}$')
    ax2.set_ylim(1, 1e18)
    ax2.set_yticks(y_ticks)
    ax2.set_yticklabels(y_tick_labels)
    ax2.plot(mass_dark, mesonDecay_dark[0], label = r'$\pi^{0}\to\gamma\chi\overline{\chi}$')
    ax2.plot(mass_dark, mesonDecay_dark[1], label = r'$\eta\to\gamma\chi\overline{\chi}$')
    ax2.plot(mass_dark, mesonDecay_dark[2], label = r'$J/\psi\to\chi\overline{\chi}$')
    ax2.plot(mass_dark, mesonDecay_dark[3], label = r'$\Upsilon\to\chi\overline{\chi}$')

    # ax2.plot(mass_dark, mesonDecay_dark[2], label = r'$J/\psi\to\chi\overline{\chi}$')
    # ax2.plot(mass_dark, mesonDecay_dark[3], label = r'$\Upsilon\to\chi\overline{\chi}$')
    # ax2.plot(mass_dark, mesonDecay_dark[4], label = r'$\Omega\to\chi\overline{\chi}$')

    # ax2.plot(mass_dark, dy_dark,            label = r'$q\overline{q}\to\gamma^{*}\to\chi\overline{\chi}$')
    # ax2.plot(mass_dark, total_dark, color='black', linestyle='solid', linewidth=2, label='Total')
    ax2.legend(loc = 'lower left', ncol = 2)
    ax2.text(0.55, 0.93, '$1$ Year at DUNE ($10^{21}$ POT) \n $574$ m, $1$ m radius cylindrical detector', transform=ax2.transAxes, fontsize=10, verticalalignment='center', multialignment = 'left', bbox=textbox_props)



    # plt.show()
    # save plot as pgf format
    fig.savefig('MCP-production.png')