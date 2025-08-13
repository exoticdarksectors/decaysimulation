import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from matplotlib import rc
from matplotlib.ticker import MaxNLocator
from matplotlib.pyplot import figure
from matplotlib.colors import LinearSegmentedColormap, LogNorm
import pandas as pd

#
def mesonDecayProduction(fileNamePi, fileNameEta, fileNameRho, fileNamePhi, fileNameOmega, fileNameJpsi, NPOT, mass, charge, a, fileMass, N_gamma=2.5e6): # num of layers
    # EM constant
    alpha = 1.0 / 137.0

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
    c_pi = 3.98
    c_eta = 0.51
    c_rho = 0.24
    c_omega = 0.24
    c_phi = 4.9e-03
    c_jpsi = 3.81e-5
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
    mass_upsilon = fileMass[:]

    # Separate efficiency and mass
    ageo_rho = ageo_datarho[:,1]
    ageo_omega = ageo_dataomega[:,1]
    ageo_phi = ageo_dataphi[:,1]
    ageo_pi = ageo_datapi[:,1]
    if fileNameEta == '/Users/leobailloeul/Documents/coding/decaysimulation/plotting/sensitivity-plot/ageo_darkquest_0.5m.txt':
        ageo_eta = ageo_dataeta[:,2]
    else:
        ageo_eta = ageo_dataeta[:,1]
    ageo_jpsi = ageo_datajpsi[:,1]
    ageo_upsilon = np.full(mass_upsilon.size, 0.010127)

    list_ageo = [ageo_rho, ageo_omega, ageo_phi, ageo_pi, ageo_eta, ageo_jpsi, ageo_upsilon]
    print(fileMass.size)
    print(charge.size)
    #padding
    for i in range(len(list_ageo)):
        difference = fileMass.size - np.size(list_ageo[i])
        list_ageo[i] = np.pad(list_ageo[i], (0, difference), 'constant',  constant_values=0)
        list_ageo[i] = np.tile(list_ageo[i], (int(charge.size / fileMass.size), 1))

    # Re-assign back to the original variables
    ageo_rho, ageo_omega, ageo_phi, ageo_pi, ageo_eta, ageo_jpsi, ageo_upsilon = list_ageo

    # define phase space integrals
    def I2(x, y):
        return ((1 + 2 * x) * np.sqrt(1 - 4 * x)) / ((1 + 2 * y) * np.sqrt(1 - 4 * y))

    def I3_integrand(z, x):
        return 2 / (3 * 3.14) * np.sqrt(1 - 4 * x / z)  * ((1 - z) ** 3) * (2 * x + z) / (z ** 2)

    def I3(x):
        return quad(I3_integrand, 4 * x, 1, args=(x))[0]  # integrate
    v_I3 = np.vectorize(I3) # vectorize

    pi = np.zeros(np.shape(masses))
    eta = np.zeros(np.shape(masses))
    phi = np.zeros(np.shape(masses))
    omega = np.zeros(np.shape(masses))
    jpsi = np.zeros(np.shape(masses))
    rho = np.zeros(np.shape(masses))
    upsilon = np.zeros(np.shape(masses))

    # Dalitz decay
    pi[mass < m_pi/2] = ageo_pi[mass < m_pi/2] * 2 * c_pi * branch_pi * alpha * v_I3(mass[mass < m_pi/2] ** 2 / m_pi ** 2) * (1 - np.exp(-N_gamma * charge[mass < m_pi / 2] ** 2)) ** a * charge[mass < m_pi / 2] ** 2

    eta[mass < m_eta/2] = ageo_eta[mass < m_eta/2] * 2 * c_eta * branch_eta * alpha * v_I3(mass[mass < m_eta/2] ** 2 / m_eta ** 2) * (1 - np.exp(-N_gamma * charge[mass < m_eta / 2] ** 2)) ** a * charge[mass < m_eta / 2] ** 2

    # Direct decay
    jpsi[mass < m_jpsi/2] = ageo_jpsi[mass < m_jpsi/2] * 2  * c_jpsi * branch_jpsi * I2( mass[mass < m_jpsi/2] ** 2 / m_jpsi ** 2, m_e ** 2 / m_jpsi ** 2  ) * (1 - np.exp(-N_gamma * charge[mass < m_jpsi / 2] ** 2)) ** a * charge[mass < m_jpsi / 2] ** 2

    rho[mass < m_rho/2] = ageo_rho[mass < m_rho/2] * 2  * c_rho * branch_rho * I2( mass[mass < m_rho/2] ** 2 / m_rho ** 2, m_e ** 2 / m_rho ** 2  ) * (1 - np.exp(-N_gamma * charge[mass < m_rho / 2] ** 2)) ** a * charge[mass < m_rho / 2] ** 2

    omega[mass < m_omega/2] = ageo_omega[mass < m_omega/2] * 2  * c_omega * branch_omega * I2( mass[mass < m_omega/2] ** 2 / m_omega ** 2, m_e ** 2 / m_omega ** 2  ) * (1 - np.exp(-N_gamma * charge[mass < m_omega / 2] ** 2)) ** a * charge[mass < m_omega / 2] ** 2

    phi[mass < m_phi/2] = ageo_phi[mass < m_phi/2] * 2  * c_phi * branch_phi * I2( mass[mass < m_phi/2] ** 2 / m_phi ** 2, m_e ** 2 / m_phi ** 2  ) * (1 - np.exp(-N_gamma * charge[mass < m_phi / 2] ** 2)) ** a * charge[mass < m_phi / 2] ** 2

    upsilon[mass < m_upsilon/2] = ageo_upsilon[mass < m_upsilon/2] * 2  * c_upsilon * branch_upsilon * I2( mass[mass < m_upsilon/2] ** 2 / m_upsilon ** 2, m_e ** 2 / m_upsilon ** 2  )* (1 - np.exp(-N_gamma * charge[mass < m_upsilon / 2] ** 2)) ** a * charge[mass < m_upsilon / 2] ** 2

    sensitivity = pi + eta + jpsi + rho + omega + phi + upsilon

    return sensitivity * NPOT

def dyProduction(fileNamecross, filenamedy, NPOT, charge, a, N_gamma=2.5e6):

    totalCrossSection = 300e-3

    # Import Data from efficiency files
    ageo_datacross = np.loadtxt(fileNamecross)
    ageo_datady = np.loadtxt(filenamedy)

    cross_dy = ageo_datacross[:,1] * 1e-12 # pico barn
    ageo_dy = ageo_datady[:]

    return NPOT * cross_dy / totalCrossSection * ageo_dy * (1 - np.exp(-N_gamma * charge ** 2)) ** a * charge ** 2

base = '/Users/leobailloeul/Documents/coding/decaysimulation/plotting/sensitivity-plot/'
NPOT = 4*1e19
mass = np.loadtxt(base + 'mass_ship.txt')
charge = np.logspace(-5, 1, 500)
masses, charges = np.meshgrid(mass, charge)
print(charges.size)

mesonDecay_dune = list(mesonDecayProduction(base + 'total_efficiency_outputdecayPion1mSpin.txt',
                              base + 'total_efficiency_outputdecayEta1mSpin.txt',
                              base + 'total_efficiency_output2bodydecay-rho1mSpin.txt',
                              base + 'total_efficiency_output2bodydecay-phi1mSpin.txt',
                              base + 'total_efficiency_output2bodydecay-omega1mSpin.txt',
                              base + 'total_efficiency_output2bodydecay-jsi1mSpin.txt',
                              NPOT, masses, charges, 2, mass))
# mesonDecay_dune = list(mesonDecayProduction(base + 'total_efficiency_outputdecayPion.txt',
#                                                 base + 'ageo_darkquest_0.5m.txt',
#                                                 base + 'total_efficiency_output2bodydecay-rho.txt',
#                                                 base + 'total_efficiency_output2bodydecay-phi.txt',
#                                                 base + 'total_efficiency_output2bodydecay-omega.txt',
#                                                 base + 'total_efficiency_output2bodydecay-jsi.txt',
#                                                 NPOT, masses, charges, 2, mass))

# mesonDecay_dune = list(mesonDecayProduction(base + 'total_efficiency_outputdecayPion1m.txt',
#                                             base + 'total_efficiency_outputdecayEta1m.txt',
#                                             base + 'total_efficiency_output2bodydecay-rho1m.txt',
#                                             base + 'total_efficiency_output2bodydecay-phi1m.txt',
#                                             base + 'total_efficiency_output2bodydecay-omega1m.txt',
#                                             base + 'total_efficiency_output2bodydecay-jsi1m.txt',
#                                             NPOT, masses, charges, 2))

# get DY production
# dy_dune = dyProduction(base+'efficiency_dy_dark.txt', base+'dy_cross_dark.txt', NPOT, charges, 2)
dy_dune = dyProduction(base+'efficiency_dy_ship_1m.txt', base+'dy_cross_ship.txt', NPOT, charges, 2)
total_dune = dy_dune + mesonDecay_dune

print("Generation completed")

figure(figsize=(22, 16), dpi=500)
plt.tick_params(axis='both', which='both', labelsize=25)

plt.contour(mass, charge, total_dune, levels=[10], colors = 'greenyellow')
plt.plot(10, 5, color='greenyellow', label='2 layers, 80 bars/layer, bkg = 20')

plt.contour(mass, charge, total_dune, levels=[21], colors = 'tomato')
plt.plot(10, 5, color='tomato', label='2 layers, 80 bars/layer, bkg = 100')

plt.contour(mass, charge, total_dune, levels=[29], colors = 'skyblue')
plt.plot(10, 5, color='skyblue', label='2 layers, 80 bars/layer, bkg = 200')

plt.contour(mass, charge, total_dune, levels=[45], colors = 'cyan')
plt.plot(10, 5, color='cyan', label='2 layers, 80 bars/layer, bkg = 500')

plt.contour(mass, charge, total_dune, levels=[63], colors='orange')
plt.plot(10, 5, color='orange', label='2 layers, 80 bars/layer, bkg = 1000')

slac = pd.read_csv('../experiment-contours/slac.csv',header=None)
colliders = pd.read_csv('../experiment-contours/colliders.csv',header=None)
bebc = pd.read_csv('../experiment-contours/bebc.csv',header=None)
charmii = pd.read_csv('../experiment-contours/charmii.csv',header=None)
mq_demo = pd.read_csv('../experiment-contours/mq_demonstrator_sort.csv',header=None)
argoneut = pd.read_csv('../experiment-contours/argoneut_sort.csv',header=None)
lsnd = pd.read_csv('../experiment-contours/LSND.csv', header=None)
lsnd.sort_values(by=lsnd.columns[0])

colliders_plot = plt.fill_between(colliders[0], colliders[1], 2, label = 'Colliders', alpha=0.5, color = 'mediumseagreen') #'green')
slac_plot = plt.fill_between(slac[0], slac[1], 2, label = 'SLAC', alpha=0.5, color = 'gold') #'orange')
bebc_plot = plt.fill_between(bebc[0], bebc[1], 2, label = 'BEBC', alpha=0.5, color = 'gray')
lsnd_plot = plt.fill_between(lsnd[0]/1000, lsnd[1], 2, label='LSND', alpha=0.5, color = 'skyblue')
charmii_plot = plt.fill_between(charmii[0], charmii[1], 2, label = 'Charm II', alpha=0.5, color = 'lightgray')
argonuet_plot = plt.fill_between(argoneut[0], argoneut[1], 2, label = 'ArgoNeuT', alpha=0.5, color = 'cornflowerblue') #color = 'mediumorchid')
mq_demo_plot = plt.fill_between(mq_demo[0], mq_demo[1], 2, label = 'MilliQan demonstrator', alpha=0.5, color = 'lightcoral')

plt.xscale('log')
plt.yscale('log')
plt.xlabel('$m_{\chi}$ [$\mathrm{GeV}/\mathrm{c}^2$]', fontsize = 35)
plt.xticks(fontsize = 35)
plt.ylabel('$\epsilon=Q/e$', fontsize = 35)
plt.yticks(fontsize = 35)
plt.xlim(0.01, 10)
plt.ylim(0.000005, 1)
plt.legend(loc='upper left', fontsize=18.5)
plt.savefig('limit-plot-LongQuest-' + str(NPOT) + '.png')
print("limit plot drawn")

# Corrected second plot for heatmap
figure(figsize=(22, 16), dpi=500)
plt.tick_params(axis='both', which='both', labelsize=25)

# Display the heatmap
cmap = LinearSegmentedColormap.from_list('white_to_blue', ['white', 'darkblue'])
plt.pcolormesh(mass, charge, total_dune, norm=LogNorm(vmin=1e0, vmax=1e16), shading='auto', cmap=cmap)
plt.colorbar()

# Adding the contour lines
plt.contour(mass, charge, total_dune, levels=[10], colors='greenyellow')
plt.plot(10, 5, color='greenyellow', label='2 layers, 80 bars/layer, bkg = 20')

plt.contour(mass, charge, total_dune, levels=[21], colors='tomato')
plt.plot(10, 5, color='tomato', label='2 layers, 80 bars/layer, bkg = 100')

plt.contour(mass, charge, total_dune, levels=[29], colors='skyblue')
plt.plot(10, 5, color='skyblue', label='2 layers, 80 bars/layer, bkg = 200')

plt.contour(mass, charge, total_dune, levels=[45], colors='cyan')
plt.plot(10, 5, color='cyan', label='2 layers, 80 bars/layer, bkg = 500')

plt.contour(mass, charge, total_dune, levels=[63], colors='orange')
plt.plot(10, 5, color='orange', label='2 layers, 80 bars/layer, bkg = 1000')

# Setting logarithmic scale for both axes
plt.xscale('log')
plt.yscale('log')

# Setting labels and limits
plt.xlabel('$m_{\chi}$ [$\mathrm{GeV}/\mathrm{c}^2$]', fontsize=35)
plt.xticks(fontsize=35)
plt.ylabel('$\epsilon=Q/e$', fontsize=35)
plt.yticks(fontsize=35)
plt.xlim(0.01, 8)
plt.ylim(0.000005, 1)

# Adding the legend
plt.legend(loc='lower right', fontsize=19.5)

# Saving the heatmap figure
plt.savefig('heatmap_longquest.png')
print("heatmap_longquest.png drawn")
