import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from matplotlib import rc
from matplotlib.ticker import MaxNLocator
from matplotlib.pyplot import figure
from matplotlib.colors import LinearSegmentedColormap, LogNorm
import pandas as pd

def mesonDecayProduction(eff_file, NPOT, masses, charges, a,
                         det_col=1, N_gamma=2.5e5):
    # constants
    alpha  = 1.0/137.0
    m_pi   = 0.135  # GeV
    c_pi   = 3.98
    branch = 0.98

    # load [mchi, eff_det1, eff_det2]
    data      = np.loadtxt(eff_file, skiprows=1)
    m_file    = data[:,0]
    eff       = data[:, det_col]

    # pad or truncate so that eff.size == number of mass‐bins in your grid
    n_mass = masses.shape[1]
    if eff.size < n_mass:
        eff = np.pad(eff, (0, n_mass - eff.size), 'constant', constant_values=0)
    elif eff.size > n_mass:
        eff = eff[:n_mass]

    # now tile into a (n_charge × n_mass) array
    ageo_pi = np.tile(eff, (charges.shape[0], 1))

    # Dalitz PS integral
    def I3_integrand(z, x):
        return 2/(3*np.pi)*np.sqrt(1 - 4*x/z)*(1 - z)**3*(2*x + z)/z**2
    def I3(x):
        return quad(I3_integrand, 4*x, 1, args=(x,))[0]
    vI3 = np.vectorize(I3)

    # build sensitivity
    pi   = np.zeros_like(masses)
    mask = (masses < m_pi/2)

    xvar    = masses[mask]**2 / m_pi**2
    damping = (1 - np.exp(-N_gamma * charges[mask]**2))**a * charges[mask]**2

    pi[mask] = (
            ageo_pi[mask]
            * 2*c_pi*branch*alpha
            * vI3(xvar)
            * damping
    )

    return pi * NPOT


NPOT = 5.9e22

# masses SHIP vs Dark
# 1. load your mchi grid and charge grid
m_file   = np.loadtxt('/Users/leobailloeul/Documents/coding/decaysimulation/decay/output-data/mchi_values.txt')        # shape (M,)
charge_1d= np.logspace(-5, 1, 500)                      # shape (N,)
masses, charges = np.meshgrid(m_file, charge_1d)        # shapes both (N,M)

NPOT = 5.9e22
a    = 2

# 2. compute two sensitivity arrays
sens_det1 = mesonDecayProduction(
    '/Users/leobailloeul/Documents/coding/decaysimulation/decay/efficiencies_all_column.txt',
    NPOT, masses, charges, a, det_col=1
)
sens_det2 = mesonDecayProduction(
    '/Users/leobailloeul/Documents/coding/decaysimulation/decay/efficiencies_all_column.txt',
    NPOT, masses, charges, a, det_col=2
)

print("Generation completed")

# ——— first limit‐contour plot ———
figure(figsize=(22, 16), dpi=500)
plt.tick_params(axis='both', which='both', labelsize=25)

# LANSCE-mQ @ ER1 (Detector 1)
plt.contour(
    m_file, charge_1d, sens_det1,
    levels=[45],
    colors='greenyellow'
)
plt.plot(10, 5, color='greenyellow', label='LANSCE-mQ @ ER1')

# LANSCE-mQ @ ER2 (Detector 2)
plt.contour(
    m_file, charge_1d, sens_det2,
    levels=[45],
    colors='tomato'
)
plt.plot(10, 5, color='tomato', label='LANSCE-mQ @ ER2')

# — keep all your existing experimental contours —
slac      = pd.read_csv('../experiment-contours/slac.csv',header=None)
colliders = pd.read_csv('../experiment-contours/colliders.csv',header=None)
bebc      = pd.read_csv('../experiment-contours/bebc.csv',header=None)
charmii   = pd.read_csv('../experiment-contours/charmii.csv',header=None)
mq_demo   = pd.read_csv('../experiment-contours/mq_demonstrator_sort.csv',header=None)
argoneut  = pd.read_csv('../experiment-contours/argoneut_sort.csv',header=None)
lsnd      = pd.read_csv('../experiment-contours/LSND.csv', header=None)

plt.fill_between(colliders[0],    colliders[1],   2, label='Colliders', alpha=0.5, color='mediumseagreen')
plt.fill_between(slac[0],         slac[1],        2, label='SLAC',      alpha=0.5, color='gold')
plt.fill_between(bebc[0],         bebc[1],        2, label='BEBC',      alpha=0.5, color='gray')
plt.fill_between(lsnd[0]/1000,    lsnd[1],        2, label='LSND',      alpha=0.5, color='skyblue')
plt.fill_between(charmii[0],      charmii[1],     2, label='Charm II',  alpha=0.5, color='lightgray')
plt.fill_between(argoneut[0],     argoneut[1],    2, label='ArgoNeuT',  alpha=0.5, color='cornflowerblue')
plt.fill_between(mq_demo[0],      mq_demo[1],     2, label='MilliQan demonstrator', alpha=0.5, color='lightcoral')

plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$m_{\chi}$ [$\mathrm{GeV}/c^2$]', fontsize=35)
plt.xticks(fontsize=35)
plt.ylabel(r'$\epsilon=Q/e$',          fontsize=35)
plt.yticks(fontsize=35)
plt.xlim(0.01, 10)
plt.ylim(0.000005, 1)
plt.legend(loc='upper left', fontsize=18.5)
plt.savefig(f'limit-plot-LANSCE-mQ-{NPOT:.1e}.png')
print("limit plot drawn")


# ——— second heatmap plot ———
figure(figsize=(22, 16), dpi=500)
plt.tick_params(axis='both', which='both', labelsize=25)

cmap = LinearSegmentedColormap.from_list('white_to_blue', ['white', 'darkblue'])
plt.pcolormesh(
    m_file, charge_1d, sens_det1,  # change to sens_det2 if you want ER2
    norm=LogNorm(vmin=1e0, vmax=1e16),
    shading='auto',
    cmap=cmap
)
plt.colorbar()

# overlay the same two contours on the heatmap
plt.contour(m_file, charge_1d, sens_det1, levels=[900], colors='greenyellow')
plt.plot(10, 5, color='greenyellow', label='LANSCE-mQ @ ER1')

plt.contour(m_file, charge_1d, sens_det2, levels=[900], colors='tomato')
plt.plot(10, 5, color='tomato', label='LANSCE-mQ @ ER2')

plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$m_{\chi}$ [$\mathrm{GeV}/c^2$]', fontsize=35)
plt.xticks(fontsize=35)
plt.ylabel(r'$\epsilon=Q/e$',          fontsize=35)
plt.yticks(fontsize=35)
plt.xlim(0.01, 8)
plt.ylim(0.000005, 1)
plt.legend(loc='lower right', fontsize=19.5)
plt.savefig('heatmap_LANSCE-mQ.png')
print("heatmap drawn")