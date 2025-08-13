import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from matplotlib import rc
from matplotlib.ticker import MaxNLocator
from matplotlib.pyplot import figure
from matplotlib.colors import LinearSegmentedColormap, LogNorm
import pandas as pd

background = [20, 100, 200, 500, 1000]
levels = [10, 21, 29, 45, 63]
i = 0
#New Plot bkg = 20, 100, 200, 500, 1000]
for i in range(0, len(levels)):
    bkg = background[i]
    figure(figsize=(22, 16), dpi=500)
    plt.tick_params(axis='both', which='both', labelsize=25)
    plt.title('', fontsize=25)

    slac = pd.read_csv('slac.csv',header=None)
    colliders = pd.read_csv('colliders.csv',header=None)
    bebc = pd.read_csv('bebc.csv',header=None)
    charmii = pd.read_csv('charmii.csv',header=None)
    mq_demo = pd.read_csv('mq_demonstrator_sort.csv',header=None)
    argoneut = pd.read_csv('argoneut_sort.csv',header=None)
    lsnd = pd.read_csv('LSND.csv', header=None)
    sensei = pd.read_csv('SENSEI_MINOS.csv', header=None)
    FORMOSA = pd.read_csv('FORMOSA.csv', header=None)
    SUBMET = pd.read_csv('SUBMET.csv', header=None)
    LANCSE = pd.read_csv('LANSCE-mQ_ER1.csv', header=None)
    FERMINI = pd.read_csv('FERMINI.csv', header=None)
    LongQuest = pd.read_csv('LongQuest.csv', header=None)
    CMS = pd.read_csv('CMS.csv', header=None)

    lsnd.sort_values(by=lsnd.columns[0])

    # Plotting FORMOSA with a deep teal color
    plt.plot(FORMOSA[0], FORMOSA[1], linewidth=4, color='#008080', label='FORMOSA')  # Deep Teal

    # Plotting LANCSE with a bright orange-red color
    plt.plot(LANCSE[0], LANCSE[1], linewidth=4, color='#FF5733', label='LANSCE-mQ')  # Bright Orange-Red

    # Plotting SUBMET with a royal blue color
    plt.plot(SUBMET[0], SUBMET[1], linewidth=4, color='#005EB8', label='SUBMET')  # Royal Blue

    # Plotting FERMINI with a vibrant purple color
    plt.plot(FERMINI[0], FERMINI[1], linewidth=4, color='#9B59B6', label='FerMINI')  # Purple

    # Plotting LongQuest with a bright cyan color
    plt.plot(LongQuest[0], LongQuest[1], linewidth=4, color='#00BCD4', label='LongQuest')  # Bright Cyan

    colliders_plot = plt.fill_between(colliders[0], colliders[1], 2, label = 'Colliders', alpha=0.5, color = 'mediumseagreen') #'green')
    slac_plot = plt.fill_between(slac[0], slac[1], 2, label = 'SLAC', alpha=0.5, color = 'gold') #'orange')
    bebc_plot = plt.fill_between(bebc[0], bebc[1], 2, label = 'BEBC', alpha=0.5, color = 'gray')
    sensei_plot = plt.fill_between(sensei[0]/1000, sensei[1], 2, label = 'SENSEI@MINOS', alpha=0.5, color = 'darksalmon')
    lsnd_plot = plt.fill_between(lsnd[0]/1000, lsnd[1], 2, label='LSND', alpha=0.5, color = 'skyblue')
    charmii_plot = plt.fill_between(charmii[0], charmii[1], 2, label = 'Charm II', alpha=0.5, color = 'lightgray')
    argonuet_plot = plt.fill_between(argoneut[0], argoneut[1], 2, label = 'ArgoNeuT', alpha=0.5, color = 'cornflowerblue') #color = 'mediumorchid')
    mq_demo_plot = plt.fill_between(mq_demo[0], mq_demo[1], 2, label = 'MilliQan demonstrator', alpha=0.5, color = 'lightcoral')
    CMS_plot = plt.fill_between(CMS[0], CMS[1], 2, label = 'CMS', alpha=0.5, color = 'darkgreen')

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$m_{\chi}$ [$\mathrm{GeV}/\mathrm{c}^2$]', fontsize = 35)
    plt.xticks(fontsize = 35)
    plt.ylabel('$\epsilon=Q/e$', fontsize = 35)
    plt.yticks(fontsize = 35)
    plt.xlim(0.02, 100)
    plt.ylim(0.00005, 0.5)
    plt.legend(loc='upper left', fontsize=18.5)
    plt.savefig('limit-plot-LongQuest-bkg=' + str(bkg) + '.png')
    print("limit plot drawn for bkg = " + str(bkg))
    i = i+1