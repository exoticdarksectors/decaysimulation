import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.integrate import quad
import os, glob

def dyProduction(fileNamecross, filenamedy, NPOT, totalCrossSection):

    # totalCrossSection = 300e-3

    # Import Data from efficiency files
    ageo_datacross = np.loadtxt(fileNamecross)
    ageo_datady = np.loadtxt(filenamedy)


    cross_dy = ageo_datacross[:,1] # pico barn
    ageo_dy = ageo_datady[:]

    return NPOT * (1.42e-12) * cross_dy * ageo_dy

if __name__ == '__main__':
    NPOT_dark = 1e21
    base = '/Users/leobailloeul/Documents/coding/decaysimulation-archive/plotting/sensitivity-plot/'
    dy_ship = dyProduction(base+'ageo_dy_test.txt', base+'dy_cross_test.txt', NPOT_dark, 300e-3)
    mass_ship = np.loadtxt('/Users/leobailloeul/Documents/coding/decaysimulation/plotting/sensitivity-plot/mship_values.txt')
    fileNamecross = base+'ageo_dy_test.txt'
    ageo_datacross = np.loadtxt(fileNamecross)
    mass_ship2 = mass_ship
    mass_ship = ageo_datacross[:, 0]
    fig, ax = plt.subplots(figsize=(8.09, 5))
    dy_ship2 = dyProduction(base+'ageo_dy_dune.txt', base+'dy_cross_dune.txt', NPOT_dark, 300e-3)


    # axes + scales
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.margins(x=0, y=0)
    ax.set_xlabel(r'$m_{\chi}\, [\mathrm{GeV}/c^{2}]$')
    ax.set_ylabel(r'$N_{\chi}/\epsilon^{2}$')

    # nice ticks (powers of 10 with sparse labels)
    x_ticks = [10**i for i in range(-2, 2)]
    x_tick_labels = [r'$10^{-2}$', r'$10^{-1}$', r'$1$', r'$10^{1}$']
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_tick_labels)

    y_ticks = [10**i for i in range(0, 20)]
    y_tick_labels = [r'$1$'] + [ (rf'$10^{{{i}}}$' if i % 2 == 0 else ' ') for i in range(1, 20) ]
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_tick_labels)
    ax.set_ylim(1e0, 1e16)



    # thin grid (optional, helps readability)
    ax.grid(True, which='both', linewidth=0.5, alpha=0.25)

    ax.plot(mass_ship, dy_ship, lw=3, label = r'$q\overline{q}\to\gamma^{*}\to\chi\overline{\chi}$', marker='o')
    ax.plot(mass_ship2, dy_ship2, lw=3, label = r'$q\overline{q}\to\gamma^{*}\to\chi\overline{\chi}$', color='red')

    # legend styled like your example
    leg = ax.legend(loc='lower left', ncol=2, frameon=True, fancybox=True, framealpha=0.85)
    for line in leg.get_lines():
        line.set_linewidth(2.0)

    # top-right annotation (two lines)
    textbox_props = dict(boxstyle='round', facecolor='white', edgecolor='none', alpha=0.8)
    ax.text(0.62, 0.93,
            r'$1$ Year at DUNE ($10^{21}$ POT)' '\n'
            r'$574$ m, $0.56$ m radius cylindrical detector',
            transform=ax.transAxes, ha='left', va='center', bbox=textbox_props)

    # save full-page-ish
    fig.tight_layout()
    fig.savefig('MCP-production-dy-DUNE.pdf', dpi=300, bbox_inches='tight')
    print(dy_ship2.shape)
    print(mass_ship2[154], dy_ship2[154])
    print(mass_ship[2], dy_ship[2])
    # for i in range(0,199):
        # print(i, mass_ship2[i], dy_ship2[i])
    