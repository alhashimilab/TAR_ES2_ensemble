'''adapt from Bei's scirpt'''
from uncertainties import ufloat, umath
import uncertainties.unumpy as unumpy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rcParams
from uncertainties.umath import *
from scipy.optimize import curve_fit
from scipy.constants import k as kB  # Boltzmann constant
from scipy.constants import h as hC  # Plancks constant
from scipy.constants import R  # gas constant
import os
import math
import sys

rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams['font.family'] = 'Arial'
rcParams['legend.title_fontsize'] = 30
sns.set_context('poster')

path = './'

R = R * 0.000239006  # gas constant in kcal/mol

bol = kB * 0.000239006  # boltzman constant in kcal/K

plank = hC * 0.000239006  # planks constant kcal*s


# vanthoff equation without error propogation
def vanhoff(x, a, b, tharm):
    result = (bol * x / plank) * np.exp(-a / (R * tharm) - (b / R) * (1 / x - 1 / tharm))
    return result


def get_parm_for_linear_k(kon_ss, koff_ss, T_ss, tharm):
    '''
    linear fitting for kon and koff using harmonic mean temperature equation
    kon_ss, koff_ss are on and off rate constants, in ufloat format
    T_ss is the experimental temperatures, in ufloat format
    tharm is the harmonic mean of T_ss

    return the raw fitting results + cofactors for kon and koff linear equations
    '''
    # set x, y and y error
    y_on_ss = unumpy.nominal_values(kon_ss)
    y_on_error_ss = unumpy.std_devs(kon_ss)

    y_off_ss = unumpy.nominal_values(koff_ss)
    y_off_error_ss = unumpy.std_devs(koff_ss)

    # linear fitting fo kon and koff
    fiton_result, fiton_err = curve_fit(lambda x, a, b: vanhoff(x, a, b, tharm), T_ss, y_on_ss, sigma=y_on_error_ss)

    fitoff_result, fitoff_err = curve_fit(lambda x, a, b: vanhoff(x, a, b, tharm), T_ss, y_off_ss, sigma=y_off_error_ss)

    # fitting results
    fiton_a, fiton_b = fiton_result

    onerr_a, onerr_b = np.sqrt(np.diag(fiton_err))

    fitoff_a, fitoff_b = fitoff_result

    offerr_a, offerr_b = np.sqrt(np.diag(fitoff_err))

    # fitting results in float
    on_a_fl = ufloat(fiton_a, onerr_a)
    on_b_fl = ufloat(fiton_b, onerr_b)

    off_a_fl = ufloat(fitoff_a, offerr_a)
    off_b_fl = ufloat(fitoff_b, offerr_b)

    return fiton_result, fitoff_result, on_a_fl, on_b_fl, off_a_fl, off_b_fl


def kon_koff_temp_plot(kon_ss, koff_ss, T_ss, fiton_result, fitoff_result, tharm, target_temp, outputname,
                       ylim1=[-5, 0], ylim2=[-5, 0], t_range=[1, 30]):
    '''
    make temperature dependent plot for kon/koff
    kon_ss, koff_ss are on and off rate constants, in ufloat format
    T_ss is the experimental temperatures, in ufloat format
    tharm is the harmonic mean of T_ss

    target temperature is the vertical line at target temperature, input unit is C
    outputname: name of the output pdf

    ylim1/2 are the ylim of plot 1 and 2
    '''
    # get the value and error for log(kon/T), log(koff/T)
    y_on_ufloat_ss = [log(z[0] / z[1]) for z in zip(kon_ss, T_ss)]
    y_off_ufloat_ss = [log(z[0] / z[1]) for z in zip(koff_ss, T_ss)]

    # set x, y and y error
    y_on_ss = unumpy.nominal_values(y_on_ufloat_ss)
    y_on_error_ss = unumpy.std_devs(y_on_ufloat_ss)

    y_off_ss = unumpy.nominal_values(y_off_ufloat_ss)
    y_off_error_ss = unumpy.std_devs(y_off_ufloat_ss)

    r_square_on = np.round(r_squred_new(T_ss, y_on_ss, vanhoff, np.append(fiton_result, tharm)), 3)
    r_square_off = np.round(r_squred_new(T_ss, y_off_ss, vanhoff, np.append(fitoff_result, tharm)), 3)

    # make the plot
    fig, ax = plt.subplots(1, 2, figsize=(15, 5))
    sim_x = np.linspace(0 + 273, 100 + 273, 10)
    ax[0].errorbar((1 / T_ss - 1 / tharm), y_on_ss, yerr=y_on_error_ss, fmt='ok', ecolor='k', ms=12, capsize=5)
    ax[0].plot((1 / sim_x - 1 / tharm), np.log(vanhoff(sim_x, *np.append(fiton_result, tharm)) / sim_x), 'k-',
               label='r$^2$={:.2f}'.format(r_square_on))
    ax[0].set_ylim(ylim1)
    ax[0].set_xlim(-0.0001, 0.0001)
    #     ax[0].set_xlim((1/(t_range[1]+273) - 1/tharm),(1/(t_range[0]+273) - 1/tharm))
    ax[0].set_ylabel('ln(k$_1$/T)')
    #     ax[0].axvline(x=1.0/(273.15 + target_temp)- 1/tharm, color='k', linestyle='dashed', linewidth=2)
    leg = ax[0].legend(frameon=False)
    for item in leg.legendHandles:
        item.set_visible(False)

    ax[1].errorbar((1 / T_ss - 1 / tharm), y_off_ss, yerr=y_off_error_ss, fmt='ok', ecolor='k', ms=12, capsize=5)
    ax[1].plot((1 / sim_x - 1 / tharm), np.log(vanhoff(sim_x, *np.append(fitoff_result, tharm)) / sim_x), 'k-',
               label='r$^2$={:.2f}'.format(r_square_off))
    ax[1].set_ylim(ylim2)
    ax[1].set_ylabel('ln(k$_{-1}$/T)')
    ax[1].set_xlim(-0.0001, 0.0001)
    #     ax[1].set_xlim((1/(t_range[1]+273) - 1/tharm),(1/(t_range[0]+273) - 1/tharm))
    #     ax[1].axvline(x=1.0/(273.15 + target_temp)- 1/tharm, color='k', linestyle='dashed', linewidth=2)
    leg = ax[1].legend(frameon=False)

    for item in leg.legendHandles:
        item.set_visible(False)

    fig.text(0.52, -0.03, '1/T-1/T$_{hm}$ (K$^{-1}$)', ha='center', fontsize=25)
    # fig.text(0.50, 0.95, 'Temperature dependence of methyl rotation rates', ha='center', fontsize = 20)
    plt.tight_layout()
    #     plt.savefig('methyl_temp_extra.pdf', dpi = 300, transparent=True)
    if outputname:
        plt.savefig(outputname + '.pdf', dpi=300, transparent=True)


# calculate R squared, for the euqation using harmonic mean of temperature
def r_squred_new(x, y, func, parms):
    residuals = y - np.log(func(x, *parms) / x)

    ss_res = np.sum(residuals ** 2)

    ss_tot = np.sum((y - np.mean(y)) ** 2)

    return 1 - ss_res / ss_tot


def get_data_for_energy_diagram(T, kex, pB, lim1=[-5, 0], lim2=[0, 4], kon_inp=None, koff_inp=None, trange=[1, 30],
                                outputname='1'):
    '''
    calculate energetics needed for energy diagram
    '''
    # harmonic mean of tempertaures
    tharm_ss = (1 / (np.sum(1 / (T + 273.15)))) * len(T)

    # manually input kon and koff
    if kon_inp:
        kon_ss, koff_ss = kon_inp, koff_inp
    else:
        kon_ss = [kex[i] * pB[i] / 100 for i in range(len(T))]

        koff_ss = [kex[i] * (1 - pB[i] / 100) for i in range(len(T))]

    # get a, b in linear equation for kon and koff
    fiton_result_ss, fitoff_result_ss, on_a_fl_ss, on_b_fl_ss, off_a_fl_ss, off_b_fl_ss = get_parm_for_linear_k(kon_ss,
                                                                                                                koff_ss,
                                                                                                                T + 273.15,
                                                                                                         tharm_ss)

    # calculate dS, dH for on rate
    dS_on = (on_b_fl_ss - on_a_fl_ss) / tharm_ss
    # calculate dS, dH for off rate
    dS_off = (off_b_fl_ss - off_a_fl_ss) / tharm_ss

    dH_on = on_b_fl_ss
    dH_off = off_b_fl_ss
    print(dH_off-dH_on)
    print(dS_on * (273.15 + 25) - dS_off * (273.15 + 25))
    dG_on_37C = -dS_on * (273.15 + 25) + dH_on

    dG_off_37C = -dS_off * (273.15 + 25) + dH_off

    Energies = [0, dG_on_37C, dG_on_37C - dG_off_37C, 0, dS_on * (273.15 + 25),
                dS_on * (273.15 + 25) - dS_off * (273.15 + 25), 0, dH_on, dH_on - dH_off]
    print(Energies)
    data = np.array(Energies)
    data = data.reshape(3, -1)
    data = np.concatenate((data[:1], data[2:3], data[-2:-1]))

    kon_koff_temp_plot(kon_ss, koff_ss, T + 273.15, fiton_result_ss, fitoff_result_ss, tharm_ss, 37, outputname,
                       lim1, lim2, trange)

    return data

# #### Free fitting U38-H3
# construct = 'wtTAR'
# pB = unumpy.uarray([0.238, 0.264, 0.359, 0.398], [0.042, 0.019, 0.042, 0.039])
# dw = unumpy.uarray([-2.73, -2.81, -2.81, -2.84], [0.01, 0.01, 0.03, 0.17])
# kex = unumpy.uarray([4417.04951, 2369.35368, 1006.48744, 472.8092], [1044.80359, 156.06627, 79.42719, 45.98247])
# output_name = '4T_wtTAR_energy_diagram(U38_free).pdf'
#
#### Fixed dw=-2.84 U38-H3
# construct = 'wtTAR'
# pB = unumpy.uarray([0.225, 0.259, 0.350, 0.397], [0.041, 0.016, 0.043, 0.041])
# dw = unumpy.uarray([-2.84, -2.84, -2.84, -2.84], [0.01, 0.01, 0.03, 0.17])
# kex = unumpy.uarray([4578.42994, 2393.42749, 1026.03225, 473.68584], [1085.77913, 157.82382, 82.21814, 46.06552])
# output_name = '4T_wtTAR_energy_diagram(U38_fix_dw).pdf'
#
# T = np.array([40, 35, 30, 25])
# data1 = get_data_for_energy_diagram(T, kex, pB, lim1 = [-6,-2], lim2 = [-2,4], outputname='wtTAR')

#### Free U38-N3
construct = 'wtTAR'
pB = unumpy.uarray([0.232, 0.246, 0.311], [0.01, 0.022, 0.026])
kex = unumpy.uarray([1653, 835.79792, 315.8486], [84,67.37823, 25.30473])
T = np.array([35, 30, 25])
data1 = get_data_for_energy_diagram(T, kex, pB, lim1 = [-6,-2], lim2 = [-2,4], outputname='wtTAR_15N_U38')
output_name = '4T_wtTAR_energy_diagram(U38_fix_dw).pdf'

#### Free fitting TARES2 U38-H3
# construct = 'TARES2'
# pB = unumpy.uarray([18.96, 14.71], [21.86,8.01])
# kex = unumpy.uarray([272.9, 160.2], [19.44, 12.26])
# output_name = '2T_TARES2_energy_diagram(U38_free).pdf'

# #### Free fitting TARES2 G54-H1
# construct = 'TARES2'
# pB = unumpy.uarray([31.444, 14.383], [29.7,6.180])
# kex = unumpy.uarray([93.01731, 117.66359], [9.04958, 13.29264])
# output_name = '2T_TARES2_energy_diagram(G54_free).pdf'

# T = np.array([30, 25])
# data1 = get_data_for_energy_diagram(T, kex, pB, lim1 = [-6,-2], lim2 = [-2,4], outputname='uucgES2_U38H3')
#

print(data1)

def energy_rank(data, marker_width=.5, color='blue', linest='solid'):
    y_data = np.repeat(data, 2)
    x_data = np.empty_like(y_data)
    x_data[0::2] = np.arange(1, len(data) + 1) - (marker_width / 2)
    x_data[1::2] = np.arange(1, len(data) + 1) + (marker_width / 2)
    lines = []
    for x in range(0, len(data) * 2, 2):
        lines.append(plt.Line2D(x_data[x:x + 2], y_data[x:x + 2], lw=8, linestyle=linest, color=color))
    lines.append(plt.Line2D(x_data, y_data, lw=3, linestyle=linest, color=color))
    
    return lines


artists1 = []
for row, color, l, st in zip(unumpy.nominal_values(data1), ('k', 'gray', 'gray'),
                             ('$\Delta$G', '$\Delta$dS', '$\Delta$dH'), ('solid', 'solid', '--')):
    artists1.extend(energy_rank(row, color=color, linest=st))

fig, ax = plt.subplots(1, 1, figsize=(9, 6))
for i in range(len(artists1)):
    ax.add_artist(artists1[i])
colors = ['k', 'gray', 'gray']
lst = ['-', '-', '--']
lines = [plt.Line2D([0], [0], color=c, linewidth=5, linestyle=l) for c, l in zip(colors, lst)]
labels = ['G', 'H', 'TS']
data_total = [data1]
ax.set_ylim([-20, 40])

ax.legend(lines, labels, frameon=False, fontsize=20)
ax.set_title(['wtTAR Proton CEST'][0])
ax.set_xlim([.5, 3.5])
ax.set_xticks([1, 2, 3])
ax.set_xticklabels(['GS', 'TS', 'ES2'])

ax.set_ylabel('Energy (kcal/mol)')

plt.savefig(os.path.join(path,output_name), dpi=300, transparent=True)
plt.show()
