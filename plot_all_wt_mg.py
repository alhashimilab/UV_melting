import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import pandas as pd
import math
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

mpl.rcParams['axes.linewidth'] = 2.0
fs=18
R = 8.314/4184
Ct = 3 * math.pow(10, -6)

def remove_spaces(line):
    line_temp = []
    for ele in line:
        if ele != "":
            line_temp.append(ele)
    return line_temp

def read_data(filename):
    # Read the UV data file
    f = open(filename, "r")
    temperature = []
    absorbance = []
    # Now read the rest of the lines
    for line in f.readlines():
        if line[0] == 'D' or len(line) <= 2:
            break
        line = line.strip("\n")
        line = line.split("\t")
        line = remove_spaces(line)
        line[0] = float(line[0])
        line[1] = float(line[1])
        temperature.append(line[0])
        absorbance.append(line[1]) 
    return [np.array(temperature), np.array(absorbance)]

def uv_hairpin(t, mds, bds, mss, bss, delH, Tm):
    # Gas constant in kcal/mol
    R = 8.314 / 4184
    # Get the exponential factor
    exp_factor = np.exp(((1/Tm) - (1/(t+273.16))) * delH / R)
    # Get the fraction of hairpin that is intact
    f = exp_factor / (1 + exp_factor)
    # Use f to get absorbance of hairpin/unfolded species
    absorbance = (((mds*t)+bds)*f) + (((mss*t)+bss)*(1-f))
    return absorbance 


def uv_hairpin_es(t, mss, bss, mds, bds, mds2, bds2, delHconf, delSconf, delHm, delSm):
    global R, Ct

    # Get Kconf at temperature of interest
    Kconf = np.exp(np.divide((delHconf - (t*delSconf)), -R * t))
    Km = np.exp(np.divide((delHm - (t*delSm)), -R * t))

    # here, duplex refers to hairpin species
    duplex = np.divide(Ct, (1 + Kconf + Km))
    duplex_conf = np.multiply(Kconf, duplex)
    duplexm = np.multiply(Km, duplex)

    # Get the extinction coefficients
    epsilon_ss = bss + (mss*(t-273.16))
    epsilon_ds = bds + (mds*(t-273.16))
    epsilon_ds_conf = bds2 + (mds2*(t-273.16))

    # Get the absorbance
    abs_duplex = epsilon_ds*duplex
    abs_duplex_conf = epsilon_ds_conf*duplex_conf
    abs_ss = epsilon_ss*duplexm
    abs_net = (abs_duplex + abs_duplex_conf + abs_ss)/Ct

    return abs_net

def plot_uv(ax, dirname, filename, yaxislabel, ylim, yaxisticks, xlim, xaxisticks, xaxisticks_status, fit_type):
    ''' plot a 2state data + fit curve '''
    bf = pd.read_csv(dirname + filename)
    colors_plot = ['k', 'b', 'r']
    t_dummy = np.linspace(10.0, 100.0, 500)   
    counter = 0

    for dummy in range(len(bf['filename'])-2):
        # Read the data
        dataname = bf['filename'].iloc[dummy]
        [t, a] = read_data(dirname + dataname)
        # If the fit is 2-state
        if fit_type == '2state':
            #def uv_hairpin(t, mds, bds, mss, bss, delH, Tm):
            sim_a = uv_hairpin(t, bf['mds'].iloc[dummy], bf['bds'].iloc[dummy], bf['mss'].iloc[dummy], bf['bss'].iloc[dummy], bf['delH'].iloc[dummy], bf['Tm(K)'].iloc[dummy])

            # compute chi2 for comparison
            error_absorbance = np.std(a[:20])
            chi2 = np.sum(np.square((a - sim_a)/error_absorbance))
            redchi2 = chi2 / (len(a)-6.0)
            print redchi2, bf['redchi2'].iloc[dummy]

            sim_a = uv_hairpin(t_dummy, bf['mds'].iloc[dummy], bf['bds'].iloc[dummy], bf['mss'].iloc[dummy], bf['bss'].iloc[dummy], bf['delH'].iloc[dummy], bf['Tm(K)'].iloc[dummy])
            ax[counter].plot(t_dummy, sim_a, color='b', linewidth=3.0) 

            # Inset
            axins = inset_axes(ax[counter], width="75%", height="75%", bbox_to_anchor=[0.05, 0.44, 0.55, 0.55], bbox_transform=ax[counter].transAxes, borderpad=0.2)
            axins.plot(t_dummy, sim_a, color='b', linewidth=3.0)
            axins.plot(t, a, color='k', linewidth=0.0, marker='o', markersize=3) 
            axins.tick_params(axis='x', direction='out', width=2, length=3, pad=0.0, top=False)
            axins.tick_params(axis='y', direction='out', width=2, length=3, pad=0.0, right=False)
            if "wt3" in dataname:
                axins.set_ylim([0.895, 0.93])
                axins.set_yticks([0.90, 0.92])
                axins.set_yticklabels(['%3.2f'%ele for ele in [0.90, 0.92]], fontsize=fs)
                axins.set_xlim([30.0, 70.0])
                axins.set_xticks([30., 50., 70.])
                axins.set_xticklabels(['%d'%ele for ele in [30., 50., 70.]], fontsize=fs)
            else:
                axins.set_xlim([15.0, 70.0])
                axins.set_xticks([20., 40., 60.])
                axins.set_xticklabels(['%d'%ele for ele in [20., 40., 60.]], fontsize=fs)
                axins.set_ylim([0.86, 0.90])
                axins.set_yticks([0.87, 0.89])
                axins.set_yticklabels(['%3.2f'%ele for ele in [0.87, 0.89]], fontsize=fs)
 
        elif fit_type == "3state_constrain":
            sim_a = uv_hairpin_es(t+273.16, bf['mss'].iloc[dummy], bf['bss'].iloc[dummy], bf['mds'].iloc[dummy], bf['bds'].iloc[dummy], bf['mds2'].iloc[dummy], bf['bds2'].iloc[dummy], bf['delHconf'].iloc[dummy], bf['delSconf'].iloc[dummy], -1.0*bf['delHm'].iloc[dummy], -1.0 * bf['delSm'].iloc[dummy])

            # compute chi2 for comparison
            error_absorbance = np.std(a[:20])
            chi2 = np.sum(np.square((a - sim_a)/error_absorbance))
            redchi2 = chi2 / (len(a)-9.0)
            print redchi2, bf['redchi2'].iloc[dummy]

            sim_a = uv_hairpin_es(t_dummy+273.16, bf['mss'].iloc[dummy], bf['bss'].iloc[dummy], bf['mds'].iloc[dummy], bf['bds'].iloc[dummy], bf['mds2'].iloc[dummy], bf['bds2'].iloc[dummy], bf['delHconf'].iloc[dummy], bf['delSconf'].iloc[dummy], -1.0*bf['delHm'].iloc[dummy], -1.0 * bf['delSm'].iloc[dummy])
            ax[counter].plot(t_dummy, sim_a, color='b', linewidth=3.0) 
 
            # inset
            axins = inset_axes(ax[counter], width="75%", height="75%", bbox_to_anchor=[0.05, 0.44, 0.55, 0.55], bbox_transform=ax[counter].transAxes, borderpad=0.2)
            axins.plot(t_dummy, sim_a, color='b', linewidth=3.0)
            axins.plot(t, a, color='k', linewidth=0.0, marker='o', markersize=3) 
            axins.tick_params(axis='x', direction='out', width=2, length=3, pad=0.0, top=False)
            axins.tick_params(axis='y', direction='out', width=2, length=3, pad=0.0, right=False)
            if "wt3" in dataname:
                axins.set_ylim([0.895, 0.93])
                axins.set_yticks([0.90, 0.92])
                axins.set_yticklabels(['%3.2f'%ele for ele in [0.90, 0.92]], fontsize=fs)
                axins.set_xlim([30.0, 70.0])
                axins.set_xticks([30., 50., 70.])
                axins.set_xticklabels(['%d'%ele for ele in [30., 50., 70.]], fontsize=fs)
            else:
                axins.set_xlim([15.0, 70.0])
                axins.set_xticks([20., 40., 60.])
                axins.set_xticklabels(['%d'%ele for ele in [20., 40., 60.]], fontsize=fs)
                axins.set_ylim([0.86, 0.90])
                axins.set_yticks([0.87, 0.89])
                axins.set_yticklabels(['%3.2f'%ele for ele in [0.87, 0.89]], fontsize=fs)

        elif fit_type == "3state_fixtherm":
            sim_a = uv_hairpin_es(t+273.16, bf['mss'].iloc[dummy], bf['bss'].iloc[dummy], bf['mds'].iloc[dummy], bf['bds'].iloc[dummy], bf['mds2'].iloc[dummy], bf['bds2'].iloc[dummy], bf['delHconf'].iloc[dummy], bf['delSconf'].iloc[dummy], -1.0*bf['delHm'].iloc[dummy], -1.0 * bf['delSm'].iloc[dummy])

            # compute chi2 for comparison
            error_absorbance = np.std(a[:20])
            chi2 = np.sum(np.square((a - sim_a)/error_absorbance))
            redchi2 = chi2 / (len(a)-8.0)
            print redchi2, bf['redchi2'].iloc[dummy]

            sim_a = uv_hairpin_es(t_dummy+273.16, bf['mss'].iloc[dummy], bf['bss'].iloc[dummy], bf['mds'].iloc[dummy], bf['bds'].iloc[dummy], bf['mds2'].iloc[dummy], bf['bds2'].iloc[dummy], bf['delHconf'].iloc[dummy], bf['delSconf'].iloc[dummy], -1.0*bf['delHm'].iloc[dummy], -1.0 * bf['delSm'].iloc[dummy])
            ax[counter].plot(t_dummy, sim_a, color='b', linewidth=3.0) 

            # inset
            axins = inset_axes(ax[counter], width="75%", height="75%", bbox_to_anchor=[0.05, 0.44, 0.55, 0.55], bbox_transform=ax[counter].transAxes, borderpad=0.2)
            axins.plot(t_dummy, sim_a, color='b', linewidth=3.0)
            axins.plot(t, a, color='k', linewidth=0.0, marker='o', markersize=3) 
            axins.set_xlim([15.0, 60.0])
            axins.set_ylim([0.86, 0.90])
            axins.tick_params(axis='x', direction='out', width=2, length=3, pad=0.0, top=False)
            axins.tick_params(axis='y', direction='out', width=2, length=3, pad=0.0, right=False)
            axins.set_yticks([0.87, 0.89])
            axins.set_yticklabels(['%3.2f'%ele for ele in [0.87, 0.89]], fontsize=fs)
            axins.set_xticks([20., 40., 60.])
            axins.set_xticklabels(['%d'%ele for ele in [20., 40., 60.]], fontsize=fs)

        elif fit_type == "3state_free":
            sim_a = uv_hairpin_es(t+273.16, bf['mss'].iloc[dummy], bf['bss'].iloc[dummy], bf['mds'].iloc[dummy], bf['bds'].iloc[dummy], bf['mds2'].iloc[dummy], bf['bds2'].iloc[dummy], bf['delHconf'].iloc[dummy], bf['delSconf'].iloc[dummy], -1.0*bf['delHm'].iloc[dummy], -1.0 * bf['delSm'].iloc[dummy])

            # compute chi2 for comparison
            error_absorbance = np.std(a[:20])
            chi2 = np.sum(np.square((a - sim_a)/error_absorbance))
            redchi2 = chi2 / (len(a)-10.0)
            print redchi2, bf['redchi2'].iloc[dummy]

            sim_a = uv_hairpin_es(t_dummy+273.16, bf['mss'].iloc[dummy], bf['bss'].iloc[dummy], bf['mds'].iloc[dummy], bf['bds'].iloc[dummy], bf['mds2'].iloc[dummy], bf['bds2'].iloc[dummy], bf['delHconf'].iloc[dummy], bf['delSconf'].iloc[dummy], -1.0*bf['delHm'].iloc[dummy], -1.0 * bf['delSm'].iloc[dummy])
            ax[counter].plot(t_dummy, sim_a, color='b', linewidth=3.0) 

            # inset
            axins = inset_axes(ax[counter], width="75%", height="75%", bbox_to_anchor=[0.05, 0.44, 0.55, 0.55], bbox_transform=ax[counter].transAxes, borderpad=0.2)
            axins.plot(t_dummy, sim_a, color='b', linewidth=3.0)
            axins.plot(t, a, color='k', linewidth=0.0, marker='o', markersize=3) 
            if "wt3" in dataname:
                axins.set_ylim([0.895, 0.93])
                axins.set_yticks([0.90, 0.92])
                axins.set_yticklabels(['%3.2f'%ele for ele in [0.90, 0.92]], fontsize=fs)
                axins.set_xlim([30.0, 70.0])
                axins.set_xticks([30., 50., 70.])
                axins.set_xticklabels(['%d'%ele for ele in [30., 50., 70.]], fontsize=fs)
            else:
                axins.set_xlim([15.0, 70.0])
                axins.set_xticks([20., 40., 60.])
                axins.set_xticklabels(['%d'%ele for ele in [20., 40., 60.]], fontsize=fs)
                axins.set_ylim([0.86, 0.90])
                axins.set_yticks([0.87, 0.89])
                axins.set_yticklabels(['%3.2f'%ele for ele in [0.87, 0.89]], fontsize=fs)
            axins.tick_params(axis='x', direction='out', width=2, length=3, pad=0.0, top=False)
            axins.tick_params(axis='y', direction='out', width=2, length=3, pad=0.0, right=False)

        # data points
        ax[counter].plot(t, a, color='k', linewidth=0.0, marker='o', markersize=3) 
        if counter == 0:
            ax[counter].set_ylabel(yaxislabel, fontsize=fs)

        ax[counter].tick_params(axis='x', direction='out', width=2, length=5, top=False)
        ax[counter].tick_params(axis='y', direction='out', width=2, length=5, right=False)
        ax[counter].set_title('r' + '$\chi^{2}$' + ' ' + '%.2f'%redchi2, fontsize=fs)
        ax[counter].set_yticks(yaxisticks)
        #if counter == 0:
        ax[counter].set_yticklabels(['%3.2f'%ele for ele in yaxisticks], fontsize=fs)
        #else:
        #    ax[counter].set_yticklabels([])
        ax[counter].set_xticks(xaxisticks)
        if xaxisticks_status == 1:
            ax[counter].set_xticklabels(['%d'%ele for ele in xaxisticks], fontsize=fs)
        else:
            ax[counter].set_xticklabels([])
        if "wt3" in dataname:
            ax[counter].set_xlim([30.0, 100.0])
            ax[counter].set_ylim([0.89, 1.05])
        else:
            print "here"
            ax[counter].set_xlim(xlim)
            ax[counter].set_ylim(ylim)
 
        counter = counter + 1
    print
   
fig, axarr = plt.subplots(4, 3, figsize=(16, 16.))
#fig, axarr = plt.subplots(5, 3)

plot_uv(axarr[0, :], 'wt_mg/2state_final/', '2state.csv', 'TAR + 3mM Mg' + '$^{2+}$' + '\n2-state\nA' + '$_{260}$', [0.86, 1.05], np.arange(0.90, 1.01, 0.05), [15.0, 100.0], np.arange(20.0, 101.0, 20), 0, '2state')
plot_uv(axarr[1, :], 'wt_mg/3state_free_final/', '3state_free_final.csv', 'TAR + 3mM Mg' + '$^{2+}$' + '\n3-state free\nA' + '$_{260}$', [0.86, 1.05], np.arange(0.90, 1.01, 0.05), [15.0, 100.0], np.arange(20.0, 101.0, 20), 0, '3state_free')
plot_uv(axarr[2, :], 'wt_mg/3state_constrain/', '3state_constrain_7pct.csv', 'TAR + 3mM Mg' + '$^{2+}$' + '\n3-state p' + '$_{B}$' + ' 7%\nA' + '$_{260}$', [0.86, 1.05], np.arange(0.90, 1.01, 0.05), [15.0, 100.0], np.arange(20.0, 101.0, 20), 0, '3state_constrain')
plot_uv(axarr[3, :], 'wt_mg/3state_constrain/', '3state_constrain_14pct.csv', 'TAR + 3mM Mg' + '$^{2+}$' + '\n3-state p' + '$_{B}$' + ' 14%\nA' + '$_{260}$', [0.86, 1.05], np.arange(0.90, 1.01, 0.05), [15.0, 100.0], np.arange(20.0, 101.0, 20), 1, '3state_constrain')

for dummy in range(3):
    axarr[3, dummy].set_xlabel('Temperature (' + '$^\circ$' + 'C)', fontsize=fs)
#plt.show()
plt.savefig('wt_mg.pdf')

