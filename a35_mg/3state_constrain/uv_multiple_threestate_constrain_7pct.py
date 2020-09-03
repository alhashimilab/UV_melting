# Curve fitting for melting curves
import matplotlib.pyplot as plt
import numpy as np
import sys
import math
from scipy.optimize import curve_fit
import os

fl = 18
tl = 14

def uv_duplex(t, mds, bds, mss, bss, delH, Tm):
    global R
    # Get the exponential factor
    exp_factor = np.exp(((1/Tm) - (1/t)) * delH / R)
    # Get the fraction of double stranded species for a duplex
    f = ((1 + (4*exp_factor)) - np.sqrt(1 + (8*exp_factor)))/(4*exp_factor)
    # Use f to get absorbance of duplex/ss mixture
    absorbance = (((mds*t)+bds)*f) + (((mss*t)+bss)*(1-f))
    return absorbance

def uv_hairpin(t, mds, bds, mss, bss, delH, Tm):
    global R
    # Get the exponential factor
    exp_factor = np.exp(((1/Tm) - (1/t)) * delH / R)
    # Get the fraction of hairpin that is intact
    f = exp_factor / (1 + exp_factor)
    # Use f to get absorbance of hairpin/unfolded species
    absorbance = (((mds*t)+bds)*f) + (((mss*t)+bss)*(1-f))
    return absorbance

def uv_duplex_es(t, mss, bss, mds, bds, mds2, bds2, delHconf, delHm, delSm):
    global R, delG_constrain, T_constrain
    print R, Ct
    print mss, bss, mds, bds, mds2, bds2, delHconf, delSconf, delHm, delSm

    # Get delSconf
    # delG_conf = delH_conf - T * delS_conf
    delSconf = (delHconf - delG_constrain) / (273.16 + T_constrain)

    # Get Kconf at temperature of interest
    Kconf = np.exp(np.divide((delHconf - (t*delSconf)), -R * t))
    Km = np.exp(np.divide((delHm - (t*delSm)), -R * t))
 
    # Get AB concentration
    # The negative root is the right answer
    alpha = 1 + Kconf
    discriminant = np.sqrt(np.multiply(Km, Km) + (np.multiply(alpha, Km) * 4 * Ct))
    pre_term = (Km + (2 * Ct * alpha))
    denominator = 2 * np.multiply(alpha, alpha)

    # Choose the negative root as it is the only one consistent with 
    # what one would observe
    duplex = np.divide((pre_term - discriminant), denominator)
    duplex_conf = np.multiply(Kconf, duplex)
    duplexm = np.sqrt(np.multiply(Km, duplex))
 
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

def uv_hairpin_es(t, mss, bss, mds, bds, mds2, bds2, delHconf, delHm, delSm):
    global R, Ct
    global delG_constrain, T_constrain

    # Get delSconf
    # delG_conf = delH_conf - T * delS_conf
    delSconf = (delHconf - delG_constrain) / (273.16 + T_constrain)

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

def fit_data(filename, temperature, absorbance):
    global K
    global mode, Ct
    global R
    t_dummy = np.linspace(0,100,500) + 273.16
    fig, ax = plt.subplots(1, 3)

    if mode == 'hairpin': 
        print "Hairpin"
        #def uv_hairpin_es(t, mss, bss, mds, bds, mds2, bds2, delHconf, delHm, delSm):
        #p0 = np.array([0.000417, 0.984, 0.000407, 0.869, 0.000407, 0.869, 3.54, 75.1, 0.21398])
        p0 = np.array([0.000417, 0.984, 0.000407, 0.869, 0.000407, 0.869, 0.9, 75.1, 0.21398])
        #p0 = np.array([0.000417, 0.984, 0.000407, 0.869, 0.000407, 0.869, -3.54, 75.1, 0.21398])
        #p0 = np.array([0.000417, 0.984, 0.000407, 0.869, 0.000407, 0.869, -27.13, 75.1, 0.21398])
        #p0 = np.array([0.000417, 0.984, 0.000407, 0.869, 0.000407, 0.869, -0.9, 75.1, 0.21398])
        #p0 = np.array([0.000417, 0.984, 0.000407, 0.869, 0.000407, 0.869, 27.13, 75.1, 0.21398])
        [popt, pcov] = curve_fit(uv_hairpin_es, temperature, absorbance, p0, maxfev=150000)

        # Statistics
        N = float(len(absorbance))
        ideal_absorbance = uv_hairpin_es(temperature, *popt)
        error_absorbance = np.std(absorbance[:20])
        chi2 = np.sum(np.square((ideal_absorbance - absorbance)/error_absorbance))
        redchi2 = chi2 / (N - K)
        residual = np.sum(np.square(ideal_absorbance - absorbance))
 
        # Compute AIC and BIC values
        # First AIC
        if (N / K) >= 40:
            aic = (N * math.log(residual/N)) + (2 * K)
        else:
            aic = (N * math.log(residual/N)) + (2 * K) + ((2 * K * (K+1))/(N-K-1))
        bic = (N * math.log(residual/N)) + (K*math.log(N))

        # Unpack popt array
        [mss, bss, mds, bds, mds2, bds2, delHconf, delHm, delSm] = popt
        # Get delSconf
        # delG_conf = delH_conf - T * delS_conf
        delSconf = (delHconf - delG_constrain) / (273.16 + T_constrain)

        # Get concentrations from popt
        # Get Kconf at temperature of interest
        Kconf = np.exp(np.divide((delHconf - (t_dummy*delSconf)), -R * t_dummy))
        Km = np.exp(np.divide((delHm - (t_dummy*delSm)), -R * t_dummy))
 
        # here, duplex refers to hairpin species
        hairpin = np.divide(Ct, (1 + Kconf + Km))
        hairpin_conf = np.multiply(Kconf, hairpin)
        hairpinm = np.multiply(Km, hairpin)

        # Get the extinction coefficients
        epsilon_ss = bss + (mss*(t_dummy-273.16))
        epsilon_ds = bds + (mds*(t_dummy-273.16))
        epsilon_ds_conf = bds2 + (mds2*(t_dummy-273.16))
        min_epsilon = min(list(epsilon_ss) + list(epsilon_ds) + list(epsilon_ds_conf))
        max_epsilon = max(list(epsilon_ss) + list(epsilon_ds) + list(epsilon_ds_conf))
        min_epsilon = min_epsilon - ((min_epsilon%0.050))
        max_epsilon = max_epsilon + (0.050 - (min_epsilon%0.050))
        min_absorbance = min(absorbance)
        max_absorbance = max(absorbance)
        min_absorbance = min_absorbance - ((min_absorbance%0.05))
        max_absorbance = max_absorbance + (0.05 - (min_absorbance%0.05))

        # Get the absorbance
        abs_hairpin = epsilon_ds*hairpin
        abs_hairpin_conf = epsilon_ds_conf*hairpin_conf
        abs_ss = epsilon_ss*hairpinm
        abs_net = (abs_hairpin + abs_hairpin_conf + abs_ss)/Ct
       
        # Plot concentrations
        ax[0].set_ylabel("Population", fontsize=fl)
        ax[0].set_xlabel("Temperature (" + "$^\circ$" + "C)", fontsize=fl) 
        ax[0].plot(t_dummy-273.16, hairpin/Ct, color='k', label='A')
        ax[0].plot(t_dummy-273.16, hairpin_conf/Ct, color='r', label="A'")
        ax[0].plot(t_dummy-273.16, hairpinm/Ct, color='b', label="Amelt")
        ax[0].set_yticks(np.arange(0.0, 1.1, 0.2))
        ax[0].set_yticklabels(["%0.1f"%ele for ele in np.arange(0.0, 1.1, 0.2)], fontsize=tl)
        ax[0].set_xticks(np.arange(0.0, 100.1, 20))
        ax[0].set_xticklabels(["%2.0f"%ele for ele in np.arange(0.0, 100.1, 20)], fontsize=tl)
        ax[0].set_xlim([t_dummy[0]-273.16, t_dummy[-1]-273.16])
        ax[0].set_ylim([0.0, 1.0])
        ax[0].plot([T_constrain, T_constrain], [0., 1.], linewidth=2, color='k', linestyle='--')
        ax[0].plot([t_dummy[0]-273.16, t_dummy[-1]-273.16], [pB_constrain, pB_constrain], linewidth=2, color='r', linestyle='--')
        ax[0].legend(fontsize=fl)
 
        ax[1].set_ylabel("Extinction coefficient (" + "$\epsilon$" + ")", fontsize=fl)
        ax[1].plot(t_dummy-273.16, epsilon_ds, color='k', label='A')
        ax[1].plot(t_dummy-273.16, epsilon_ds_conf, color='r', label="A'")
        ax[1].plot(t_dummy-273.16, epsilon_ss, color='b', label="Amelt")
        ax[1].set_xticks(np.arange(0.0, 100.1, 20))
        ax[1].set_xticklabels(["%2.0f"%ele for ele in np.arange(0.0, 100.1, 20)], fontsize=tl)
        ax[1].set_xlim([t_dummy[0]-273.16, t_dummy[-1]-273.16])
        ax[1].set_yticks(np.linspace(min_epsilon, max_epsilon, 10))    
        ax[1].set_yticklabels(["%1.2f"%ele for ele in np.linspace(min_epsilon, max_epsilon, 10)], fontsize=tl)    
        ax[1].set_xlabel("Temperature (" + "$^\circ$" + "C)", fontsize=fl) 
        ax[1].legend(fontsize=fl)

        ax[2].plot(temperature-273.16, absorbance, color='k', linewidth=0.0, marker='o', markersize=1, label='data')
        ax[2].plot(t_dummy-273.16, abs_net, color='b', label='fit')
        ax[2].set_xticks(np.arange(0.0, 100.1, 20))
        ax[2].set_xticklabels(["%2.0f"%ele for ele in np.arange(0.0, 100.1, 20)], fontsize=tl)
        ax[2].set_xlim([t_dummy[0]-273.16, t_dummy[-1]-273.16])
        ax[2].set_title('$\chi^{2}$' + ' ' + str('%.3f'%redchi2) + '\n')
        ax[2].set_ylabel("Absorbance", fontsize=fl)
        ax[2].set_xlabel("Temperature (" + "$^\circ$" + "C)", fontsize=fl) 
        ax[2].legend(fontsize=fl)
        fig.suptitle(filename)
        plt.show()

        # Flip the sign of delHm and delSm to correspond to annealing
        # Compute delG
        popt[7] = popt[7] * -1.0
        popt[8] = popt[8] * -1.0
        [mss, bss, mds, bds, mds2, bds2, delHconf, delHm, delSm] = popt
        # Get delSconf
        # delG_conf = delH_conf - T * delS_conf
        delSconf = (delHconf - delG_constrain) / (273.16 + T_constrain)

        delGm_25 = delHm - ((273.16+25)*delSm)
        delGm_37 = delHm - ((273.16+37)*delSm)
        delGconf_25 = delHconf - ((273.16+25)*delSconf)
        delGconf_37 = delHconf - ((273.16+37)*delSconf)
        #print "delGm_25C ", delGm_25
        #print "delGm_37C ", delGm_37

        # Now, if delGconf_25C is netagive, we have an issue, 
        # this happens because the GS and ES are swapped
        if delGconf_25 < 0.0:
            print "*** SWAPPING ***"
            delHm_new = delHm + delHconf
            delSm_new = delSm + delSconf
            delHconf_new = -1.0  * delHconf
            delSconf_new = -1.0 * delSconf
 
            delHm = delHm_new
            delSm = delSm_new
            delHconf = delHconf_new 
            delSconf = delSconf_new

            mds_new = mds
            bds_new = bds
            mds = mds2
            bds = bds2
            mds2 = mds_new
            bds2 = bds_new

            popt = [mss, bss, mds, bds, mds2, bds2, delHconf, delHm, delSm]

        delGm_25 = delHm - ((273.16+25)*delSm)
        delGm_37 = delHm - ((273.16+37)*delSm)
        delGconf_25 = delHconf - ((273.16+25)*delSconf)
        delGconf_37 = delHconf - ((273.16+37)*delSconf)
        print "delGm_25C ", delGm_25
        print "delGm_37C ", delGm_37

        # Create results array
        print_result(popt, pcov, ['mss', 'bss', 'mds', 'bds', 'mds2', 'bds2', 'delHconf', 'delHm', 'delSm'])
        labels = ['mss', 'bss', 'mds', 'bds', 'mds2', 'bds2', 'delHconf', 'delSconf', 'pB_constrain', 'T_constrain', 'delHm', 'delSm', 'delGconf_25C', 'delGconf_37C', 'delGm_25C', 'delGm_37C', 'Ct', 'redchi2', 'aic', 'bic']
        result = np.array(list(popt[:7]) + [delSconf, pB_constrain, T_constrain] + list(popt[7:]) + [delGconf_25, delGconf_37, delGm_25, delGm_37, Ct, redchi2, aic, bic])

    elif mode == 'duplex':
        print "Duplex"
        #def uv_duplex_es(t, mss, bss, mds, bds, mds2, bds2, delHconf, delHm, delSm):
        p0 = np.array([0.0004288, 0.74594, 0.000548, 0.6555, 0.000348, 0.677, -27.1, 83.1, 0.243])
        [popt, pcov] = curve_fit(uv_duplex_es, temperature, absorbance, p0)

        # Statistics
        N = float(len(absorbance))
        ideal_absorbance = uv_duplex_es(temperature, *popt)
        error_absorbance = np.std(absorbance[:20])
        chi2 = np.sum(np.square((ideal_absorbance - absorbance)/error_absorbance))
        redchi2 = chi2 / (N - K)
        residual = np.sum(np.square(ideal_absorbance - absorbance))
 
        # Compute AIC and BIC values
        # First AIC
        if (N / K) >= 40:
            aic = (N * math.log(residual/N)) + (2 * K)
        else:
            aic = (N * math.log(residual/N)) + (2 * K) + ((2 * K * (K+1))/(N-K-1))
        bic = (N * math.log(residual/N)) + (K*math.log(N))

        # Unpack popt array
        [mss, bss, mds, bds, mds2, bds2, delHconf, delHm, delSm] = popt
        # Get delSconf
        # delG_conf = delH_conf - T * delS_conf
        delSconf = (delHconf - delG_constrain) / (273.16 + T_constrain)

        # Get Kconf at temperature of interest
        Kconf = np.exp(np.divide((delHconf - (t_dummy*delSconf)), -R * t_dummy))
        Km = np.exp(np.divide((delHm - (t_dummy*delSm)), -R * t_dummy))
 
        # Get AB concentration
        # The negative root is the right answer
        alpha = 1 + Kconf
        discriminant = np.sqrt(np.multiply(Km, Km) + (np.multiply(alpha, Km) * 4 * Ct))
        pre_term = (Km + (2 * Ct * alpha))
        denominator = 2 * np.multiply(alpha, alpha)

        # Choose the negative root as it is the only one consistent with 
        # what one would observe
        duplex = np.divide((pre_term - discriminant), denominator)
        duplex_conf = np.multiply(Kconf, duplex)
        duplexm = np.sqrt(np.multiply(Km, duplex))
  
        # Get the extinction coefficients
        epsilon_ss = bss + (mss*(t_dummy-273.16))
        epsilon_ds = bds + (mds*(t_dummy-273.16))
        epsilon_ds_conf = bds2 + (mds2*(t_dummy-273.16))
        min_epsilon = min(list(epsilon_ss) + list(epsilon_ds) + list(epsilon_ds_conf))
        max_epsilon = max(list(epsilon_ss) + list(epsilon_ds) + list(epsilon_ds_conf))
        min_epsilon = min_epsilon - ((min_epsilon%0.050))
        max_epsilon = max_epsilon + (0.050 - (min_epsilon%0.050))

        # Get the absorbance
        abs_duplex = epsilon_ds*duplex
        abs_duplex_conf = epsilon_ds_conf*duplex_conf 
        abs_ss = epsilon_ss*duplexm
        abs_net = (abs_duplex + abs_duplex_conf + abs_ss)/Ct
        min_absorbance = min(list(absorbance))
        max_absorbance = max(list(absorbance))
        min_absorbance = min_absorbance - ((min_absorbance%0.05))
        max_absorbance = max_absorbance + (0.05 - (min_absorbance%0.05))

        ax[0].set_ylabel("Population", fontsize=fl)
        ax[0].set_xlabel("Temperature (" + "$^\circ$" + "C)", fontsize=fl) 
        ax[0].plot(t_dummy-273.16, duplex/Ct, color='k', label='AB')
        ax[0].plot(t_dummy-273.16, duplex_conf/Ct, color='r', label="AB'")
        ax[0].plot(t_dummy-273.16, duplexm/Ct, color='b', label="A")
        ax[0].set_yticks(np.arange(0.0, 1.1, 0.2))
        ax[0].set_yticklabels(["%0.1f"%ele for ele in np.arange(0.0, 1.1, 0.2)], fontsize=tl)
        ax[0].set_xticks(np.arange(0.0, 100.1, 20))
        ax[0].set_xticklabels(["%2.0f"%ele for ele in np.arange(0.0, 100.1, 20)], fontsize=tl)
        ax[0].set_xlim([t_dummy[0]-273.16, t_dummy[-1]-273.16])
        ax[0].set_ylim([0.0, 1.0])
        ax[0].plot([T_constrain, T_constrain], [0., 1.], linewidth=2, color='k', linestyle='--')
        ax[0].plot([t_dummy[0]-273.16, t_dummy[-1]-273.16], [pB_constrain, pB_constrain], linewidth=2, color='r', linestyle='--')
        ax[0].legend(fontsize=fl)

        ax[1].set_ylabel("Extinction coefficient (" + "$\epsilon$" + ")", fontsize=fl)
        ax[1].plot(t_dummy-273.16, epsilon_ds, color='k', label='AB')
        ax[1].plot(t_dummy-273.16, epsilon_ds_conf, color='r', label="AB'")
        ax[1].plot(t_dummy-273.16, epsilon_ss, color='b', label="A")
        ax[1].set_xlim([t_dummy[0]-273.16, t_dummy[-1]-273.16])
        ax[1].set_xticks(np.arange(0.0, 100.1, 20))
        ax[1].set_xticklabels(["%2.0f"%ele for ele in np.arange(0.0, 100.1, 20)], fontsize=tl)
        ax[1].set_yticks(np.linspace(min_epsilon, max_epsilon, 10))    
        ax[1].set_yticklabels(["%1.2f"%ele for ele in np.linspace(min_epsilon, max_epsilon, 10)], fontsize=tl)    
        ax[1].set_xlabel("Temperature (" + "$^\circ$" + "C)", fontsize=fl) 
        ax[1].legend(fontsize=fl)
       
        ax[2].plot(temperature-273.16, absorbance, color='k', linewidth=0.0, marker='o', markersize=2)
        ax[2].plot(t_dummy-273.16, abs_net, color='b')
        ax[2].set_xlim([t_dummy[0]-273.16, t_dummy[-1]-273.16])
        ax[2].set_xticks(np.arange(0.0, 100.1, 20))
        ax[2].set_xticklabels(["%2.0f"%ele for ele in np.arange(0.0, 100.1, 20)], fontsize=tl)
        ax[2].set_title('$\chi^{2}$' + ' ' + str('.3f'%redchi2) + '\n')
        fig.suptitle(filename)
        plt.show()
   
        # Flip the sign of delHm and delSm to correspond to annealing
        # Compute delG
        popt[7] = popt[7] * -1.0
        popt[8] = popt[8] * -1.0
        [mss, bss, mds, bds, mds2, bds2, delHconf, delHm, delSm] = popt
        # Get delSconf
        # delG_conf = delH_conf - T * delS_conf
        delSconf = (delHconf - delG_constrain) / (273.16 + T_constrain)

        delGm_25 = delHm - ((273.16+25)*delSm)
        delGm_37 = delHm - ((273.16+37)*delSm)
        delGconf_25 = delHconf - ((273.16+25)*delSconf)
        delGconf_37 = delHconf - ((273.16+37)*delSconf)
        #print "delGm_25C ", delGm_25
        #print "delGm_37C ", delGm_37

        # Now, if delGconf_25 is negative 
        # this happens because the GS and ES are swapped
        if delGconf_25 < 0.0:
            print "*** SWAPPING ***"
            delHm_new = delHm + delHconf
            delSm_new = delSm + delSconf
            delHconf_new = -1.0  * delHconf
            delSconf_new = -1.0 * delSconf

            mds_new = mds
            bds_new = bds
            mds = mds2
            bds = bds2
            mds2 = mds_new
            bds2 = bds_new
 
            delHm = delHm_new
            delSm = delSm_new
            delHconf = delHconf_new 
            delSconf = delSconf_new
            popt = [mss, bss, mds, bds, mds2, bds2, delHconf, delHm, delSm]

        delGm_25 = delHm - ((273.16+25)*delSm)
        delGm_37 = delHm - ((273.16+37)*delSm)
        delGconf_25 = delHconf - ((273.16+25)*delSconf)
        delGconf_37 = delHconf - ((273.16+37)*delSconf)
        print "delGm_25C ", delGm_25
        print "delGm_37C ", delGm_37
   
        print_result(popt, pcov, ['mss', 'bss', 'mds', 'bds', 'mds2', 'bds2', 'delHconf', 'delHm', 'delSm'])
        labels = ['mss', 'bss', 'mds', 'bds', 'mds2', 'bds2', 'delHconf', 'delSconf', 'pB_constrain', 'T_constrain', 'delHm', 'delSm', 'delGconf_25', 'delGconf_37', 'delGm_25C', 'delGm_37C', 'Ct', 'redchi2', 'aic', 'bic']
        result = np.array(list(popt[:7]) + [delSconf, pB_constrain, T_constrain] + list(popt[7:]) + [delGconf_25, delGconf_37, delGm_25, delGm_37, Ct, redchi2, aic, bic])
 
    return [labels, result]

def print_result(parameters, errors, labels):
    ''' Print the results of the fitting '''
    for dummy in range(len(labels)):
        print labels[dummy], parameters[dummy]

def remove_spaces(line):
    line_temp = []
    for ele in line:
        if ele != "":
            line_temp.append(ele)
    return line_temp

def read_data(filename):
    global mode, Ct

    # Read the UV data file
    f = open(filename, "r")
    temperature = []
    absorbance = []
    # Read the 1st two lines
    #line = f.readline()
    #line = f.readline()
    # Now read the rest of the lines
    for line in f.readlines():
        if line[0] == 'D' or len(line) <= 2:
            break
        line = line.strip("\n")
        line = line.split("\t")
        line = remove_spaces(line)
        line[0] = float(line[0])
        line[1] = float(line[1])
        # Convert temperature to Kelvin
        temperature.append(line[0] + 273.16)
        absorbance.append(line[1])
    # Close the file
    f.close()

    # Now call the functions
    temperature = np.array(temperature)
    absorbance = np.array(absorbance)
    #return [temperature, absorbance]
    [labels, result] = fit_data(filename, temperature, absorbance)
    
    return [labels, result]

def caller():
    ''' Call the fitting function and op stats '''
    global output_filename

    result_net = []
    labels_net = []
    filenames_net = []
    for filename in os.listdir('.'):
        if filename[-3:] == "txt":
            print "Reading ", filename
            filenames_net.append(filename)
            [labels, result] = read_data(filename)
            labels_net = labels
            if result_net == []:
	        result_net = np.array([list(result)])
	    else:
		result_net = np.concatenate((result_net, [list(result)]), axis=0)

    # Now, compute stats
    result_avg = np.average(result_net, axis=0)
    result_std = np.std(result_net, axis=0)

    #print result_avg
    #print result_std
    #print result_net

    # Now, o/p to file
    f = open(output_filename + ".csv", "w")

    # Write header
    labels.append('mode')
    header = "filename," + labels[0] + ''.join(["," + labels[dummy] for dummy in range(1, len(labels))]) + '\n'
    print header
    f.write(header)

    # Now, write the data
    for dummy in range(np.shape(result_net)[0]):
        data_line = str(filenames_net[dummy]) + "," + str(result_net[dummy][0]) + ''.join(["," + str(result_net[dummy][dummy_int]) for dummy_int in range(1, np.shape(result_net)[1])]) + ',' + str(mode) + '\n'
	f.write(data_line)
  
    # Now, write the stats
    avg_line = "avg," + str(result_avg[0]) + ''.join(["," + str(result_avg[dummy_int]) for dummy_int in range(1, np.shape(result_net)[1])]) + ',' + str(mode) + '\n'
    f.write(avg_line)
    
    std_line = "std," + str(result_std[0]) + ''.join(["," + str(result_std[dummy_int]) for dummy_int in range(1, np.shape(result_net)[1])]) + ',0\n'
    f.write(std_line)

    f.close()


# Gas constant in kcal/mol
R = 8.314/4184
mode = str(sys.argv[1])
Ct = float(sys.argv[2]) * math.pow(10, -6)
output_filename = str(sys.argv[3])
K = 9.0
T_constrain = 25.0
pB_constrain = 7.0/100.0
delG_constrain = -R * (273.16 + T_constrain) * math.log(pB_constrain/(1-pB_constrain))

caller()
