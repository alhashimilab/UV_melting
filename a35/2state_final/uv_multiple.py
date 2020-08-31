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
    # Gas constant in kcal/mol
    R = 8.314 / 4184
    # Get the exponential factor
    exp_factor = np.exp(((1/Tm) - (1/(t+273.16))) * delH / R)
    # Get the fraction of double stranded species for a duplex
    f = ((1 + (4*exp_factor)) - np.sqrt(1 + (8*exp_factor)))/(4*exp_factor)
    # Use f to get absorbance of duplex/ss mixture
    absorbance = (((mds*t)+bds)*f) + (((mss*t)+bss)*(1-f))
    return absorbance

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

def fit_data(filename, temperature, absorbance, mode, Ct):
    # Call the appropriate fitting function
    global K
    R = 8.314 / 4184
    Ct = float(Ct)
    t_dummy = np.linspace(0, 100, 500) + 273.16

    if mode == 'hairpin': 
        print "Hairpin"
        p0 = np.array([0.0003, 0.85, 0.0006, 0.95, -64, 340])
        [popt, pcov] = curve_fit(uv_hairpin, temperature, absorbance, p0)

        print_result(popt, pcov, ['mds', 'bds', 'mss', 'bss', 'delH', 'Tm'])
        delS = (popt[4] / popt[5]) 
        delG_25 = popt[4] - ((273.16+25)*delS)
        delG_37 = popt[4] - ((273.16+37)*delS)
        print "delS ", delS
        print "delG_25C ", delG_25
        print "delG_37C ", delG_37

        # populations
        [mds_opt, bds_opt, mss_opt, bss_opt, delH_opt, Tm_opt] = popt
        exp_factor = np.exp(((1/Tm_opt) - (1/(t_dummy))) * delH_opt / R)
        # Get the fraction of hairpin that is intact
        f_hairpin = (exp_factor / (1 + exp_factor))
        f_ss = 1. - f_hairpin

        # extinction coefficients
        epsilon_ss = bss_opt + (mss_opt*(t_dummy-273.16))
        epsilon_ds = bds_opt + (mds_opt*(t_dummy-273.16))
        min_epsilon = min(list(epsilon_ss) + list(epsilon_ds))
        max_epsilon = max(list(epsilon_ss) + list(epsilon_ds))
        min_epsilon = min_epsilon - ((min_epsilon%0.050))
        max_epsilon = max_epsilon + (0.050 - (min_epsilon%0.050))

        # statistics
        # define error as given by std. deviation of 1st 20 points
        N = float(len(absorbance))
        ideal_absorbance = uv_hairpin(temperature, *popt)
        error_absorbance = np.std(absorbance[:20])
        chi2 = np.sum(np.square((absorbance - ideal_absorbance)/error_absorbance)) 
        redchi2 = chi2 / (N - K)
        residual = np.sum(np.square(absorbance - ideal_absorbance))
 
        # Compute AIC and BIC values
        # First AIC
        if (N / K) >= 40:
            aic = (N * math.log(residual/N)) + (2 * K)
        else:
            aic = (N * math.log(residual/N)) + (2 * K) + ((2 * K * (K+1))/(N-K-1))
        bic = (N * math.log(residual/N)) + (K*math.log(N))
 
        # Plotting    
        fig, ax = plt.subplots(1, 3)
        fig.suptitle(str(filename), fontsize=fl)
        ax[0].set_ylabel("Population", fontsize=fl)
        ax[0].set_xlabel("Temperature (" + "$^\circ$" + "C)", fontsize=fl)
        ax[0].plot(t_dummy - 273.16, f_hairpin, color='k', label='Afold')
        ax[0].plot(t_dummy - 273.16, f_ss, color='r', label='Amelt')
        ax[0].set_yticks(np.arange(0.0, 1.1, 0.2))
        ax[0].set_yticklabels(["%0.1f"%ele for ele in np.arange(0.0, 1.1, 0.2)], fontsize=tl)
        ax[0].set_xticks(np.arange(0.0, 100.1, 20))
        ax[0].set_xticklabels(["%2.0f"%ele for ele in np.arange(0.0, 100.1, 20)], fontsize=tl)
        ax[0].set_xlim([t_dummy[0]-273.16, t_dummy[-1]-273.16])
        ax[0].set_ylim([0.0, 1.0])
        ax[0].legend(fontsize=fl)

        ax[1].set_ylabel("Extinction coefficient (" + "$\epsilon$" + ")", fontsize=fl)
        ax[1].plot(t_dummy-273.16, epsilon_ds, color='k', label='Afold')
        ax[1].plot(t_dummy-273.16, epsilon_ss, color='r', label="Amelt")
        ax[1].set_xticks(np.arange(0.0, 100.1, 20))
        ax[1].set_xticklabels(["%2.0f"%ele for ele in np.arange(0.0, 100.1, 20)], fontsize=tl)
        ax[1].set_xlim([t_dummy[0]-273.16, t_dummy[-1]-273.16])
        #ax[1].set_yticks(np.linspace(min_epsilon, max_epsilon, 10))
        #ax[1].set_yticklabels(["%1.2f"%ele for ele in np.linspace(min_epsilon, max_epsilon, 10)], fontsize=tl)
        ax[1].set_xlabel("Temperature (" + "$^\circ$" + "C)", fontsize=fl)
        ax[1].legend(fontsize=fl)

        ax[2].set_title('$\chi^{2}$' + ' ' + str('%.3f'%redchi2))
        ax[2].set_ylabel('Absorbance', fontsize=fl)
        ax[2].set_xlabel("Temperature (" + "$^\circ$" + "C)", fontsize=fl)
        ax[2].plot(temperature, absorbance, color='k', linewidth=3, label='data')
        ax[2].plot(np.linspace(5,100,500), uv_hairpin(np.linspace(5,100,500), *popt), color='b', label='fit')
        ax[2].set_xticks(np.arange(0.0, 100.1, 20))
        ax[2].set_xticklabels(["%2.0f"%ele for ele in np.arange(0.0, 100.1, 20)], fontsize=tl)
        ax[2].set_xlim([t_dummy[0]-273.16, t_dummy[-1]-273.16])
        ax[2].legend(fontsize=fl)
        plt.show()

        labels = ['mds', 'bds', 'mss', 'bss', 'delH', 'Tm(K)', 'Tm(C)','delS', 'delG_25C', 'delG_37C', 'Ct', 'redchi2', 'aic', 'bic']
        result = np.array(list(popt) + [popt[5]-273.16,delS, delG_25, delG_37, Ct, redchi2, aic, bic])

    elif mode == 'duplex':
        print "Duplex"
        p0 = np.array([0.001, 0.5, 0.001, 0.7, -70, 310])
        [popt, pcov] = curve_fit(uv_duplex, temperature, absorbance, p0)
        print_result(popt, pcov, ['mds', 'bds', 'mss', 'bss', 'delH', 'Tm'])
        delS = (popt[4] / popt[5]) - (R * math.log(Ct * 0.000001 / 2))
        delG_25 = popt[4] - ((273.16+25)*delS)
        delG_37 = popt[4] - ((273.16+37)*delS)
        print "delS ", delS
        print "delG_25C ", delG_25
        print "delG_37C ", delG_37

        [mds_opt, bds_opt, mss_opt, bss_opt, delH_opt, Tm_opt] = popt
        exp_factor = np.exp(((1/Tm_opt) - (1/(t_dummy))) * delH_opt / R)
        # Get the fraction of double stranded species for a duplex
        f_duplex = ((1 + (4*exp_factor)) - np.sqrt(1 + (8*exp_factor)))/(4*exp_factor)
        f_ss = 1 - f_duplex

        # extinction coefficients
        epsilon_ss = bss_opt + (mss_opt*(t_dummy-273.16))
        epsilon_ds = bds_opt + (mds_opt*(t_dummy-273.16))
        min_epsilon = min(list(epsilon_ss) + list(epsilon_ds))
        max_epsilon = max(list(epsilon_ss) + list(epsilon_ds))
        min_epsilon = min_epsilon - ((min_epsilon%0.050))
        max_epsilon = max_epsilon + (0.050 - (min_epsilon%0.050))

        # statistics
        N = float(len(absorbance))
        error_absorbance = np.std(absorbance[:20])
        ideal_absorbance = uv_duplex(temperature, *popt)
        chi2 = np.sum(np.square((absorbance - ideal_absorbance)/error_absorbance)) 
        redchi2 = chi2 / (N - K)
        residual = np.sum(np.square(absorbance - ideal_absorbance))

        # Compute AIC and BIC values
        # First AIC
        if (N / K) >= 40:
            aic = (N * math.log(residual/N)) + (2 * K)
        else:
            aic = (N * math.log(residual/N)) + (2 * K) + ((2 * K * (K+1))/(N-K-1))
        bic = (N * math.log(residual/N)) + (K*math.log(N))
        
        # Plotting    
        fig, ax = plt.subplots(1, 3)
        fig.suptitle(str(filename), fontsize=fl)
        ax[0].set_ylabel("Population", fontsize=fl)
        ax[0].set_xlabel("Temperature (" + "$^\circ$" + "C)", fontsize=fl)
        ax[0].plot(t_dummy - 273.16, f_duplex, color='k', label='Afold')
        ax[0].plot(t_dummy - 273.16, f_ss, color='r', label='Amelt')
        ax[0].set_yticks(np.arange(0.0, 1.1, 0.2))
        ax[0].set_yticklabels(["%0.1f"%ele for ele in np.arange(0.0, 1.1, 0.2)], fontsize=tl)
        ax[0].set_xticks(np.arange(0.0, 100.1, 20))
        ax[0].set_xticklabels(["%2.0f"%ele for ele in np.arange(0.0, 100.1, 20)], fontsize=tl)
        ax[0].set_xlim([t_dummy[0]-273.16, t_dummy[-1]-273.16])
        ax[0].set_ylim([0.0, 1.0])
        ax[0].legend(fontsize=fl)

        ax[1].set_ylabel("Extinction coefficient (" + "$\epsilon$" + ")", fontsize=fl)
        ax[1].plot(t_dummy-273.16, epsilon_ds, color='k', label='Afold')
        ax[1].plot(t_dummy-273.16, epsilon_ss, color='r', label="Amelt")
        ax[1].set_xticks(np.arange(0.0, 100.1, 20))
        ax[1].set_xticklabels(["%2.0f"%ele for ele in np.arange(0.0, 100.1, 20)], fontsize=tl)
        ax[1].set_xlim([t_dummy[0]-273.16, t_dummy[-1]-273.16])
        #ax[1].set_yticks(np.linspace(min_epsilon, max_epsilon, 10))
        #ax[1].set_yticklabels(["%1.2f"%ele for ele in np.linspace(min_epsilon, max_epsilon, 10)], fontsize=tl)
        ax[1].set_xlabel("Temperature (" + "$^\circ$" + "C)", fontsize=fl)
        ax[1].legend(fontsize=fl)

        # absorbance
        ax[2].set_title('$\chi^{2}$' + ' ' + str('%.3f'%redchi2))
        ax[2].set_ylabel('Absorbance', fontsize=fl)
        ax[2].set_xlabel("Temperature (" + "$^\circ$" + "C)", fontsize=fl)
        ax[2].plot(temperature, absorbance, color='k', linewidth=3, label='data')
        ax[2].plot(np.linspace(5,100,500), uv_duplex(np.linspace(5,100,500), *popt), color='b', label='fit')
        ax[2].set_xticks(np.arange(0.0, 100.1, 20))
        ax[2].set_xticklabels(["%2.0f"%ele for ele in np.arange(0.0, 100.1, 20)], fontsize=tl)
        ax[2].set_xlim([t_dummy[0]-273.16, t_dummy[-1]-273.16])
        ax[2].legend(fontsize=fl)

        plt.show()
   
        labels = ['mds', 'bds', 'mss', 'bss', 'delH', 'Tm(K)', 'Tm(C)','delS', 'delG_25C', 'delG_37C', 'Ct', 'redchi2', 'aic', 'bic']
        result = np.array(list(popt) + [popt[5]-273.16,delS, delG_25, delG_37, Ct, redchi2, aic, bic])
 
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

def read_data(filename, mode, Ct):
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
        temperature.append(line[0])
        absorbance.append(line[1])
    # Now call the functions
    temperature = np.array(temperature)
    absorbance = np.array(absorbance)
    [labels, result] = fit_data(filename, temperature, absorbance, mode, Ct)
 
    # Close the file
    f.close()
    
    return [labels, result]

def caller():
    ''' Call the fitting function and op stats '''
    result_net = []
    labels_net = []
    filenames_net = []
    for filename in os.listdir('.'):
        if filename[-3:] == "txt":
            print "Reading ", filename
            [labels, result] = read_data(filename, sys.argv[1], sys.argv[2])
            labels_net = labels
            if result_net == []:
	        result_net = np.array([list(result)])
	    else:
		result_net = np.concatenate((result_net, [list(result)]), axis=0)
            filenames_net.append(filename)

    # Now, compute stats
    result_avg = np.average(result_net, axis=0)
    result_std = np.std(result_net, axis=0)

    # Now, add some more columns to the output
    labels.append('mode')
     
    #print result_avg
    #print result_std
    #print result_net

    # Now, o/p to file
    f = open(str(sys.argv[3]) + ".csv", "w")

    # Write header
    header = 'filename,' + labels[0] + ''.join(["," + labels[dummy] for dummy in range(1, len(labels))]) + '\n'
    print header
    f.write(header)

    # Now, write the data
    for dummy in range(np.shape(result_net)[0]):
        data_line = str(filenames_net[dummy]) + "," + str(result_net[dummy][0]) + ''.join(["," + str(result_net[dummy][dummy_int]) for dummy_int in range(1, np.shape(result_net)[1])]) + ',' + str(sys.argv[1]) + '\n'
	f.write(data_line)
  
    # Now, write the stats
    avg_line = "avg," + str(result_avg[0]) + ''.join(["," + str(result_avg[dummy_int]) for dummy_int in range(1, np.shape(result_net)[1])]) + ',' + str(sys.argv[1]) + '\n'
    f.write(avg_line)
    
    std_line = "std," + str(result_std[0]) + ''.join(["," + str(result_std[dummy_int]) for dummy_int in range(1, np.shape(result_net)[1])]) + ',0\n'
    f.write(std_line)

    f.close()

K = 6.0
caller()
