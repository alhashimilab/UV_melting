uv_multiple.py
- performs 2-state fitting of UV data

uv_multiple_threestate.py 
- performs 3-state fitting of UV data

uv_multiple_threestate_constrain.py
- performs 3-state fitting of UV data
- while constraining the population of the excited state 
- the temperature and population of the measurement is specified by the T_constrain and pB_constrain variables at the end of the script

Instructions for usage
- place .py scripts in folder along with .txt files containing melting data. THese .txt files should contain temperature, absorbance data separated by spaces, with each temperature, absorbance pair on a new line
- execute by typing "python script_name.py sample_type concentration output_csv"
  - where sample_type is either 'duplex' or 'hairpin'
  - concentration is the concentration of the sample in uM
  - output_csv is the name of the output csv file containing the results
- it is advisable to keep the initial guesses for the curve fitting for the 3-state fits equal to that obtained from the 2-state fit. These are specified by the p0 variables prior to the calls to the curve_fit function

To perform statistical tests on fit_outputs
- type pyton statistical_test.py output1.csv output2.csv
- where output1.csv and output2.csv are the fit outputs from two different fits on the same UV melting data
- the program will print out the AIC/BIC values for each of the individual melting curves  

plot_all*.py - they are the scripts for plotting the fits to the UV data for each type of sample
- type python script_name.py to generate a pdf file containing the plots
