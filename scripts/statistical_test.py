import pandas as pd
import numpy as np
import sys

number_models = len(sys.argv) - 1
files = [pd.read_csv(sys.argv[dummy+1]) for dummy in range(number_models)]
number_samples = len(files[0]) - 2
sample_names = list(files[0]['filename'][:-2])

# Loop over samples
for dummy in range(number_samples):
    sample = sample_names[dummy]
    print "Sample ", sample
    aic_values = []
    bic_values = []

    # Loop over models
    for dummy_model in range(number_models):
        model_data = files[dummy_model].loc[files[dummy_model]['filename'] == sample]
        aic_values.append(model_data['aic'].iloc[0])
        bic_values.append(model_data['bic'].iloc[0])

    # Compute model weights
    aic_values = np.array(aic_values)
    bic_values = np.array(bic_values)
    aic_values = aic_values - np.amin(aic_values)
    bic_values = bic_values - np.amin(bic_values)
    #print aic_values
    #print bic_values

    weights_aic = np.exp(-0.5 * aic_values)      
    weights_bic = np.exp(-0.5 * bic_values)      
    #print weights_aic
    #print weights_bic
    weights_aic = weights_aic / np.sum(weights_aic)
    weights_bic = weights_bic / np.sum(weights_bic)
    print "AIC weights = ", weights_aic
    print "BIC weights = ", weights_bic
    print "AIC: Favored model is ", np.argmax(weights_aic) + 1, " with probability ", np.amax(weights_aic)
    print "BIC: Favored model is ", np.argmax(weights_bic) + 1, " with probability ", np.amax(weights_bic)
 

    print
            
