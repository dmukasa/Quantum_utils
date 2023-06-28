import joblib
import datetime
import time
import psutil
import glob
from tqdm import tqdm
from joblib import Parallel, delayed
import numpy as np
from sklearn.kernel_ridge import KernelRidge
from sklearn.model_selection import train_test_split
from scipy import stats
from matplotlib import pyplot as plt

########################### functions ########################### 
def H_kJ_mol(E):
    """Converts Hartree to kJ/mol"""
    return E*27.211*96.485 #kJ/mol
########################## import labels ########################
print("Loading labels")
Energy = np.load("E_arr_gen_06-12-2023.npz", mmap_mode='r')["arr_0"] # 0.001 Gb
print("Done!")

########################## import features ########################
print("Loading features")
feature_arr = np.load('features_arr_gen_06-11-2023.npz', mmap_mode='r')["arr_0"] # 200 Gb

############### Remove and unprocessed rows (i.e. errors) ##########################
non_zero_mask = np.any(feature_arr != 0, axis = (1,2,3))
feature_arr = feature_arr[non_zero_mask]
Energy = Energy[non_zero_mask]

# TESTING: SHORTEN DATASET
#feature_arr = feature_arr[:20000]
#Energy = Energy[:20000]
print("Done!")

##TESTING: SAVE SUBSET OF DATA
#np.savez_compressed('feature_arr_5K', feature_arr)
#np.savez_compressed('Energy_5K', Energy)
########################## split the data ########################
print("split data")
feature_arr = np.reshape(feature_arr, (len(feature_arr), -1))
train_x, test_x, train_y, test_y = train_test_split(feature_arr, Energy, test_size=0.1, random_state=42) #200 + 200 = 400 Gb

# Clear memory
feature_arr = None # 400 - 200 = 200 Gb
print("Done!")

# Print the shapes of train and test sets
print("Train set shapes - train_x: {}, train_y: {}".format(train_x.shape, train_y.shape))
print("Test set shapes - test_x: {}, test_y: {}".format(test_x.shape, test_y.shape))

########################## train KRR ########################
print("Training KRR")
# Define the KRR model
krr_Energy = KernelRidge(alpha=0.001)
krr_Energy.fit(train_x, train_y) # 200 + 180 Gb = 380 Gb
print("Done!")

# Save the trained model to a file
current_datetime = datetime.datetime.now()
formatted_datetime = current_datetime.strftime("%d-%m-%Y_%H:%MPST")
filename = 'krr_' + str(len(Energy)) + '_datapoints_' + formatted_datetime + '_time.joblib'
joblib.dump(krr_Energy, filename)

# Evaluate the model
test_labels = np.reshape(test_y, (len(test_y), 1))

# Use the predict method on the testing data
test_predictions = np.reshape(krr_Energy.predict(test_x), (len(test_x), 1))

# Calculate the absolute errors
average_absolute_test_err = np.average(np.abs(test_predictions - test_labels))
print(average_absolute_test_err)
print("MAE = " + str(average_absolute_test_err), "kJ/mol")

########################## plot the results ########################

# Fit the testing data
test_slope, test_intercept, test_r_value, test_p_value, test_std_err = stats.linregress(np.array(test_predictions).flatten(),np.array(test_labels).flatten())
test_fit_line = lambda x2 : test_slope*x2 + test_intercept

# Plot the test data fit
plt.figure(figsize=(5,5))
arial_font = {'fontname':'Arial'}
#Fit data to predictions
plt.scatter(test_predictions, test_labels, marker = 'o', label = str(test_r_value*test_r_value))

#Define axis and title
plt.xlabel("ML Predicted Energy")
plt.ticklabel_format(style = 'scientific',scilimits=(-5,3))
plt.ylabel("DFT Computed Energy")
plt.xticks()
plt.yticks()
plt.legend(loc='upper left')
plt.show()
# Save the plot to a file
plt.savefig(filename[:-7] + '.png')

import csv

# Float value to save
value = average_absolute_test_err

# CSV file path
csv_file = filename[:-7] + '.csv'

# Save the float value to the CSV file
with open(csv_file, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow([value])
