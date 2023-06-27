import time
import psutil
import glob
from tqdm import tqdm
from joblib import Parallel, delayed
import numpy as np

####################################### Failure mode detection functions #######################################
def check_for_tmp(directory):
    """Function to find directory with files with the .tmp extension or a
    .tmp.gz extension. This indicates that the orca .out file did not finish 
    it's calculation (i.e. indicating an error)"""

    # Create a pattern to match files with the extensions .tmp.gz or .tmp
    pattern1 = f"{directory}/*.{'tmp.gz' if directory.endswith('/') else 'tmp.gz'}"
    pattern2 = f"{directory}/*.{'tmp' if directory.endswith('/') else 'tmp'}"


    # Use glob to find matching files
    file1 = glob.glob(pattern1)
    file2 = glob.glob(pattern2)

    # Check if any files were found
    if file1 or file2:
        return True
    else:
        return False    

def check_for_gbw(directory):
    """Function to find directory with files with the .gbw extension or a
    .gbw.gz extension. This indicates that the orca file finished its 
    calculation. The lack of this file indicates some failure to run."""

    # Create a pattern to match files with the extensions .tmp.gz or .tmp
    pattern1 = f"{directory}/*.{'gbw.gz' if directory.endswith('/') else 'gbw.gz'}"
    pattern2 = f"{directory}/*.{'gbw' if directory.endswith('/') else 'gbw'}"

    # Use glob to find matching files
    file1 = glob.glob(pattern1)
    file2 = glob.glob(pattern2)

    # Check if any files were found
    if file1 or file2:
        return False
    else:
        return True

##############################################################################################################

#def ERR_lst_loop_tmp(i):
#    identifier = str(i).zfill(6)
#    directory = "/groups/GaoGroup/dmukasa/ML_DFT/QM9_PM3_copy/dsgdb9nsd_" + identifier + ".xyz"
#    if check_for_tmp(directory) or check_for_gbw(directory):
#        index = int(identifier)-1
#        print("ERROR IN: ", directory)
#        print("Index = ", index)
#        ERR_lst_tmp.append(index)
#
#startTime = time.time()
#
#ERR_lst_tmp = []
#num_files = 133885
#
#print(check_directory_for_files("/groups/GaoGroup/dmukasa/ML_DFT/QM9_PM3_copy/dsgdb9nsd_033702.xyz"))
#Parallel(n_jobs=-1,backend="threading")(delayed(ERR_lst_loop_tmp)(i) for i in tqdm(range(1,num_files+1)))
#
#executionTime = (time.time() - startTime)
#print('Time to make names list in minutes: ' + str(executionTime/60))
#
#print(ERR_lst_tmp)
#np.savez_compressed('ERR_lst_tmp.npz', ERR_lst_tmp)
