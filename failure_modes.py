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

