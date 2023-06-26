import time
import psutil
import glob
from tqdm import tqdm
from joblib import Parallel, delayed
import numpy as np
from Matrix_utils import *

def gunzip_file(input_file, output_file):
    with gzip.open(input_file, "rb") as f_in:
        with open(output_file, "wb") as f_out:
            f_out.write(f_in.read())

ERR_lst = np.load('ERR_lst_total.npz')["arr_0"] 
ERR_lst_tmp = np.load('ERR_lst_tmp.npz')["arr_0"]
ERR_lst = np.append(ERR_lst, ERR_lst_tmp)
ERR_lst += 1

print(ERR_lst)

def main_loop(i):
    if i not in ERR_lst:
        identifier = str(i).zfill(6)
        path_to_out = "/groups/GaoGroup/dmukasa/ML_DFT/QM9_PM3_copy/dsgdb9nsd_" + identifier + ".xyz/"
        path_to_out += "dsgdb9nsd_" + identifier + ".out"

        # If the .out file isnt unzipped unzip it
        if not os.path.isfile(path_to_out):
            print("Gunziping " + path_to_out)
            gunzip_file(path_to_out + ".gz", path_to_out)
        print("Generating matrix"+str(identifier))
        #generate the approriate matricies
        one_e = Get_one_electron(path_to_out, max_atoms)
        overlap = Get_overlap(path_to_out, max_atoms)
        S12 = Get_S12(path_to_out, max_atoms)
        Density = Get_Density(path_to_out, max_atoms)
        Fock = Get_Fock(path_to_out, max_atoms)

        #add them to the features_arr
        features_arr[i-1] += np.array(([one_e], [overlap], [S12], [Density], [Fock])).reshape(max_atoms,max_atoms,5)



#count the number of files in QM9
path = "/groups/GaoGroup/dmukasa/ML_DFT/QM9"
num_files = 0
for filename in os.listdir(path):
    if filename != ".DS_Store" and filename != "._.DS_Store":
        num_files += 1

# Remove files with errors
num_files -= len(ERR_lst)
#choose the maximum size of the system
max_atoms = 200

#make array to hold all labels
#features_arr = np.zeros((num_files, max_atoms, max_atoms, 5))

#Parallel(n_jobs=-1,backend="threading")(delayed(main_loop)(i) for i in tqdm(range(1,num_files+1)))

# Save and compress the feature array
savez_compressed('/groups/GaoGroup/dmukasa/ML_DFT/QM9_editor/features_arr.npz', features_arr)
