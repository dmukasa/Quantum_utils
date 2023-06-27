import time
import psutil
import glob
from tqdm import tqdm
from joblib import Parallel, delayed
import numpy as np
from ERR_lst_loop import *
from Matrix_utils import *

############################### Nessecary functions ################################################################################
def gunzip_file(input_file, output_file):
    with gzip.open(input_file, "rb") as f_in:
        with open(output_file, "wb") as f_out:
            f_out.write(f_in.read())

def H_kJ_mol(E):
    """Converts Hartree to kJ/mol"""
    return E*27.211*96.485 #kJ/mol
################################ Generate the Energy array and ferature array at the same time #############################################
# Do this while monitoring for features activley

def Gen_loop(i):
    # Loop through all indicies in our loop rage 
    # Define the identifier id used iun the directory name
    identifier = str(i).zfill(6)
    # Define the index of the array we are at for this i
    index = i-1
    # Save the directory name
    directory = "/groups/GaoGroup/dmukasa/ML_DFT/QM9_PM3_copy/dsgdb9nsd_" + identifier + ".xyz"

    # If failure modes are detected to be true skip this value
    if check_for_tmp(directory) or check_for_gbw(directory):
        print("ERROR IN: ", directory, "... PASSING")
        pass

    # If no failure modes are detected then contiue with adding this value to the array
    else:
        # Define the path to the out file
        path_to_out = directory + "/dsgdb9nsd_" + identifier + ".out"

        # If the .out file isn't unzipped unzip it
        if not os.path.isfile(path_to_out):
            print("Gunziping " + path_to_out)
            gunzip_file(path_to_out + ".gz", path_to_out)

        # Generate the approriate matricies
        one_e = Get_one_electron(path_to_out, max_atoms)
        overlap = Get_overlap(path_to_out, max_atoms)
        S12 = Get_S12(path_to_out, max_atoms)
        Density = Get_Density(path_to_out, max_atoms)
        Fock = Get_Fock(path_to_out, max_atoms)

        # Add them to the features_arr
        features_arr[index] += np.array(([one_e], [overlap], [S12], [Density], [Fock])).reshape(5, max_atoms,max_atoms)
  
        # Add the Energies to E_arr
        E_arr[index] += Energy[index]

# Count the number of files in QM9
path = "/groups/GaoGroup/dmukasa/ML_DFT/QM9_PM3_copy"
num_files = 0
for filename in os.listdir(path):
    if filename[0] != ".":
        num_files += 1

# Load in the labels
# collumns = ['tag', 'identifier','Rot_const1', 'Rot_const2','Rot_const3', 
#          'Dipole moment', 'polarizability', 'E_Homo',
#          'E_Lumo', 'E_Gap', 'Electronic spatial extent', 
#          'ZPE','E(0K)', 'E(298.15K)', 'H', 'G', 'C_v']
QM9_lables = np.genfromtxt("QM9_labels.csv", delimiter=',')

# Delete the first collum, this just contains a tag gdb
QM9_lables = QM9_lables[:,1:]

#Order by identifier
QM9_lables = QM9_lables[QM9_lables[:, 0].argsort()]

# Define energy array
Energy = H_kJ_mol(QM9_lables[:num_files,11])

#make array to hold all features
max_atoms = 200
features_arr = np.zeros((num_files, 5, max_atoms, max_atoms))
E_arr = np.zeros((num_files,1))

# Generate the arrays
Parallel(n_jobs=-1,backend="threading")(delayed(Gen_loop)(i) for i in tqdm(range(1,num_files+1)))

# Save and compress the feature array
savez_compressed('/groups/GaoGroup/dmukasa/ML_DFT/QM9_editor/features_arr_gen_06-11-2023.npz', features_arr)
savez_compressed('/groups/GaoGroup/dmukasa/ML_DFT/QM9_editor/E_arr_gen_06-11-2023.npz', E_arr)
