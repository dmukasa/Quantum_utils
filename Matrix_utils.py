import numpy as np
import csv
import os
from numpy import savez_compressed
import gzip
import shutil
import re

def H_kJ_mol(E):
    """Converts Hartree to kJ/mol"""
    return E*27.211*96.485 #kJ/mol

def MAPE(Y_actual,Y_Predicted):
    mape = np.mean(np.abs((Y_actual - Y_Predicted)/Y_actual))*100
    return mape

def MAE(Y_actual,Y_Predicted):
    err = Y_actual - Y_Predicted
    mae = np.average(np.abs(err))
    return mae

def pad_array(arr, max_atoms=200):
    #Pads array with zeros to hte right and bottom 
    # output is max_atoms x max_atoms array
    diff = max_atoms-len(arr)
    return np.pad(arr, ((0,diff),(0,diff)), 'constant')
#import functions for feature extraction
def Get_one_electron(Path_to_out_file, max_atoms = 100):
    ## Meant to grab to one electron
    ## matrix from the .out file assuming proper 
    ## print settings. Output is an nxn array
    
    #open the file
    file = open(Path_to_out_file, 'r')
    Lines = file.readlines()
    
    #Check the file format. Should be .out
    
    if re.search(" O   R   C   A ", Lines[2]):
        #print('this is an orca output')
        we_continue = True
    else :
        print('dunno this output format!')
        we_continue = False
    
    #If the format's good continue 
    if we_continue:
        # Find the line which say "ONE ELECTRON MATRIX (AU)"
        # Mark 3 lines after this as where the matrix starts
        start = 0
        for line in Lines:
            if line == "ONE ELECTRON MATRIX (AU)\n":
                break
            ####### MAJOR CHANGE FOR SYSTEMS WITH LESS THAN 6 
            ####### ORGBITALS, COUNT THIS AND SAVE IT FOR LATER
            elif line.split()[:3] == ['Basis', 'Dimension', 'Dim']:
                n = int(line.split()[-1])
                start += 1
            else:
                start += 1
        start += 3
        
        # Now find the end whcih will be the first line after
        # start which is not an integrer
        end = start
        while True:
            try:
                for i in range(start, len(Lines)-1):
                    if len(Lines[i].split()) > 0:
                        test = int(Lines[i].split()[0])
                        end += 1
            except ValueError:
                end -= 1
                break
                
        #Identify where the breks are that indicate the collum 
        breaks = []
        breaks_index = []
        for i in range(start,end):
            num_spaces =  len(Lines[i]) - len(Lines[i].lstrip())
            if num_spaces > 10:
                breaks.append(Lines[i].strip())
                breaks_index.append(i)
        # Define the final matrix size, it will be an (n x n) matrix
        # where n is the last number in the last break line
#         n = int(breaks[-1].split()[-1]) +1 #MAJOR CHANGE
        one_electron_matrix = np.zeros((n,n))
    
        if n > 6: #MAJOR CHANGE
            # From range(start,breaks_index[0]) is the first n rows
            # and up to 6 collums. Then from 
            # range(breaks_index[0],breaks_index[1]) is the next n rows
            # and up to 6 collums. Add these to one_electron_matrix
            prev_col_len = []
            for i in range(0,len(breaks_index)+1):
                # Condition for start block
                if i == 0:
                    starting_line = start
                    ending_line = breaks_index[i]
                    count = 0

                    for j in range(starting_line,ending_line):
                        test = np.array(Lines[j].split()[1:]).astype(float)
                        one_electron_matrix[count][0:len(test)] += test
                        count += 1
                    prev_col_len.append(len(test))
                # Condition for middle blocks
                elif  0 < i < len(breaks_index):#-1:
                    starting_line = breaks_index[i-1] + 1
                    ending_line = breaks_index[i]
                    count = 0

                    for j in range(starting_line,ending_line):
                        test = np.array(Lines[j].split()[1:]).astype(float)
                        one_electron_matrix[count][prev_col_len[-1]:len(test)+prev_col_len[-1]] += test
                        count += 1
                    prev_col_len.append(len(test)+prev_col_len[-1])
                # Condition for end block
                elif  i == len(breaks_index):
                    starting_line = breaks_index[i-1] + 1
                    ending_line = end + 1
                    count = 0

                    for j in range(starting_line,ending_line):
                        test = np.array(Lines[j].split()[1:]).astype(float)
                        one_electron_matrix[count][prev_col_len[-1]:len(test)+prev_col_len[-1]] += test
                        count += 1
        else: #MAJOR CHANGE
            prev_col_len = []
            for i in range(0,len(breaks_index)+1):
                # Condition for start block
                if i == 0:
                    starting_line = start
                    ending_line = start+n-1#breaks_index[i]
                    count = 0
                    
                    for j in range(starting_line,ending_line):
                        test = np.array(Lines[j].split()[1:]).astype(float)
                        one_electron_matrix[count][0:len(test)] += test
                        count += 1
                    prev_col_len.append(len(test))
#         print(np.array_str(one_electron_matrix, precision=6, suppress_small=True))
    #pad the final array to be 170x170
    return pad_array(one_electron_matrix, max_atoms)

def Get_overlap(Path_to_out_file, max_atoms = 100):
    ## Meant to grab to Overlap
    ## matrix from the .out file assuming proper 
    ## print settings. Output is an nxn array
    
    #open the file
    file = open(Path_to_out_file, 'r')
    Lines = file.readlines()
    
    #Check the file format. Should be .out
    
    if re.search(" O   R   C   A ", Lines[2]):
        #print('this is an orca output')
        we_continue = True
    else :
        print('dunno this output format!')
        we_continue = False
    
    #If the format's good continue 
    if we_continue:
        # Find the line which say "OVERLAP MATRIX (AU)"
        # Mark 3 lines after this as where the matrix starts
        start = 0
        for line in Lines:
            if line == "OVERLAP MATRIX\n":
                break
            ####### MAJOR CHANGE FOR SYSTEMS WITH LESS THAN 6 
            ####### ORGBITALS, COUNT THIS AND SAVE IT FOR LATER
            elif line.split()[:3] == ['Basis', 'Dimension', 'Dim']:
                n = int(line.split()[-1])
                start += 1
            else:
                start += 1
        start += 3
        
        # Now find the end whcih will be the first line after
        # start which is not an integrer
        end = start
        while True:
            try:
                for i in range(start, len(Lines)-1):
                    if len(Lines[i].split()) > 0:
                        test = int(Lines[i].split()[0])
                        end += 1
            except ValueError:
                end -= 1
                break
                
        #Identify where the breks are that indicate the collum 
        breaks = []
        breaks_index = []
        for i in range(start,end):
            num_spaces =  len(Lines[i]) - len(Lines[i].lstrip())
            if num_spaces > 10:
                breaks.append(Lines[i].strip())
                breaks_index.append(i)
        # Define the final matrix size, it will be an (n x n) matrix
        # where n is the last number in the last break line
#         n = int(breaks[-1].split()[-1]) +1 #MAJOR CHANGE
        overlap_matrix = np.zeros((n,n))
    
        
        if n > 6: #MAJOR CHANGE
            # From range(start,breaks_index[0]) is the first n rows
            # and up to 6 collums. Then from 
            # range(breaks_index[0],breaks_index[1]) is the next n rows
            # and up to 6 collums. Add these to overlap_matrix
            prev_col_len = []
            for i in range(0,len(breaks_index)+1):
                # Condition for start block
                if i == 0:
                    starting_line = start
                    ending_line = breaks_index[i]
                    count = 0

                    for j in range(starting_line,ending_line):
                        test = np.array(Lines[j].split()[1:]).astype(float)
                        overlap_matrix[count][0:len(test)] += test
                        count += 1
                    prev_col_len.append(len(test))
                # Condition for middle blocks
                elif  0 < i < len(breaks_index):#-1:
                    starting_line = breaks_index[i-1] + 1
                    ending_line = breaks_index[i]
                    count = 0

                    for j in range(starting_line,ending_line):
                        test = np.array(Lines[j].split()[1:]).astype(float)
                        overlap_matrix[count][prev_col_len[-1]:len(test)+prev_col_len[-1]] += test
                        count += 1
                    prev_col_len.append(len(test)+prev_col_len[-1])
                # Condition for end block
                elif  i == len(breaks_index):
                    starting_line = breaks_index[i-1] + 1
                    ending_line = end + 1
                    count = 0

                    for j in range(starting_line,ending_line):
                        test = np.array(Lines[j].split()[1:]).astype(float)
                        overlap_matrix[count][prev_col_len[-1]:len(test)+prev_col_len[-1]] += test
                        count += 1
        else:
            # From range(start,breaks_index[0]) is the first n rows
            # and up to 6 collums. Then from 
            # range(breaks_index[0],breaks_index[1]) is the next n rows
            # and up to 6 collums. Add these to one_electron_matrix
            prev_col_len = []
            for i in range(0,len(breaks_index)+1):
                # Condition for start block
                if i == 0:
                    starting_line = start
                    ending_line = start+n-1#breaks_index[i]
                    count = 0
                    
                    for j in range(starting_line,ending_line):
                        test = np.array(Lines[j].split()[1:]).astype(float)
                        overlap_matrix[count][0:len(test)] += test
                        count += 1
                    prev_col_len.append(len(test))
    return pad_array(overlap_matrix, max_atoms)

def Get_S12(Path_to_out_file, max_atoms = 100):
    ## Meant to grab to Spin^-1/2
    ## matrix from the .out file assuming proper 
    ## print settings. Output is an nxn array
    
    #open the file
    file = open(Path_to_out_file, 'r')
    Lines = file.readlines()
    
    #Check the file format. Should be .out
    
    if re.search(" O   R   C   A ", Lines[2]):
        #print('this is an orca output')
        we_continue = True
    else :
        print('dunno this output format!')
        we_continue = False
    
    #If the format's good continue 
    if we_continue:
        # Find the line which say "S**(-1/2) MATRIX"
        # Mark 3 lines after this as where the matrix starts
        start = 0
        for line in Lines:
            if line == "S**(-1/2) MATRIX\n":
                break
            ####### MAJOR CHANGE FOR SYSTEMS WITH LESS THAN 6 
            ####### ORGBITALS, COUNT THIS AND SAVE IT FOR LATER
            elif line.split()[:3] == ['Basis', 'Dimension', 'Dim']:
                n = int(line.split()[-1])
                start += 1
            else:
                start += 1
        start += 3
        
        # Now find the end whcih will be the first line after
        # start which is not an integrer
        end = start
        while True:
            try:
                for i in range(start, len(Lines)-1):
                    if len(Lines[i].split()) > 0:
                        test = int(Lines[i].split()[0])
                        end += 1
            except ValueError:
                end -= 1
                break
                
        #Identify where the breks are that indicate the collum 
        breaks = []
        breaks_index = []
        for i in range(start,end):
            num_spaces =  len(Lines[i]) - len(Lines[i].lstrip())
            if num_spaces > 10:
                breaks.append(Lines[i].strip())
                breaks_index.append(i)
        # Define the final matrix size, it will be an (n x n) matrix
        # where n is the last number in the last break line
#         n = int(breaks[-1].split()[-1]) +1 #MAJOR CHANGE 
        S12_matrix = np.zeros((n,n))
    
        if n > 6: #MAJOR CHANGE
            # From range(start,breaks_index[0]) is the first n rows
            # and up to 6 collums. Then from 
            # range(breaks_index[0],breaks_index[1]) is the next n rows
            # and up to 6 collums. Add these to S12_matrix
            prev_col_len = []
            for i in range(0,len(breaks_index)+1):
                # Condition for start block
                if i == 0:
                    starting_line = start
                    ending_line = breaks_index[i]
                    count = 0

                    for j in range(starting_line,ending_line):
                        test = np.array(Lines[j].split()[1:]).astype(float)
                        S12_matrix[count][0:len(test)] += test
                        count += 1
                    prev_col_len.append(len(test))
                # Condition for middle blocks
                elif  0 < i < len(breaks_index):#-1:
                    starting_line = breaks_index[i-1] + 1
                    ending_line = breaks_index[i]
                    count = 0

                    for j in range(starting_line,ending_line):
                        test = np.array(Lines[j].split()[1:]).astype(float)
                        S12_matrix[count][prev_col_len[-1]:len(test)+prev_col_len[-1]] += test
                        count += 1
                    prev_col_len.append(len(test)+prev_col_len[-1])
                # Condition for end block
                elif  i == len(breaks_index):
                    starting_line = breaks_index[i-1] + 1
                    ending_line = end + 1
                    count = 0

                    for j in range(starting_line,ending_line):
                        test = np.array(Lines[j].split()[1:]).astype(float)
                        S12_matrix[count][prev_col_len[-1]:len(test)+prev_col_len[-1]] += test
                        count += 1
        else:
            # From range(start,breaks_index[0]) is the first n rows
            # and up to 6 collums. Then from 
            # range(breaks_index[0],breaks_index[1]) is the next n rows
            # and up to 6 collums. Add these to one_electron_matrix
            prev_col_len = []
            for i in range(0,len(breaks_index)+1):
                # Condition for start block
                if i == 0:
                    starting_line = start
                    ending_line = start+n-1#breaks_index[i]
                    count = 0
                    
                    for j in range(starting_line,ending_line):
                        test = np.array(Lines[j].split()[1:]).astype(float)
                        S12_matrix[count][0:len(test)] += test
                        count += 1
                    prev_col_len.append(len(test))
    #pad the final array to be 170x170
    return pad_array(S12_matrix, max_atoms)

def Get_Density(Path_to_out_file, max_atoms = 100):
    ## Meant to grab to Spin^-1/2
    ## matrix from the .out file assuming proper 
    ## print settings. Output is an nxn array
    
    #open the file
    file = open(Path_to_out_file, 'r')
    Lines = file.readlines()
    
    #Check the file format. Should be .out
    
    if re.search(" O   R   C   A ", Lines[2]):
        #print('this is an orca output')
        we_continue = True
    else :
        print('dunno this output format!')
        we_continue = False
    
    #If the format's good continue 
    if we_continue:
        # Find the line which say "DENSITY"
        # Mark 3 lines after this as where the matrix starts
        start = 0
        for line in Lines:
            if line == "DENSITY\n":
                break
            ####### MAJOR CHANGE FOR SYSTEMS WITH LESS THAN 6 
            ####### ORGBITALS, COUNT THIS AND SAVE IT FOR LATER
            elif line.split()[:3] == ['Basis', 'Dimension', 'Dim']:
                n = int(line.split()[-1])
                start += 1
            else:
                start += 1
        start += 3
        
        # Now find the end whcih will be the first line after
        # start which is not an integrer
        end = start
        while True:
            try:
                for i in range(start, len(Lines)-1):
                    if len(Lines[i].split()) > 0:
                        test = int(Lines[i].split()[0])
                        end += 1
            except ValueError:
                end -= 1
                break
                
        #Identify where the breks are that indicate the collum 
        breaks = []
        breaks_index = []
        for i in range(start,end):
            num_spaces =  len(Lines[i]) - len(Lines[i].lstrip())
            if num_spaces > 10:
                breaks.append(Lines[i].strip())
                breaks_index.append(i)
        # Define the final matrix size, it will be an (n x n) matrix
        # where n is the last number in the last break line
#         n = int(breaks[-1].split()[-1]) +1 #MAJOR CHANGE
        DENSITY_matrix = np.zeros((n,n))
    
        if n > 6: #MAJOR CHANGE
            # From range(start,breaks_index[0]) is the first n rows
            # and up to 6 collums. Then from 
            # range(breaks_index[0],breaks_index[1]) is the next n rows
            # and up to 6 collums. Add these to DENSITY_matrix
            prev_col_len = []
            for i in range(0,len(breaks_index)+1):
                # Condition for start block
                if i == 0:
                    starting_line = start
                    ending_line = breaks_index[i]
                    count = 0

                    for j in range(starting_line,ending_line):
                        test = np.array(Lines[j].split()[1:]).astype(float)
                        DENSITY_matrix[count][0:len(test)] += test
                        count += 1
                    prev_col_len.append(len(test))
                # Condition for middle blocks
                elif  0 < i < len(breaks_index):#-1:
                    starting_line = breaks_index[i-1] + 1
                    ending_line = breaks_index[i]
                    count = 0

                    for j in range(starting_line,ending_line):
                        test = np.array(Lines[j].split()[1:]).astype(float)
                        DENSITY_matrix[count][prev_col_len[-1]:len(test)+prev_col_len[-1]] += test
                        count += 1
                    prev_col_len.append(len(test)+prev_col_len[-1])
                # Condition for end block
                elif  i == len(breaks_index):
                    starting_line = breaks_index[i-1] + 1
                    ending_line = end + 1
                    count = 0

                    for j in range(starting_line,ending_line):
                        test = np.array(Lines[j].split()[1:]).astype(float)
                        DENSITY_matrix[count][prev_col_len[-1]:len(test)+prev_col_len[-1]] += test
                        count += 1
        else: #MAJOR CHANGE
            # From range(start,breaks_index[0]) is the first n rows
            # and up to 6 collums. Then from 
            # range(breaks_index[0],breaks_index[1]) is the next n rows
            # and up to 6 collums. Add these to one_electron_matrix
            prev_col_len = []
            for i in range(0,len(breaks_index)+1):
                # Condition for start block
                if i == 0:
                    starting_line = start
                    ending_line = start+n#breaks_index[i]
                    count = 0
                    
                    for j in range(starting_line,ending_line):
                        test = np.array(Lines[j].split()[1:]).astype(float)
                        DENSITY_matrix[count][0:len(test)] += test
                        count += 1
                    prev_col_len.append(len(test))
    #pad the final array to be 170x170
    return pad_array(DENSITY_matrix, max_atoms)

def Get_Fock(Path_to_out_file, max_atoms = 100):
    ## Meant to grab to Fock
    ## matrix from the .out file assuming proper 
    ## print settings. Output is an nxn array
    
    #open the file
    file = open(Path_to_out_file, 'r')
    Lines = file.readlines()
    
    #Check the file format. Should be .out
    
    if re.search(" O   R   C   A ", Lines[2]):
        #print('this is an orca output')
        we_continue = True
    else :
        print('dunno this output format!')
        we_continue = False
    
    #If the format's good continue 
    if we_continue:
        #Make a reversed Lines
        #input list
        lst=[10, 11, 12, 13, 14, 15]
        # the above input can also be given as
        # lst=list(map(int,input().split()))
        Revered_lines=[] # empty list
        for i in Lines:
            #reversing the list
            Revered_lines.insert(0,i)
            
        # Find the line which say "Fock matrix for operator "
        # Mark 3 lines after this as where the matrix starts
        start = len(Lines)
        
        for line in Revered_lines:
            tmp = line[:-2]
            if tmp == "Fock matrix for operator ":
                break
            ####### MAJOR CHANGE FOR SYSTEMS WITH LESS THAN 6 
            ####### ORGBITALS, COUNT THIS AND SAVE IT FOR LATER
            elif line.split()[:3] == ['Basis', 'Dimension', 'Dim']:
                n = int(line.split()[-1])
                start += 1
            else:
                start -= 1
        start += 1
        # Now find the end whcih will be the first line after
        # start which is not an integrer
        end = start
        while True:
            try:
                for i in range(start, len(Lines)-1):
                    if len(Lines[i].split()) > 0:
                        test = int(Lines[i].split()[0])
                        end += 1
            except ValueError:
                end -= 1
                break
                
        #find n with a forward pass
        for line in Lines:
            ####### MAJOR CHANGE FOR SYSTEMS WITH LESS THAN 6 
            ####### ORGBITALS, COUNT THIS AND SAVE IT FOR LATER
            if line.split()[:3] == ['Basis', 'Dimension', 'Dim']:
                n = int(line.split()[-1])
                
        #Identify where the breks are that indicate the collum 
        breaks = []
        breaks_index = []
        for i in range(start,end):
            num_spaces =  len(Lines[i]) - len(Lines[i].lstrip())
            if num_spaces > 10:
                breaks.append(Lines[i].strip())
                breaks_index.append(i)
        # Define the final matrix size, it will be an (n x n) matrix
        # where n is the last number in the last break line
#         n = int(breaks[-1].split()[-1]) +1 #MAJOR CHANGE
        Fock_matrix = np.zeros((n,n))
    
        if n > 6: #MAJOR CHANGE
            # From range(start,breaks_index[0]) is the first n rows
            # and up to 6 collums. Then from 
            # range(breaks_index[0],breaks_index[1]) is the next n rows
            # and up to 6 collums. Add these to Fock_matrix
            prev_col_len = []
            for i in range(0,len(breaks_index)+1):
                # Condition for start block
                if i == 0:
                    starting_line = start
                    ending_line = breaks_index[i]
                    count = 0

                    for j in range(starting_line,ending_line):
                        test = np.array(Lines[j].split()[1:]).astype(float)
                        Fock_matrix[count][0:len(test)] += test
                        count += 1
                    prev_col_len.append(len(test))
                # Condition for middle blocks
                elif  0 < i < len(breaks_index):#-1:
                    starting_line = breaks_index[i-1] + 1
                    ending_line = breaks_index[i]
                    count = 0

                    for j in range(starting_line,ending_line):
                        test = np.array(Lines[j].split()[1:]).astype(float)
                        Fock_matrix[count][prev_col_len[-1]:len(test)+prev_col_len[-1]] += test
                        count += 1
                    prev_col_len.append(len(test)+prev_col_len[-1])
                # Condition for end block
                elif  i == len(breaks_index):
                    starting_line = breaks_index[i-1] + 1
                    ending_line = end + 1
                    count = 0

                    for j in range(starting_line,ending_line):
                        test = np.array(Lines[j].split()[1:]).astype(float)
                        Fock_matrix[count][prev_col_len[-1]:len(test)+prev_col_len[-1]] += test
                        count += 1
        else: #MAJOR CHANGE
            # From range(start,breaks_index[0]) is the first n rows
            # and up to 6 collums. Then from 
            # range(breaks_index[0],breaks_index[1]) is the next n rows
            # and up to 6 collums. Add these to one_electron_matrix
            prev_col_len = []
            for i in range(0,len(breaks_index)+1):
                # Condition for start block
                if i == 0:
                    starting_line = start
                    ending_line = start+n#breaks_index[i]
                    count = 0
                    
                    for j in range(starting_line,ending_line):
                        test = np.array(Lines[j].split()[1:]).astype(float)
                        Fock_matrix[count][0:len(test)] += test
                        count += 1
                    prev_col_len.append(len(test))
    #pad the final array to be 170x170
    return pad_array(Fock_matrix, max_atoms)


def Get_CM(Path_to_xyz_file, max_atoms = 56):
    cm = CoulombMatrix(n_atoms_max = max_atoms,)
    
#     path = "/Volumes/MIP_DESIGN/QM9/QM9_PM3/dsgdb9nsd_"
#     identifier = str(i).zfill(6)
#     xyz_path = path + identifier + ".xyz/dsgdb9nsd_" + identifier + ".xyz"
#     # We need to make a new xyz FILENAME_dscribe.xyz
#     # this file will have an added line for the second line
#     dscribe_path = path + identifier + ".xyz/dsgdb9nsd_" + identifier + "_dscribe" + ".xyz"
    dscribe_path = Path_to_xyz_file
    
    # Check the length of the file, the first two lines are commented out
    reading = open(dscribe_path, "r")
    contents = reading.readlines()
    num_atoms = int(contents[0])
    len_file = len(contents)
    reading.close()
    if num_atoms + 1 == len_file:
        #add a new line to the second line
        write = open(dscribe_path, "w")
        contents.insert(1, '\n')
        write.writelines(contents)
        write.close()
            

    
    if exists(dscribe_path):
        while True:
            try:
                structure = read(dscribe_path)
#                 QM9_10K_CMs[i-1] +=  cm.create(structure) 
                return np.reshape(cm.create(structure), (max_atoms, max_atoms))
                break
            except ValueError:
                break
    else: 
        #make a copy of the xyz file to the dscribe path
        file_dscribe = open(dscribe_path, "w")
        file = open(xyz_path)
        Lines = file.readlines()
        # Insert an end line in the comment line
        Lines.insert(1, '\n')
        # Write this to the dscribe file
        file_dscribe.writelines(Lines)
        file_dscribe.close()
        file.close()    
    # Now check again
    if exists(dscribe_path):
        while True:
            try:
                structure = read(dscribe_path)
#                 QM9_10K_CMs[i-1] +=  cm.create(structure) 
                return np.reshape(cm.create(structure), (max_atoms, max_atoms))
                break
            except ValueError:
                break
    else:
        print("FATAL ERROR IN ", identifier)
        
