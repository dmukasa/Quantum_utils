import sys
import subprocess as sp
import re

def Run_Opt_inp(FILENAME, WD):
    """ Function to make FILENAME.inp from FILENAME.xyz:

        FILENAME: The afforementioned FILENAME, no .xyz

        WD: The working directory path containing FILENAME.xyz
        For math path/to/current/directory
        no slash at the end
    """
    INP_commands = """! PM3 """
    #INP_commands = """! B3LYP D3BJ 6-31++G**"""


    INP_commands += """
%output
Print[P_OneElec]1
Print[P_Overlap]1
Print[P_KinEn]1
Print[P_S12]1
Print[P_Density]1
Print[P_Iter_F]1
Print[P_OrbEn]1
Print[P_MOs]1
Print[P_SpinDensity]1
Print[P_Fockian]1
Print[P_Mayer]1
Print[P_NatPop]1
#Print[P_Hirshfeld]1
Print[P_Mulliken]1
end

%scf
MaxIter 1000000
end


* xyzfile 0 1"""
    INP_commands += """ """ + FILENAME + """.xyz\n"""

    inp_path = WD + """/""" + FILENAME + """.inp"""

    input = open(inp_path, "w")
    #write the coordinates in it
    input.write(INP_commands)
    #Save the file
    input.close()


Run_Opt_inp(sys.argv[1], sys.argv[2])
