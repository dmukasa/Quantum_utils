start=`date +%s`
# Meant to go to /groups/GaoGroup/dmukasa/ML_DFT/Dimers_out_copy
# and run a PM3 on each file. Each file should be in a directory
# with the .xyz file alone
path_to_dir=/groups/GaoGroup/dmukasa/ML_DFT/Dimers_out_copy

#Go to the directory
cd ${path_to_dir} 
i=0

for file in *; do
    #save the path
    WD=${path_to_dir}/${file}
    #Enter the directory
    cd $WD

    #save the file name
    filename=${file}
        
    #Make the inp
    python /groups/GaoGroup/dmukasa/ML_DFT/Dimer_editor/Make_inp.py ${filename} ${WD}
    #Run Orca
    /home/dmukasa/orca/orca_4_2_1_linux_x86-64_openmpi216/orca ${WD}/${filename}.inp > ${WD}/${filename}.out &

    if (( $i % 4 == 0 ))
    then
        wait
    fi

    #FOR TESTING
    #ls
    #pwd

    #i=$((i+1))
    #if [[ $i == '300' ]]
    #then
    #    break
    #fi
done

end=`date +%s`
runtime=$((end-start))
echo $runtime
