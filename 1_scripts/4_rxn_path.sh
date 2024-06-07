#Run script in commandline. This script will run xtb on all the directories containing pair.xyz files in the p-folders.
cat /zhome/92/7/155378/.bashrc
PATH=$PATH:/zhome/92/7/155378/miniforge3/envs/crest_env/bin


############################################################
###                                                      ###
###             q - batch Version 2.1                    ###
###                 Stefan Grimme                        ###
###                 Jeroen Koopman                       ###
###                                                      ###
### q-batch script for parallel run on computer clusters ###
###         - updated for xtb reaction path              ###
###              - by Isabella Fink-Jensen               ###
############################################################

d=`pwd` #Run from the top of the directory tree that contains the fragment directories of interest.
drug_name="drug_name" #Edit to drug you're working with
input_file_path="/work3/ishof/SK/path.inp" #Edit to the name of the input file you want to use for xtb

#########################################
echo 'starting parallel xtb rxn path sampling runs'


# Initialize an empty array to store directory paths
declare -a dir_paths=()

# Search for all directories containing "pair.xyz" files
# These are the directories where we want to run xtb
while IFS= read -r dir; do
    dir_paths+=("$dir")
done < <(find . -type f -name "pair.xyz" -exec dirname {} \; | sort --unique)

# Print the list of directory paths
Total=0
for path in "${dir_paths[@]}"; do
    cd $path #now in the directory containing the pair.xyz file
    mkdir xtb_rxn #creating a directory to run xtb in
    cd xtb_rxn
    cp $input_file_path . #copying the path.inp file to the xtb_rxn directory
    let TOTAL=TOTAL+1
    cd d
    last=$path
done > directories_with_pair_xyz.txt

echo "Paths to directories containing 'pair.xyz' files have been saved to 'directories_with_pair_xyz.txt'."
echo "There are $TOTAL directories with 'pair.xyz' files."

# Now we're making a .sh script in each directory that will run the xtb path search
# This script will be submitted to the cluster in the next step
#########################################
# outer loop (runs until all is done)
#########################################
 for path in "${dir_paths[@]}"; do
    cd $path/xtb_rxn
    xtb_DIR=`pwd` #path to the xtb directory
    #get the name of the directory
    basename="${path##*/}"  # This will give us 'p95'
    number="${basename#p}"  # This will give us '95'
    #if there is a file starting with MSinput_ in the directory, get the name of the file
    experiment=$(basename $(find . -maxdepth 1 -type f -name 'MSinput_*'))
    exp_name="${experiment#MSinput_}"  # This will give us 'etemp2500', unless it's empty
    # Put everything following this line into the "command file" until EOF
    cat > rxnpath_${drug_name}_${exp_name}_${number}.sh <<EOF 
#!/bin/sh 
### General options 
### -- specify queue -- 
#BSUB -q hpc
### -- set the job Name -- 
#BSUB -J rxnpath_${drug_name}_${exp_name}_${number}
### -- ask for number of cores (default: 1) -- 
#BSUB -n 16 
### -- specify that the cores must be on the same host -- 
#BSUB -R "span[hosts=1]"
### -- specify that we need 4GB of memory per core/slot -- 
#BSUB -R "rusage[mem=4GB]"
### -- specify that we want the job to get killed if it exceeds 5 GB per core/slot -- 
#BSUB -M 5GB
### -- set walltime limit: hh:mm -- 
#BSUB -W 02:00 
### -- set the email address -- 
# please uncomment the following line and put in your e-mail address,
# if you want to receive e-mail notifications on a non-default address
#BSUB -u ishof@dtu.dk
### -- send notification at start -- 
#BSUB -B 
### -- send notification at completion -- 
#BSUB -N 
### -- Specify the output and error file. %J is the job-id -- 
### -- -o and -e mean append, -oo and -eo mean overwrite -- 
#BSUB -o O_Output_%J.out 
#BSUB -e O_Output_%J.err 

cat /zhome/92/7/155378/.bashrc
PATH=$PATH

cd $xtb_DIR 

/zhome/92/7/155378/miniforge3/envs/crest_env/bin/xtb /work3/ishof/SK/CREST/${drug_name}/${drug_name}H.xyz --path ../pair.xyz --input path.inp
touch ready


EOF

#  this is the job
# VERY VERY IMPORTANT: termination of this job should produce
# a file <ready> in the subdir (e.g. by 'touch ready')
        bsub < rxnpath_${drug_name}_${exp_name}_${number}.sh  #run the command file on cluster 
#
        echo "JOB  "rxnpath_${drug_name}_${exp_name}_${number}.sh"  STARTED"
   cd d
 done  # stop outer loop
exit
