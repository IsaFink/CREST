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
###             - updated for CREST MSreact              ###
###                On multiple directories               ###
###              - by Isabella Fink-Jensen               ###
############################################################

d=`pwd` #Run from the top of the directory tree that contains the fragment directories of interest.
drug_name="drug_name" #Edit to drug you're working with
#########################################
echo 'starting parallel MSreact runs'

#Create a list of directory that contains the MSinput files
# Initialize an empty array to store directory paths
declare -a dir_paths=()

# Search for all directories containing "MSinput_*" files
while IFS= read -r dir; do
    dir_paths+=("$dir")
done < <(find . -type f -name "MSinput_*" -exec dirname {} \; | sort --unique)

# Print the list of directory paths
Total=0
for path in "${dir_paths[@]}"; do
    echo "$path"
    cd $path #now in the directory containing the MSinput file
    cp /work3/ishof/SK/path.inp .
    let TOTAL=TOTAL+1
    cd d
    last=$path
done > Optimization_directories.txt

echo "Paths to directories containing MSinput files have been saved to 'Optimization_directories.txt'."
echo "There are $TOTAL directories with MSinput files."

# Now we're making a .sh script in each fragment directory that will run the MSreact command
# This script will be submitted to the cluster in the next step
#########################################
# outer loop (runs until all is done)
#########################################

 for path in "${dir_paths[@]}"; do
    cd $path
    for file in MSinput_*; do #There is only one MSinput file in each directory
        if [ -f "$file" ]; then  # Check if it's a file
            input_name=$(basename "$file")
            echo "Directory: $path, Input Name: $input_name" 
            #if made with MSinput, the input_name will be MSinput_{name_of_experiment}
            #where {name_of_experiment} is the parameter that is varied and its value. 
            #Example: MSinput_etemp2500
            #name of experiment is input_name without MSinput_
            name_of_experiment=${input_name#MSinput_}
        fi
    done
    VZ_DIR=`pwd`
# Put everything following this line into the "command file" until EOF

    cat > MSmultiReact_${drug_name}_${name_of_experiment}.sh <<EOF 
#!/bin/sh 
### General options 
### -- specify queue -- 
#BSUB -q hpc
### -- set the job Name -- 
#BSUB -J MSmultiReact_${drug_name}_${name_of_experiment}
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
PATH=$PATH:/zhome/92/7/155378/miniforge3/envs/crest_env/bin

/zhome/92/7/155378/miniforge3/envs/crest_env/bin/crest /work3/ishof/SK/CREST/${drug_name}/${drug_name}H.xyz -msinput ${input_name} -msreact --msmolbar --mslargeprint --T 16 --chrg 1

EOF
#  this is the job

        bsub < MSmultiReact_${drug_name}_${name_of_experiment}.sh  #run the command file on cluster 
#
        echo "JOB  "MSmultiReact_${drug_name}_${name_of_experiment}.sh"  STARTED"
   cd d
 done  # stop outer loop
exit