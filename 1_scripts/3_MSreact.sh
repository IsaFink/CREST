#!/bin/sh 
### General options 
### -- specify queue -- 
#BSUB -q hpc
### -- set the job Name -- 
#BSUB -J MSReact_1
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
#path to workflow_env
PATH=$PATH:/zhome/92/7/155378/miniforge3/envs/crest_env/bin

#remember to change <drug_name> to the name of the drug you are working with
/zhome/92/7/155378/miniforge3/envs/crest_env/bin/crest /work3/ishof/SK/CREST/<drug_name>H.xyz -msreact --mslargeprint --T 16 --chrg 1
