#!/bin/sh 
### General options 
### -- specify queue -- 
#BSUB -q hpc
### -- set the job Name -- 
#BSUB -J MS_opt_MSinput_frag_dist_2
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
PATH=/zhome/92/7/155378/miniforge3/envs/crest_env/bin:/zhome/92/7/155378/miniforge3/condabin:/zhome/92/7/155378/.vscode-server/cli/servers/Stable-e170252f762678dec6ca2cc69aba1570769a5d39/server/bin/remote-cli:/zhome/92/7/155378/Desktop/crest/crest:/lsf/local/bin:/lsf/10.1/linux3.10-glibc2.17-x86_64/etc:/lsf/10.1/linux3.10-glibc2.17-x86_64/bin:/appl/latex/TexLive19/bin/x86_64-linux:/appl/Modules/4.8.0/bin:/apps/dcc/bin:/usr/bin:/bin:/usr/local/bin:/opt/puppetlabs/bin:.:/appl/steno/sw/apps/vmd/1.9.3/bin:/zhome/92/7/155378/Desktop/crest/crest:/zhome/92/7/155378/miniforge3/envs/crest_env/bin:/zhome/92/7/155378/miniforge3/condabin:/zhome/92/7/155378/.vscode-server/cli/servers/Stable-e170252f762678dec6ca2cc69aba1570769a5d39/server/bin/remote-cli:/zhome/92/7/155378/Desktop/crest/crest:/lsf/local/bin:/lsf/10.1/linux3.10-glibc2.17-x86_64/etc:/lsf/10.1/linux3.10-glibc2.17-x86_64/bin:/appl/latex/TexLive19/bin/x86_64-linux:/appl/Modules/4.8.0/bin:/apps/dcc/bin:/usr/bin:/bin:/usr/local/bin:/opt/puppetlabs/bin:.:/appl/steno/sw/apps/vmd/1.9.3/bin:/zhome/92/7/155378/.vscode-server/cli/servers/Stable-e170252f762678dec6ca2cc69aba1570769a5d39/server/bin/remote-cli:/zhome/92/7/155378/Desktop/crest/crest:/lsf/local/bin/:/lsf/10.1/linux3.10-glibc2.17-x86_64/etc:/lsf/10.1/linux3.10-glibc2.17-x86_64/bin:/appl/latex/TexLive19/bin/x86_64-linux:/appl/Modules/4.8.0/bin:/apps/dcc/bin:/usr/bin:/bin:/usr/local/bin:/opt/puppetlabs/bin:.:/appl/steno/sw/apps/vmd/1.9.3/bin:/appl/steno/sw/apps/vmd/1.9.3/bin:/zhome/92/7/155378/miniforge3/envs/crest_env/bin:/zhome/92/7/155378/miniforge3/envs/crest_env/bin

/zhome/92/7/155378/miniforge3/envs/crest_env/bin/crest /work3/ishof/SK/CREST/Caffeine/CaffeineH.xyz -msinput MSinput_frag_dist_2 -msreact --msmolbar --mslargeprint --T 16 --chrg 1

