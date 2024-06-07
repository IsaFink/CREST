import os
import numpy as np

# Setting defaults
fc_rep_default = 0.5
etemp_default = 5000
frag_dist_default = 0
distthr_attr_default = 4.0
fc_attr_default = -0.5

# Parameter values to edit as wished or outcomment, if not needed
fc_rep_values = np.array([0.5, 1, 2]) #<real>, default 0.5
etemp_values = np.array([2500, 5000, 7500]) #<real> default 5000 K
frag_dist_values = np.array([0, 2, 5]) #<real>, default 0
distthr_attr_values = np.array([3.0, 4.0, 5.0]) #<real> default 4.0 Å
fc_attr_values = np.array([-2, -1, -0.5]) #<real> default -0.5

params = ['frag_dist', 'fc_rep', 'etemp', 'attractions']

# Iterate through each parameter and generate directories and files
for param in params:
    if param == 'frag_dist':
        for i in frag_dist_values:
            directory = f'{param}/MSreact_{param}_{i}'
            os.makedirs(directory, exist_ok=True)
            filepath = os.path.join(directory, f'MSinput_{param}_{i}')
            with open(filepath, 'w') as f:
                text = f'frag_dist {i} \nfc_rep {fc_rep_default} \netemp {etemp_default} \nfc_attr {fc_attr_default} \ndistthr_attr {distthr_attr_default}'
                f.write(text)
    
    elif param == 'fc_rep':
        for i in fc_rep_values:
            directory = f'{param}/MSreact_{param}_{i}'
            os.makedirs(directory, exist_ok=True)
            filepath = os.path.join(directory, f'MSinput_{param}_{i}')
            with open(filepath, 'w') as f:
                text = f'frag_dist {frag_dist_default} \nfc_rep {i} \netemp {etemp_default} \nfc_attr {fc_attr_default} \ndistthr_attr {distthr_attr_default}'
                f.write(text)
    
    elif param == 'etemp':
        for i in etemp_values:
            directory = f'{param}/MSreact_{param}_{i}'
            os.makedirs(directory, exist_ok=True)
            filepath = os.path.join(directory, f'MSinput_{param}_{i}')
            with open(filepath, 'w') as f:
                text = f'frag_dist {frag_dist_default} \nfc_rep {fc_rep_default} \netemp {i} \nfc_attr {fc_attr_default} \ndistthr_attr {distthr_attr_default}'
                f.write(text)
    
    elif param == 'attractions':
        for k in fc_attr_values:
            for m in distthr_attr_values:
                directory = f'{param}/MSreact_fc_{k}_dist_{m}'
                os.makedirs(directory, exist_ok=True)
                filepath = os.path.join(directory, f'MSinput_fc_{k}_dist_{m}')
                with open(filepath, 'w') as f:
                    text = f'frag_dist {frag_dist_default} \nfc_rep {fc_rep_default} \netemp {etemp_default} \nfc_attr {k} \ndistthr_attr {m}'
                    f.write(text)

                
'''
figures of merit: matching peaks/number of peaks --> make a matrix of this
Correct peaks or not?
Decrease amount of H and H2 lost - more effective fragmentation runs

fc_rep  --> Breaks up fragments easier --> more fragmentation and more peaks --> higher energy transitions in theory
    0.5 (default), 1, 2

etemp --> higher temperature --> more fragmentation --> more peaks 
    --> lower temperature --> less fragmentation --> more realistic peaks? 
2500, 5000, 7500 

These two can be a pairwise full scans --> more complicated :) #keep it realistic- reduce H and H2 loss
fc_attr --> attraction between hydrogens and LMO centers --> will saturate heavy atoms? --> possibly removing half of the fragments if we get the right charges based on the hydrogen transfer
distthr_attr <real> --> increasing this will increase the stickiness of the hydrogens
# no feeling - have fun :)

frag_dist --> pushes fragments away from the centroid after everything else has been done. Is only useful for the xtb path finder.
--> try 2 (practically doubling the RMSD) and 5 (might be doing more harm than good) and see what happens. 0 probably will allow for more conformational changes.



#optimizing the reaction path finder:

$path
   nrun=1 - determines number of refinement cycles - increase to 5
        forward  barrier (kcal)  :   122.098
        backward barrier (kcal)  :   129.747
        reaction energy  (kcal)  :    -7.649
        opt. pull strength       :     0.050
        norm(g) at est. TS, point: 0.12503  14  
        # all of the above times 5 -- check in order to see if the optimization is good or not?
   
   npoint=25 - determines how many points along the path #25 is good for pretty pictures.
   anopt=10 - how many times each point is optimized along the path -> increasing this lowers energy to a minimum, but that costs time.
   
   #stuff we don't edit unless we get a lot of failed runs.
   kpush=0.003
   kpull=-0.015
   ppull=0.05
   alp=1.2
$end

 
'''