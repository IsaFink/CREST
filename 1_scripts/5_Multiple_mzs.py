
import os
home='/work3/ishof/SK/CREST'
import sys
sys.path.append(os.path.join(home,'1_scripts')) #append folder with scripts
from mZ import *
import pickle
import pandas as pd


# Add the /modules directory to the Python path

#Adjust if need be
#Mass should be round up to the nearest 10
name_of_drug, mass_of_drug = 'Verapamil', 450 
#name_of_drug, mass_of_drug = 'Caffeine', 200
'''Would you like a gif for each of the valid transitions?'''
GIF_wanted=True



work_dir=os.path.join(home,name_of_drug)
experimental_datafile= os.path.join(home,f'{name_of_drug}/{name_of_drug}_exp.csv')

'''For a drug that doesn't have the experimental datafile, set it to None'''
#experimental_datafile=None

directories={}
parameter=['fc_rep', 'etemp', 'attractions', 'frag_dist']
for param in parameter:
    directories[param]={}
    for root, dirs, files in os.walk(work_dir): 
        levels=root.split('/')
        if levels[-1].startswith('MSreact') and levels[-2]==param:
            exp_no=levels[-1]
            directories[param][exp_no]={}
            directories[param][exp_no]['root']=root

'''Can be used for a drug with only one experiment'''
if name_of_drug=='Verapamil':
    param='default'
    exp_no=0
    directories[param]={}
    directories[param][exp_no]={}
    directories[param][exp_no]['root']=work_dir

#in each root there is a file that starts with MSinput
for param in directories.keys():
    for exp_no in directories[param].keys():
        root=directories[param][exp_no]['root']
        for file in os.listdir(root):
            if file.startswith('MSinput'):
                directories[param][exp_no]['input']=file
                with open(f'{root}/{file}', 'r') as f:
                    lines=f.readlines()
                    directories[param][exp_no]['frag_dist']=lines[0].strip().split(' ')[-1]
                    directories[param][exp_no]['fc_rep']=lines[1].strip().split(' ')[-1]
                    directories[param][exp_no]['etemp']=lines[2].strip().split(' ')[-1]
                    directories[param][exp_no]['fc_attr']=lines[3].strip().split(' ')[-1]
                    directories[param][exp_no]['distthr_atrr']=lines[4].strip().split(' ')[-1]

print("An overview of the different experiments:")               
print_nested_dict(directories)


'''Saving plots in a directory called Results'''                        
try:
    os.mkdir(os.path.join(home,'2_Supervision/Results'))
except:
    pass
try:
    os.mkdir(os.path.join(home,f'2_Supervision/Results/{name_of_drug}'))
except:
    pass

for param in directories.keys():
    for exp_no in directories[param].keys():
        root=directories[param][exp_no]['root']
        if param=='default':
            name_of_experiment='Initial_run'
            print(f'Beginning of experiment {name_of_experiment} \n \n \n')
            pair_dict, mz_dict, mz_peaks, result, error_counts=dicts_w_info(root, name_of_experiment, GIF_wanted)
            print(f'{name_of_experiment} \t {len(mz_dict)} \t {len(mz_peaks)}')
            print(f'End of experiment {name_of_experiment} \n ')
            directories[param][exp_no]['name']=name_of_experiment
            directories[param][exp_no]['mz_dict']=mz_dict
            directories[param][exp_no]['mz_peaks']=mz_peaks
            directories[param][exp_no]['result']=result
            directories[param][exp_no]['error_counts']=error_counts
            if mz_dict:
                plot_mz(mz_dict, experimental_datafile, name_of_experiment, name_of_drug, mass_of_drug)
                copy_png_files(root, name_of_experiment, name_of_drug)
        elif param!='attractions':
            value=directories[param][exp_no][param]
            name_of_experiment=f'{param}_{value}'
            print(f'\n \n \n Beginning of experiment {name_of_experiment} \n')
            pair_dict, mz_dict, mz_peaks, result, error_counts=dicts_w_info(root, name_of_experiment, GIF_wanted)
            print(f'{name_of_experiment} \t {len(mz_dict)} \t {len(mz_peaks)}')
            print(f'End of experiment {name_of_experiment} \n \n \n')
            directories[param][exp_no]['name']=name_of_experiment
            directories[param][exp_no]['mz_peaks']=mz_peaks
            directories[param][exp_no]['mz_dict']=mz_dict
            directories[param][exp_no]['result']=result
            directories[param][exp_no]['error_counts']=error_counts
            if mz_dict:
                plot_mz(mz_dict, experimental_datafile ,name_of_experiment, name_of_drug, mass_of_drug)
                copy_png_files(root, name_of_experiment, name_of_drug)
        else:
            value1=directories[param][exp_no]['fc_attr']
            value2=directories[param][exp_no]['distthr_atrr']
            name_of_experiment=f'{param}_fc_attr_{value1}_distthr_atrr_{value2}'
            print(f'Beginning of experiment {name_of_experiment} \n \n \n')
            pair_dict, mz_dict, mz_peaks, result, error_counts=dicts_w_info(root, name_of_experiment, GIF_wanted)
            print(f'{name_of_experiment} \t {len(mz_dict)} \t {len(mz_peaks)}')
            print(f'End of experiment {name_of_experiment} \n ')
            directories[param][exp_no]['name']=name_of_experiment
            directories[param][exp_no]['mz_dict']=mz_dict
            directories[param][exp_no]['mz_peaks']=mz_peaks
            directories[param][exp_no]['result']=result
            directories[param][exp_no]['error_counts']=error_counts
            if mz_dict:
                plot_mz(mz_dict, experimental_datafile ,name_of_experiment, name_of_drug, mass_of_drug)
                copy_png_files(root, name_of_experiment, name_of_drug)
            #directories[param][exp_no]['pair_dict']=pair_dict
            
print('A summary of results can be found in the following overview:')
print_nested_dict(directories)

print('The results have been saved in the Results directory')
#Save the dictionary for later retrieval
with open(os.path.join(home,f'2_Supervision/Results/{name_of_drug}/directories.pkl'), 'wb') as f:
    pickle.dump(directories, f)
    
#read the dictionary file
with open(os.path.join(home,f'2_Supervision/Results/{name_of_drug}/directories.pkl'), 'rb') as f:
    directories=pickle.load(f)


#read experimental data
exp_data=pd.read_csv(experimental_datafile)
exp_peaks=exp_data['M/Z_fragment'].tolist()

pm=4 #plus minus interval to look for peaks in the range of
print(f'The matching peaks are identified within a range of {pm} g/mol from the experimental peaks.')
#Evaluating the matching peaks
for param in directories.keys():
    for exp_no in directories[param].keys():
        directories[param][exp_no]['matching_peaks']={}
        for peak in directories[param][exp_no]['mz_dict'].keys():
            for exp_peak in exp_peaks:
                if exp_peak-pm < peak < exp_peak+pm:
                    directories[param][exp_no]['matching_peaks'][peak]=directories[param][exp_no]['mz_dict'][peak]['simulations']
        directories[param][exp_no]['number_of_matching_peaks']=len(directories[param][exp_no]['matching_peaks'])
        
                    
'''Printing a summary of the matching peaks of each experiment'''
for param in directories.keys():
    for exp_no in directories[param].keys():
        if 'result' in directories[param][exp_no].keys():
            print(f'----####----\n{directories[param][exp_no]['name']}')
            print('----####----')
            print_nested_dict(directories[param][exp_no]['result'])
            print_nested_dict(directories[param][exp_no]['error_counts'])
            print(directories[param][exp_no]['mz_peaks'])
            print(f'There is/are {directories[param][exp_no]['number_of_matching_peaks']} matching peaks')
            print('Matching peaks: and their simulations')
            print_nested_dict(directories[param][exp_no]['matching_peaks'])

            

print('The directories dictionary has been overwritten, adding the matching peaks')
with open(os.path.join(home,f'2_Supervision/Results/{name_of_drug}/directories.pkl'), 'wb') as f:
    pickle.dump(directories, f)
    
#save the print_mested_dict output to a file
with open(os.path.join(home,f'2_Supervision/Results/{name_of_drug}/directories.txt'), 'w') as f:
    f.write(nested_dict_layout(directories))

    
#Now we move all relevant gifs and Energypaths to the Results directory
#first we find the relevant simulations for each mazz_peak and then we move the gifs and energy paths
for param in directories.keys():
    for exp_no in directories[param].keys():
        if 'matching_peaks' in directories[param][exp_no].keys():
            for peak in directories[param][exp_no]['matching_peaks'].keys():
                for sim in directories[param][exp_no]['matching_peaks'][peak]:
                    if directories[param][exp_no]['mz_dict'][peak].get('Excluded'):
                        if sim in directories[param][exp_no]['mz_dict'][peak]['Excluded']:
                            print(f'The simulation {sim} has been excluded from the gif transfer')
                    else:
                        root=directories[param][exp_no]['root']
                        name_of_experiment=directories[param][exp_no]['name']
                        gif_path=os.path.join(root, f'4_transitions_paths_gifs/transpath_{name_of_experiment}_{sim}.gif')
                        shutil.copy(gif_path, os.path.join(home,f'2_Supervision/Results/{name_of_drug}'))
                        print(f'The gif {sim} has been copied to the Results directory')

for param in directories.keys():
    for exp_no in directories[param].keys():
        if 'matching_peaks' in directories[param][exp_no].keys():     
            print(f'Experiment {directories[param][exp_no]['name']}') 
            print_nested_dict(directories[param][exp_no]['matching_peaks'])


