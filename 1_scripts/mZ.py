import os
import rdkit
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.rdmolfiles import MolFromXYZFile
from collections import defaultdict
from rdkit.Chem import rdmolops
from rdkit.Chem import AllChem
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem import Draw
from rdkit.Chem.Draw import MolsToGridImage
from rdkit.Chem import rdDetermineBonds
import imageio.v2 as imageio
import matplotlib.pyplot as plt
import shutil

def print_nested_dict(d, indent=0):
        #sorted by key from smallest to largest
        for key, value in sorted(d.items()):
            if isinstance(value, dict):
                print(' ' * indent + str(key) + ':')
                print_nested_dict(value, indent + 4)
            else:
                print(' ' * indent + str(key) + ': ' + str(value))


def nested_dict_layout(d, indent=0):
    result = ""
    separator = '-' * (indent + 4) + '\n' if indent > 0 else ''
    
    for key, value in sorted(d.items()):
        result += separator
        result += ' ' * indent + str(key) + ':\n'
        
        if isinstance(value, dict):
            result += nested_dict_layout(value, indent + 4)
        elif isinstance(value, pd.DataFrame):
            # Convert DataFrame to a string with an additional indent
            df_string = value.to_string().replace('\n', '\n' + ' ' * (indent + 4))
            result += ' ' * (indent + 4) + df_string + '\n'
        else:
            result += ' ' * (indent + 4) + str(value) + '\n'
    
    result += separator
    return result

#A function to check the amount of simulations included and the amount of simulations excluded:
def check_excluded_included(d, included_before, excluded_before, name_of_last_step):
    included=0
    excluded=0
    for key in d.keys():
        if d[key]['Exclude']==True:
            excluded+=1
        else:
            included+=1
    print(f'\n-----STATUS-----Amount of simulations included and excluded after {name_of_last_step}')
    print(f'The amount of simulations included before was: {included_before}')
    print(f'The amount of simulations excluded before was: {excluded_before}')
    print(f'The total amount of simulations included still is: {included}')
    print(f'The total amount of simulations excluded is: {excluded} \n')
    newly_excluded=excluded-excluded_before
    print(f'-----STATUS-----During {name_of_last_step} {newly_excluded} simulations were excluded \n')
    return included, excluded, newly_excluded

def save_molecule_to_png(molecule, output_path):
    img = Draw.MolToImage(molecule)
    img.save(output_path)

def dicts_w_info(directory, name_of_experiment, GIF_wanted):
    '''
    input: directory, the directory where the fragment pairs are stored
              name_of_experiment, the name of the experiment
    returns: pair_dict, a dictionary with all the information of the fragment pairs
                mz_dict, a dictionary with all the information of the m/z values
                mz_peaks, a list with all the m/z values
                result_overview, a dictionary with a summary of the results and the amount of simulations excluded at each function step
                counts_of_errors, a list of lists with the specific 'Error' and its count
    
    This function reads the files in the directory and extracts the information from the files.
    It calculates the m/z values for the fragments and stores them in a dictionary.
    It calculates the rate constants for the reaction paths and stores them in a dictionary.
    It plots the energy paths for the simulations to be included.
    It saves the energy paths in a folder called "Energy_paths" along with in each simulation folder.
    It saves the pair_dict and mz_dict as txt files in the directory.

    Excluding criteria: missing files from erroneus simulation runs
                        mass of fragment less than 3 g/mol
                        error determining bonds
                        no charge file found
                        not converged in CREST
                        error in xtb run
                        max energy point in last 2 points of the reaction path
    '''
    #print(rdBase.rdkitVersion)

    #change directory to input directory
    os.chdir(directory)

    pair_dict=defaultdict()

    #Coordinate file + pair no
    print('-----######-----Reading all files in the directory, finding fragment pairs and noting their simulation numbers')
    for root, dirs, files in os.walk(".", topdown=False):
        for name in files:
            if root not in pair_dict:
                if name.endswith('pair.xyz'):
                    print(root)
                    no=str(root)
                    pair_no=no.split('/')[-1]
                    #pair_no is only the numbers
                    pair_no=int(pair_no[1:])
                    pair_dict[pair_no]={'f_coord':os.path.join(root, name)}
    #Partial charges
    print('-----######-----Finding the charges files for the fragment pairs')
    for root, dirs, files in os.walk(".", topdown=False):
        for name in files:     
            if name.endswith('charges') and 'MSDIR' in root:
                no=str(root)
                pair_no=no.split('Pair_')[-1]
                pair_no=int(pair_no)
                if pair_no in pair_dict.keys():
                    path=os.path.join(root, name)
                    pair_dict[pair_no]['f_charge']=path
    #Fragment 1
    print('-----######-----Finding the mass files for the fragment pairs, fragment 1')
    for root, dirs, files in os.walk(".", topdown=False):
        for name in files:          
            if name=='mass' and root.endswith('1'):
                no=str(root)
                pair_no=no.split('/')[-1]
                pair_no=pair_no.split('p')[1]
                pair_no=pair_no.split('f')[0]
                pair_no=int(pair_no)
                #print(pair_no)
                if pair_no in pair_dict.keys():
                    path=os.path.join(root, name)
                    pair_dict[pair_no]['f_mass1']=os.path.join(root, name)
    #Fragment 2
    print('-----######-----Finding the mass files for the fragment pairs, fragment 2')
    for root, dirs, files in os.walk(".", topdown=False):
        for name in files: 
            if name=='mass' and root.endswith('2'):
                no=str(root)
                pair_no=no.split('/')[-1]
                pair_no=pair_no.split('p')[1]
                pair_no=pair_no.split('f')[0]
                pair_no=int(pair_no)
                if pair_no in pair_dict.keys():
                    pair_dict[pair_no]['f_mass2']=os.path.join(root, name)
    #Reaction path data - after xtb run
    print('-----######-----Finding the xtb reaction path files for the fragment pairs')
    for root, dirs, files in os.walk(".", topdown=False):
        for name in files:          
            if name.endswith('.out') and 'xtb_rxn' in root:
                no=str(root)
                pair_no=no.split('p')[-1]
                pair_no=pair_no.split('/')[0]
                pair_no=int(pair_no)
                if pair_no in pair_dict.keys():
                    path=os.path.join(root, name)
                    pair_dict[pair_no]['f_xtb_rxn']=path
                    pair_dict[pair_no]['xtb_path_coords']=os.path.join(root,'xtbpath.xyz')
                    

    for key in pair_dict.keys():
        pair_dict[key]['Exclude']=False

    print('\n-----######-----The following dictionary has been assembled before any exclusions:')
    print_nested_dict(pair_dict) #before any exclusions

    not_valid=[]
    for key in pair_dict.keys():
        if len(pair_dict[key])!=7: #Run initially 5, but 7 after xtb_rxn
            not_valid.append(key)
    not_valid  
    
    Fragment_pair=sorted(pair_dict.keys())
    print('\n-----######-----There are %d fragment pairs in the directory before exclusions.' % len(Fragment_pair))
    print('They can be found in the directories with the following numbers:')
    print(Fragment_pair)

    print('\n-----######-----Checking for errors in the fragment pairs due to missing files')
    '''Excluding criteria: missing files from erroneus simulation runs'''
    #Instead of deleting, we will add a key to the dictionary that contains the error and a key to exclude it
    for key in not_valid:
        error_messages = []

        if 'f_charge' not in pair_dict[key].keys():
            error_messages.append('No charges file found, not converged in CREST')
            print(f'{key} - Error: No charges file found, not converged in CREST')
        if 'f_mass1' not in pair_dict[key].keys():
            error_messages.append('No mass file found for fragment 1')
            print(f'{key} - Error: No mass file found for fragment 1')
        if 'f_mass2' not in pair_dict[key].keys():
            error_messages.append('No mass file found for fragment 2')
            print(f'{key} - Error: No mass file found for fragment 2')
        if 'f_xtb_rxn' not in pair_dict[key].keys():
            error_messages.append('No xtb reaction file found')
            print(f'{key} - Error: No xtb reaction file found')
        
        if error_messages:
            pair_dict[key]['Error'] = ', '.join(error_messages)
        else:
            pair_dict[key]['Error'] = 'Unknown error'
            print(f'{key} - Error: Unknown error')
        pair_dict[key]['Exclude'] = True

    #print_nested_dict(pair_dict) #after first exclusions due to missing files
                
    
    included, excluded, newly = check_excluded_included(pair_dict, len(Fragment_pair), 0, 'finding all files')
    included1, excluded1, newly1=included, excluded, newly
    #no_in_dictionary=0
    def mz(no_in_dictionary):
        '''
        input: no_in_dictionary, the number of the fragment pair in the dictionary
        returns: mz, the m/z values for the fragments
        
        This function reads the xyz file and the charges file for the fragment pair.
        It calculates the m/z values for the fragments and returns them.
        Furthermore it plots the fragments in the molecule and saves the image in the folder of the fragment pair.

        Excluding criteria: mass of fragment is less than 3 g/mol, 
                            error determining bonds, 
                            no charge file found/not converged in CREST
        '''
        key=Fragment_pair[no_in_dictionary]
        if pair_dict[key]['Exclude']==True:
            #print(f'{key} - {pair_dict[key]["Error"]}')
            return None
        coord_dir=pair_dict[key]['f_coord']
        '''if 'f_charge' not in pair_dict[key]:
            print(f'{key} - Error: No charges file found, not converged in CREST')
            return None'''
        charges_dir=pair_dict[key]['f_charge']
        mass0_dir=pair_dict[key]['f_mass1']
        mass1_dir=pair_dict[key]['f_mass2']
        charges_file=open(charges_dir,'r')
        charges=charges_file.readlines()
        if sum([1 for line in charges if line.strip()])==0:
            print('The coordinate file has no charge. ')
            raise ValueError('The coordinate file has no charge. CREST MSreact works with charged molecules. Remember to protonate the molecule.')
        molecule=rdkit.Chem.rdmolfiles.MolFromXYZFile(coord_dir)
        pair_dict[key]['molecule']=molecule
        mass_file0=open(mass0_dir,'r')
        mass0=mass_file0.readlines()
        mass0=float(mass0[0])
        pair_dict[key]['mass0']=mass0
        mass_file1=open(mass1_dir,'r')
        mass1=mass_file1.readlines()
        mass1=float(mass1[0])
        pair_dict[key]['mass1']=mass1
        if min(mass0,mass1)<3:
            pair_dict[key]['Error']='Mass of fragment is less than 3 g/mol'
            pair_dict[key]['Exclude']=True
            print(f'{key} - Error: Mass of fragment is less than 3 g/mol')
            return None

        # Getting the bonds 
        '''should the charge thing be removed? Would result in less errors due to "error determining bonds" '''
        try:
            rdDetermineBonds.DetermineBonds(molecule, charge=1)
            pair_dict[key]['bonds'] = molecule.GetBonds()
        except ValueError as e:
            '''Excluding criteria: error determining bonds'''
            pair_dict[key]['Error'] = 'Error determining bonds: ' + str(e)
            pair_dict[key]['Exclude'] = True
            print(f'{key} - Error: {str(e)}')
            return None

        #Find partial charges for each atom and print them
        AllChem.ComputeGasteigerCharges(molecule)
        for atom in molecule.GetAtoms():
            atom.GetProp('_GasteigerCharge')
            #print ( '%s | %s' % (atom.GetSymbol(), prop))
            
        #Overwrite all charges to zero and then overwrite them with the charges from the charges file
        i=0
        for atom in molecule.GetAtoms():
            atom.SetProp('_GasteigerCharge', str(0))
            atom.GetProp('_GasteigerCharge')
            atom.SetProp('_GasteigerCharge', str(charges[i]))
            i+=1

        output_dir = f'p{key}'
        #If we want to draw the fragments in the molecule, we can draw it with the following 
        dra=rdkit.Chem.Draw.MolToImage(molecule)
        dra.save(os.path.join(output_dir,f'{name_of_experiment}_p{key}.png'))
                    

        #Now we want to get the fragments:
        mol_frags = rdmolops.GetMolFrags(molecule, asMols = False)
        mol_frags
        #pair_dict[key]['fragments']=mol_frags
        pair_dict[key]['fragment1_ids']=[mol_frags[0]]
        pair_dict[key]['fragment2_ids']=[mol_frags[1]]
                
        #Calculating the charges for each fragment in the fragment pair
        charge0=0
        charge1=0
        chg_list0=[]
        chg_list1=[]
        for atom in molecule.GetAtoms():
            prop=atom.GetProp('_GasteigerCharge')
            id=atom.GetIdx()
            if id in mol_frags[0]:
                charge0+=float(atom.GetProp('_GasteigerCharge'))
                chg_list0.append(float(atom.GetProp('_GasteigerCharge')))
                
            elif id in mol_frags[1]:
                charge1+=float(atom.GetProp('_GasteigerCharge'))
                chg_list1.append(float(atom.GetProp('_GasteigerCharge')))

        #The charge should be 0 or 1, so we round it here and find m/z:
        charges_list=[]
        for charge in charges:
            charges_list.append(float(charge.strip()))
        charges_list
        pair_dict[key]['charges']=charges_list
        pair_dict[key]['charge0']=charge0
        pair_dict[key]['charge1']=charge1
        mz=[]
        #if charge0>charge1:
        mz.append(mass0/1)#round(charge0))
        #elif charge0<charge1:
        mz.append(mass1/1)#round(charge1))
        return mz


    #Now we loop through all of our fragment pairs, calculate the m/z values and store them in a list
    print('\n-----######-----Calculating m/z values for all valid fragment pairs')
    mzs=[]
    for i in range(len(Fragment_pair)):
        if mz(i)!= None and type(mz(i)[0])==float:
            print(f'{Fragment_pair[i]} - mass pairs:{mz(i)}')
            vals=mz(i)
            mzs.append(vals[0])
            mzs.append(vals[1])
            #handling valueerror:

    included, excluded, newly = check_excluded_included(pair_dict, included, excluded, 'finding m/z values')
    included2, excluded2, newly2=included, excluded, newly
    
    mz_peaks=np.unique(mzs)


    '''#load test file
    dir='./test.xyz'
    molecule=rdkit.Chem.rdmolfiles.MolFromXYZFile(dir)
    rdDetermineBonds.DetermineBonds(molecule,charge=0)
    dra=rdkit.Chem.Draw.MolToImage(molecule)
    dra.save('molecule.png')''' #Not sure whether this will be important.
    #key=Fragment_pair[22]

    #moving on to finding the correct intensities:
    def find_rate_const_n_plot(no_in_dictionary, GIF_wanted=False, name_of_experiment=name_of_experiment):
        '''
        input: no_in_dictionary, the number of the fragment pair in the dictionary [0:len(pair_dict)]
               GIF_wanted, a boolean that determines whether a gif of the reaction path should be created
        returns: k, the rate constant for the reaction path as float

        This function reads the xtb output file and extracts the valid energy values for the reaction path. 
        Thereafter it calculates the rate constant for the reaction path. 
        Furthermore, it plots the energy path for the simulations to be included.
        
        Lastly, if specified GIF_wanted input is 'True' it creates a gif of the reaction path for the fragment pair based on the xtb_path.xyz file.
        This however adds additional computer resources.
        
        Excluding criteria: Error in xtb run,
                            The max energy in energy path is within last 2 points
        '''
        key=Fragment_pair[no_in_dictionary]
        if pair_dict[key]['Exclude']==True:
            #print(f'{key} is excluded due to error: {pair_dict[key]["Error"]}')
            return None
        rxn_path_dir=pair_dict[key]['f_xtb_rxn']
        #readlines after the line that contains "point     drms     energy pmode ovlp pmode grad"
        file=open(rxn_path_dir,'r')
        lines=file.readlines()
        #Update error message, if need be:
        for line in lines:
            if '################################################################' in line:
                print(f'{key} contains an error in the xtb run')
                if pair_dict[key].get('Error'):
                    pair_dict[key]['Error']+=', Error in xtb run'
                else:
                    pair_dict[key]['Error']='Error in xtb run'
                    pair_dict[key]['Exclude']=True
                return None
        i=0
        for line in lines:
            if 'point     drms     energy pmode ovlp pmode grad' in line:
                break
            i+=1
            if 'reactant product RMSD' in line:
                rmsd_line=line
        i+=1 #The first line of the intensities
        #now we need to find the last line of the intensities, which is the last line after 'i' that contains more than a " "
        j=i
        while len(lines[j])>1:
            j+=1
        rmsd_line=rmsd_line.split()
        pair_dict[key]['rmsd_drug_fragments']=rmsd_line[-1]
        #i,j #last line of intensities
        #Now save these lines in a pandas dataframe
        rxnpath_data=[]
        for line in range(i,j):
            rxnpath_data.append((lines[line].split()))
        rxnpath_data = [[float(item) for item in line.split()] for line in lines[i:j]]
        #print(rxnpath_data)
        #Now we have the data in a list of lists, we can convert it to a pandas dataframe
        df=pd.DataFrame(rxnpath_data, columns=['Point', 'drms', 'Energy', 'pmode_ovlp', 'pmode_grad'])
        #pair_dict[key]['xtb_output']=df
        pair_dict[key]['transition_path']=df['Energy'].tolist()
        if max(pair_dict[key]['transition_path']) in df['Energy'][-2:].values:
            print(f'{key} is excluded due to error: max point is in the last 2 points of the reaction path')
            if pair_dict[key].get('Error'):
                pair_dict[key]['Error']+=', max energy point in last 2 points'
            else:
                pair_dict[key]['Error']='max energy point in last 2 points'
                pair_dict[key]['Exclude']=True
            return None
        pair_dict[key]['transition_energy']=max(pair_dict[key]['transition_path'])-(pair_dict[key]['transition_path'][0])
        #print_nested_dict(pair_dict)
        df.plot(x='Point', y='Energy')
        plt.xlabel('Point')
        plt.ylabel('Energy')
        #plt.xlim(0,path_data['Point'][len(path_data)-1])
        #plt.ylim(min(path_data['Energy']),max(path_data['Energy']))
        plt.title('Energy of the reaction path')
        plt.savefig(f'p{key}/'+f'Energy_path_{key}.png')
        plt.close()
        
        '''Adding the transition path as a gif'''
        if GIF_wanted==True:
            output_dir = f'p{key}'
            try:
                os.mkdir(os.path.join(output_dir,'Reaction_paths'))
            except:
                pass
            output_dir = f'p{key}/Reaction_paths'
            
            #might as well also save each of the coordinate files in the xtb_path.xyz file:
            with open(pair_dict[key]['xtb_path_coords'],'r') as f:
                lines=f.readlines()
            trajectories = {}
            current_trajectory = []
            i = 0

            for line in lines:
                stripped_line = line.strip()
                if stripped_line.isdigit():  # Identifies the start of a new trajectory
                    if current_trajectory:
                        trajectories[i] = current_trajectory
                        current_trajectory = []
                    current_trajectory.append(line)
                    i += 1
                else:
                    current_trajectory.append(line)

            # Add the last trajectory if there is one
            if current_trajectory:
                trajectories[i] = current_trajectory

            # Store each trajectory in a mol object and save to image files
            images = []
            for idx, trajectory in trajectories.items():
                # Write the current trajectory to a temporary XYZ file
                xyz_file_path = os.path.join(output_dir, f'trajectory_{idx}.xyz')
                with open(xyz_file_path, 'w') as xyz_file:
                    xyz_file.write(''.join(trajectory))
                # Read the molecule from the XYZ file
                mol = Chem.MolFromXYZFile(xyz_file_path)
                try:
                    rdDetermineBonds.DetermineBonds(mol, charge=1)
                except ValueError as e:
                    print(f'Error determining bonds with charge=1: {e}')
                    try:
                        rdDetermineBonds.DetermineBonds(mol, charge=0)
                    except ValueError as e:
                        print(f'Error determining bonds with charge=0: {e}')
                        print('Skipping bond determination due to error')

                if mol is None:
                    continue  # Skip if the molecule could not be created
                # Define the output file path
                output_path = os.path.join(output_dir, f'trajectory_{idx}.png')
                save_molecule_to_png(mol, output_path)
                # Append the output file path to the images list
                images.append(output_path)
            
            # Create a GIF from the images
            gif_name = os.path.join(output_dir, f'transpath_{name_of_experiment}_{key}.gif')
            with imageio.get_writer(gif_name, mode='I', duration=0.5) as writer:
                for image_path in images:
                    image = imageio.imread(image_path)
                    writer.append_data(image)
            try:
                os.mkdir('4_transitions_paths_gifs')
            except:
                pass
            shutil.copy(gif_name, '4_transitions_paths_gifs')
            

        #In order not to redo everything all the time, we can save the pair_dict as a pickle file
        diff=pair_dict[key]['transition_energy']
        diff=diff*4.184 #convert from kcal/mol to J/mol
        k=np.exp(-diff/(8.314*.298)) #rate constant which is proportional to the abundance
        return k #also saves the energy_path plot


    #xtb_error=[]
    #sum_ks=0
    print('-----######-----Finding rate constants for all fragment pairs and plotting the energy paths in the valid simulations')
    for no in range(len(Fragment_pair)):
        if pair_dict[Fragment_pair[no]]['Exclude']==False:
            pair_dict[Fragment_pair[no]]['rate_const']=find_rate_const_n_plot(no, GIF_wanted, name_of_experiment)
            if pair_dict[Fragment_pair[no]]['Exclude']==False:
                rate_k=pair_dict[Fragment_pair[no]]['rate_const']
                print(f'{Fragment_pair[no]} - The rate constant for the reaction path is: {rate_k} J/mol')
    print('-----######-----The following fragment pairs are excluded:')
    for no in range(len(Fragment_pair)):
        if pair_dict[Fragment_pair[no]]['Exclude']==True:
            print(f'{Fragment_pair[no]} - {pair_dict[Fragment_pair[no]]["Error"]}')
    
    included, excluded, newly = check_excluded_included(pair_dict, included, excluded, 'finding rate constants')
    included3, excluded3, newly3=included, excluded, newly
    

    print('\n-----######-----All rate constants have been calculated and the energy paths have been plotted')
    print('\n-----######-----All energy paths are copied to the folder "Energy_paths" in the current directory')
    #Make a directory called Energy_paths and save all plots from the energy paths in there
    try:
        os.mkdir(f'2_Energypaths_{name_of_experiment}')
    except:
        pass
    #copy all plots from /*/p*/Energy_path_*.png to Energy_paths/
    for root, dirs, files in os.walk(".", topdown=False):
        for name in files:
            if name.startswith('Energy_path_') and f'Energypaths_{name_of_experiment}' not in root:
                shutil.copy(os.path.join(root, name), f'2_Energypaths_{name_of_experiment}')
    

    #Make a dictionary with the mz_peaks and the corressponding information from pair_dict
    print('\n-----######-----Sorting all the m/z peaks and their corresponding information in mz_dict')
    mz_dict={}
    for mass in mz_peaks:
        '''rate constants in mz_dict are only based on valid simulations, not excluded ones'''
        mz_dict[mass]={}
        #find_keys_by_value(pair_dict, mass)
        matching_keys = [key for key, value in pair_dict.items() if value.get('mass0') == mass or value.get('mass1') == mass]
        mz_dict[mass]['simulations']=matching_keys
        mz_dict[mass]['number_of_simulations']=len(matching_keys)
        mz_dict[mass]['rate_constants']=[]
        #mz_dict[mass]['prob_contribs']=[]
        #mz_dict[mass]['log_probabilities']=[]
        mz_dict[mass]['transition_energies']=[]
        for key in matching_keys:
            if pair_dict[key]['Exclude']==True: 
                #Excluded simulations are not used to find the rate constants
                #They are however added for overview
                if mz_dict[mass].get('Excluded'):
                    mz_dict[mass]['Excluded'].append(key)
                    #remove key from 'simulations'
                    mz_dict[mass]['Error from excluded']=pair_dict[key]['Error']
                else:
                    mz_dict[mass]['Excluded']=[key]
            else:
                #Only valid simulations are used to find the rate constants
                mz_dict[mass]['rate_constants'].append(pair_dict[key]['rate_const'])
                #mz_dict[mass]['prob_contribs'].append(pair_dict[key]['prob_contrib'])
                #mz_dict[mass]['log_probabilities'].append(pair_dict[key]['log_probability'])
                mz_dict[mass]['transition_energies'].append(pair_dict[key]['transition_energy'])
        mz_dict[mass]['rate_constant_UQ']=np.unique(mz_dict[mass]['rate_constants'])
        #mz_dict[mass]['prob_contrib']=np.unique(mz_dict[mass]['prob_contribs']) #based on all rate_constants and not only unique ones - therefor this is done beneath
        #mz_dict[mass]['log_probability']=np.unique(mz_dict[mass]['log_probabilities'])
        mz_dict[mass]['transition_energy_UQ']=np.unique(mz_dict[mass]['transition_energies'])
    
    for mass in mz_peaks:
        if mz_dict[mass].get('Excluded'):
            if key in mz_dict[mass]['Excluded']:
                mz_dict[mass]['simulations'].remove(key)
                

    '''Excluding criteria from mz_dict: no valid simulations of same mass'''
    count_del=0
    for mass in mz_peaks:
        if len(mz_dict[mass]['rate_constant_UQ'])==0:
            del mz_dict[mass]
            print(f'{mass} has been removed from the dictionary due to no valid simulations')
            count_del+=1
    
    print(f'\n-----######----- In total {count_del} m/z peaks were removed from the dictionary due to only invalid simulations')

    mz_peaks=np.array(list(mz_dict.keys()))
    #mz_peaks
    #now we sum up the prob_contribs
    sum_rks=0 #all found rate constants
    for mass in mz_peaks:
        print(mass)
        for val in mz_dict[mass]['rate_constant_UQ']:
            sum_rks+=val
    sum_rks
   
    #Adding probabilities based on summed rate constants
    for mass in mz_peaks:
        mz_dict[mass]['Sum_UQ_rate_constants']=np.sum(mz_dict[mass]['rate_constant_UQ'])
        mz_dict[mass]['Prob_contrib']=mz_dict[mass]['Sum_UQ_rate_constants']/sum_rks #best one
        #mz_dict[mass]['log_prob']=np.log(np.sum(mz_dict[mass]['rate_constant_UQ']))/np.log(sum_rks)


    #Save txt files with print_nested_dict for each of the dictionaries in directory:
    # save file with print_nested_dict(pair_dict)
    #create a folder for the experiment output
    print('\n-----OUTPUT-----The resulting dictionaries are saved as txt files in the MSreact directory')
    try:
        os.mkdir(f'1_output_{name_of_experiment}')
    except:
        pass

    with open(f'1_output_{name_of_experiment}/pairdict.txt', 'w') as file:
        file.write(nested_dict_layout(pair_dict))
    
    # save file with print_nested_dict(mz_dict)
    with open(f'1_output_{name_of_experiment}/mzdict.txt', 'w') as file:
        file.write(nested_dict_layout(mz_dict))
    
    print('\n-----######-----The dictionary containing fully processed simulations (pair_dict) is here:')
    print_nested_dict(pair_dict)
    print('\n-----######-----The dictionary containing the m/z values and their corresponding information (mz_dict) is here:')
    print_nested_dict(mz_dict)
    
    print('\n-----RESUME OF ACTIONS-----')
    #total amount of simulations:
    print(f'The total amount of simulations with fragment pairs in {directory} is: {len(Fragment_pair)} \n')
    #all exclusions
    print(f'The amount of simulations excluded during finding files: {newly1}')
    print(f'The amount of simulations excluded during finding m/z values: {newly2}')
    print(f'The amount of simulations excluded during finding rate constants: {newly3} \n')
    print(f'The amount of unique mass peaks found in the simulations is: {count_del+len(mz_peaks)}')    
    print(f'The amount of mass peaks removed from the dictionary due to no valid simulations: {count_del} \n')
    
    print('The specific errors found in the simulations are: ')
    #count the amount of each error in 'Error' and state as result:
    errors=[]
    for key in pair_dict.keys():
        if pair_dict[key].get('Error'):
            errors.append(pair_dict[key]['Error'])
    errors=set(errors)
    print(f'The errors found in the simulations are: {errors} \n')
    
    counts_of_errors={}
    counts_of_errors['Types of errors']={}
    for error in errors:
        counts_of_errors['Types of errors'][error]=0
        for key in pair_dict.keys():
            if pair_dict[key].get('Error')==error:
                counts_of_errors['Types of errors'][error]+=1
    print(f'The amount of each error found in the simulations is:')
    print_nested_dict(counts_of_errors)
    
    
    print(f'\n-----RESULTS: {name_of_experiment}-----')
    print(f'The total amount of m/z peaks found in the simulations is: {len(mz_peaks)}')
    print(f'The peaks are at m/Z: {mz_peaks} \n')
    result_overview={'Simulation_overview':{'Amount of simulations': len(Fragment_pair),'Ex_finding_files': newly1, 'Ex_finding_mz_values': newly2, 'Ex_finding_rate_constants': newly3}, 'mZ_peaks': {'Total_amount': len(mz_peaks), 'Ex_Peaks': count_del}}
    

    return pair_dict, mz_dict, mz_peaks, result_overview, counts_of_errors

#print_nested_dict(mz_dict)
#print_nested_dict(pair_dict)

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot_mz(mz_dict, experimental_datafile, name_of_experiment, name_of_drug, mass_of_drug):
    try:
        os.mkdir(f'3_plots_{name_of_experiment}')
    except FileExistsError:
        pass
    plt.clf()
    x_mz = []
    #y_log_prob = []
    y_k_log = []
    y_k = []
    y_prob_contrib = []
    y_avg_trans = []
    y_avg_trans_int = []
    y_amount_sim = []
    y_amount_sim_UQ=[]
    
    amount_sim=0
    amount_sim_UQ=0
    for mass in mz_dict.keys():
        amount_sim+=len(mz_dict[mass]['rate_constants'])
        amount_sim_UQ+=len(mz_dict[mass]['rate_constant_UQ'])
    
    
    for mass in mz_dict.keys():
        x_mz.append(mass)
        #y_log_prob.append(mz_dict[mass]['log_prob'])  # log of sum of rate constants/log(of all rate constants)
        y_prob_contrib.append(mz_dict[mass]['Prob_contrib'])  # contributed unique rate constants/of all unique rate constants
        y_k.append(mz_dict[mass]['Sum_UQ_rate_constants'])  # sum of all rate constants
        y_k_log.append(-np.log(mz_dict[mass]['Sum_UQ_rate_constants']))
        y_avg_trans.append(np.average(mz_dict[mass]['transition_energy_UQ']))
        y_avg_trans_int.append(1 / y_avg_trans[-1])
        y_amount_sim.append(len(mz_dict[mass]['rate_constants'])/amount_sim)
        y_amount_sim_UQ.append(len(mz_dict[mass]['rate_constant_UQ'])/amount_sim_UQ)

    if experimental_datafile is not None:
        drug_obs = pd.read_csv(experimental_datafile)
        drug_obs_x = drug_obs['M/Z_fragment'].tolist()
        CE = drug_obs['CE']  # relatable to activation energy
        y_int_obs = 1 / CE  # intensity is inversely proportional to collision energy
        k_obs = np.exp(-CE / (8.314 * .298))
        k_obs_log = -np.log(k_obs)

    def create_stacked_barplot(ax, x, y, label, color):
        ax.bar(x, y, color=color, alpha=0.7, label=label)
        ax.set_xlabel('m/z')
        ax.set_xlim(0, mass_of_drug)
        ax.legend()

    # Plot with log rate constants
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))
    if experimental_datafile is not None:
        create_stacked_barplot(ax1, drug_obs_x, k_obs_log, 'Drug Observations', 'red')
    create_stacked_barplot(ax2, x_mz, y_k_log, name_of_experiment, 'green')
    ax1.set_ylabel('Intensity (log of rate constant)')
    ax1.set_title(f'In silico Mass spectrum for {name_of_drug}')
    ax2.set_ylabel('Intensity (log of rate constant)')
    plt.tight_layout()
    plt.savefig(f'3_plots_{name_of_experiment}/mz_{name_of_experiment}_log_rate.png')
    plt.close()
    plt.clf()
    
    # Plot with rate constants
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))
    if experimental_datafile is not None:
        create_stacked_barplot(ax1, drug_obs_x, k_obs, 'Drug Observations', 'red')
    create_stacked_barplot(ax2, x_mz, y_k, name_of_experiment, 'green')
    ax1.set_ylabel('Intensity (rate constant)')
    ax1.set_title(f'In silico Mass spectrum for {name_of_drug}')
    ax2.set_ylabel('Intensity (rate constant)')
    plt.tight_layout()
    plt.savefig(f'3_plots_{name_of_experiment}/mz_{name_of_experiment}_rate.png')
    plt.close()
    plt.clf()

    # Plot with inverse transition/collision energies
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))
    if experimental_datafile is not None:
        create_stacked_barplot(ax1, drug_obs_x, y_int_obs, 'Drug Observations', 'red')
    create_stacked_barplot(ax2, x_mz, y_avg_trans_int, name_of_experiment, 'green')
    ax1.set_ylabel('Intensity (Inverse collision energy)')
    ax1.set_title(f'In silico Mass spectrum for {name_of_drug}')
    ax2.set_ylabel('Intensity (Inverse transition energy)')
    plt.tight_layout()
    plt.savefig(f'3_plots_{name_of_experiment}/mz_{name_of_experiment}_inverse_trans.png')
    plt.close()
    plt.clf()

    # Plot with probability contributions
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))
    if experimental_datafile is not None:
        create_stacked_barplot(ax1, drug_obs_x, y_int_obs, 'Drug Observations', 'red')
    create_stacked_barplot(ax2, x_mz, y_prob_contrib, name_of_experiment, 'green')
    ax1.set_ylabel('Intensity (Inverse collision energy)')
    ax1.set_title(f'In silico Mass spectrum for  {name_of_drug}')
    ax2.set_ylabel('Intensity (Probability Contribution)')
    plt.tight_layout()
    plt.savefig(f'3_plots_{name_of_experiment}/mz_{name_of_experiment}_prob_contrib.png')
    plt.close()
    plt.clf()
    
    # Normalize based on the amount of simulations in each m/z peak
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))
    if experimental_datafile is not None:
        create_stacked_barplot(ax1, drug_obs_x, y_int_obs, 'Drug Observations', 'red')
    create_stacked_barplot(ax2, x_mz, y_amount_sim, name_of_experiment, 'green')
    ax1.set_ylabel('Intensity (Inverse collision energy)')
    ax1.set_title(f'In silico Mass spectrum for {name_of_drug}')
    ax2.set_ylabel('Intensity (normalized amount of valid rate constants)')
    plt.tight_layout()
    plt.savefig(f'3_plots_{name_of_experiment}/mz_{name_of_experiment}_amount_rates.png')
    plt.close()
    plt.clf()
    
    # Normalize based on the amount of simulations in each m/z peak
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))
    if experimental_datafile is not None:
        create_stacked_barplot(ax1, drug_obs_x, y_int_obs, 'Drug Observations', 'red')
    create_stacked_barplot(ax2, x_mz, y_amount_sim_UQ, name_of_experiment, 'green')
    ax1.set_ylabel('Intensity (Inverse collision energy)')
    ax1.set_title(f'In silico Mass spectrum for {name_of_drug}')
    ax2.set_ylabel('Intensity (normalized amount of valid unique rate constants)')
    plt.tight_layout()
    plt.savefig(f'3_plots_{name_of_experiment}/mz_{name_of_experiment}_amount_UQ_rates.png')
    plt.close()
    plt.clf()

    return

import glob
def copy_png_files(root, name_of_experiment, name_of_drug):
    source_dir = os.path.join(root, f'3_plots_{name_of_experiment}')
    dest_dir = f'/work3/ishof/SK/CREST/2_Supervision/Results/{name_of_drug}'
    
    # Ensure destination directory exists
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)
    
    # Use glob to find all png files in the source directory
    png_files = glob.glob(os.path.join(source_dir, f'mz_{name_of_experiment}*.png'))
    
    # Copy each file to the destination directory
    for file in png_files:
        shutil.copy(file, dest_dir)



