import os
import rdkit
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdmolfiles
from rdkit.Chem.rdmolfiles import MolFromXYZFile
from collections import defaultdict
from rdkit.Chem import Descriptors
from rdkit.Chem import rdmolops
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem import Draw
from rdkit.Chem.Draw import MolsToGridImage
from rdkit.Chem import rdDetermineBonds
from rdkit import rdBase
#from mpmath import mpf, mpc, mp

#from scipy.constants import Boltzmann
print(rdBase.rdkitVersion)

#Run from SK/CREST
print('Please enter the name of the drug you want to analyze:')
name_of_drug=input().lower()


pair_dict=defaultdict()

#Coordinate file + pair no
for root, dirs, files in os.walk(".", topdown=False):
    for name in files:
        if root not in pair_dict:
            if name.endswith('pair.xyz'):
                if name_of_drug in root.lower():
                    #print(root)
                    no=str(root)
                    pair_no=no.split('/')[-1]
                    #pair_no is only the numbers
                    pair_no=int(pair_no[1:])
                    pair_dict[pair_no]={'f_coord':os.path.join(root, name)}
#Partial charges
for root, dirs, files in os.walk(".", topdown=False):
    for name in files:     
        if name.endswith('charges') and 'xtb' not in root:
            if name_of_drug in root.lower():
                no=str(root)
                pair_no=no.split('Pair_')[-1]
                pair_no=int(pair_no)
                if pair_no in pair_dict.keys():
                    path=os.path.join(root, name)
                    pair_dict[pair_no]['f_charge']=path
#Fragment 1
for root, dirs, files in os.walk(".", topdown=False):
    for name in files:          
        if name=='mass' and root.endswith('1'):
            if name_of_drug in root.lower():
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
for root, dirs, files in os.walk(".", topdown=False):
    for name in files: 
        if name=='mass' and root.endswith('2'):
            if name_of_drug in root.lower():
                no=str(root)
                pair_no=no.split('/')[-1]
                pair_no=pair_no.split('p')[1]
                pair_no=pair_no.split('f')[0]
                pair_no=int(pair_no)
                if pair_no in pair_dict.keys():
                    pair_dict[pair_no]['f_mass2']=os.path.join(root, name)
#Reaction path data - after xtb run
for root, dirs, files in os.walk(".", topdown=False):
    for name in files:          
        if name.endswith('.out') and 'xtb_rxn' in root:
            if name_of_drug in root.lower():
                no=str(root)
                pair_no=no.split('p')[-1]
                pair_no=pair_no.split('/')[0]
                pair_no=int(pair_no)
                if pair_no in pair_dict.keys():
                    path=os.path.join(root, name)
                    pair_dict[pair_no]['f_xtb_rxn']=path

for key in pair_dict.keys():
    pair_dict[key]['Exclude']=False

def print_nested_dict(d, indent=0):
    for key, value in d.items():
        if isinstance(value, dict):
            print(' ' * indent + str(key) + ':')
            print_nested_dict(value, indent + 4)
        else:
            print(' ' * indent + str(key) + ': ' + str(value))

print_nested_dict(pair_dict)

not_valid=[]
for key in pair_dict.keys():
    if len(pair_dict[key])!=6: #Run initially
        not_valid.append(key)
not_valid
#for key in not_valid:
#    del pair_dict[key]    

#Instead of deleting, we will add a key to the dictionary that contains the error and a key to exclude it
for key in not_valid:
    if 'f_charge' not in pair_dict[key]:
        if pair_dict[key].get('Error'):
            pair_dict[key]['Error']+=', No charges file found, not converged in CREST'
        else:
            pair_dict[key]['Error']='No charges file found, not converged in CREST'
    if 'f_mass1' not in pair_dict[key]:
        if pair_dict[key].get('Error'):
            pair_dict[key]['Error']+=', No mass file found for fragment 1'
        else: 
            pair_dict[key]['Error']='No mass file found for fragment 1'
    if 'f_mass2' not in pair_dict[key]:
        if pair_dict[key].get('Error'):
            pair_dict[key]['Error']+=', No mass file found for fragment 2'
        else:
            pair_dict[key]['Error']='No mass file found for fragment 2'
    if 'f_xtb_rxn' not in pair_dict[key]:
        if pair_dict[key].get('Error'):
            pair_dict[key]['Error']+=', No xtb reaction file found'
        else:
            pair_dict[key]['Error']='No xtb reaction file found'
    else:
        pair_dict[key]['Error']='Unknown error'
    pair_dict[key]['Exclude']=True

print_nested_dict(pair_dict)
             
Fragment_pair=sorted(pair_dict.keys())
print('There are %d fragment pairs in the directory.' % len(Fragment_pair))
print('They can be found in the directories with the following numbers:')
Fragment_pair


def mz(no_in_dictionary):
    coord_dir=pair_dict[Fragment_pair[no_in_dictionary]]['f_coord']
    if 'f_charge' not in pair_dict[Fragment_pair[no_in_dictionary]]:
        return "Error: No charges file found, not converged in CREST"
    charges_dir=pair_dict[Fragment_pair[no_in_dictionary]]['f_charge']
    mass0_dir=pair_dict[Fragment_pair[no_in_dictionary]]['f_mass1']
    mass1_dir=pair_dict[Fragment_pair[no_in_dictionary]]['f_mass2']
    charges_file=open(charges_dir,'r')
    charges=charges_file.readlines()
    molecule=rdkit.Chem.rdmolfiles.MolFromXYZFile(coord_dir)
    pair_dict[Fragment_pair[no_in_dictionary]]['molecule']=molecule
    mass_file0=open(mass0_dir,'r')
    mass0=mass_file0.readlines()
    mass0=float(mass0[0])
    pair_dict[Fragment_pair[no_in_dictionary]]['mass0']=mass0
    mass_file1=open(mass1_dir,'r')
    mass1=mass_file1.readlines()
    mass1=float(mass1[0])
    pair_dict[Fragment_pair[no_in_dictionary]]['mass1']=mass1
    if min(mass0,mass1)<3:
        pair_dict[Fragment_pair[no_in_dictionary]]['Error']='Mass of fragment is less than 3 g/mol'
        pair_dict[Fragment_pair[no_in_dictionary]]['Exclude']=True
        return "Error: Only dihydrogen knocked off."

    #Getting the bonds
    rdDetermineBonds.DetermineBonds(molecule,charge=-1)
    pair_dict[Fragment_pair[no_in_dictionary]]['bonds']=molecule.GetBonds()

    #Find partial charges for each atom and print them
    AllChem.ComputeGasteigerCharges(molecule)
    for atom in molecule.GetAtoms():
        prop = atom.GetProp('_GasteigerCharge')
        #print ( '%s | %s' % (atom.GetSymbol(), prop))
        
    #Overwrite all charges to zero and then overwrite them with the charges from the charges file
    i=0
    for atom in molecule.GetAtoms():
        atom.SetProp('_GasteigerCharge', str(0))
        prop = atom.GetProp('_GasteigerCharge')
        prop = atom.SetProp('_GasteigerCharge', str(charges[i]))
        i+=1


    #If we want to draw the fragments in the molecule, we can draw it with the following 
    #dra=rdkit.Chem.Draw.MolToImage(molecule)
    #dra.save('molecule.png')

    #Now we want to get the fragments:
    mol_frags = rdmolops.GetMolFrags(molecule, asMols = False)
    mol_frags
    #pair_dict[Fragment_pair[no_in_dictionary]]['fragments']=mol_frags
    pair_dict[Fragment_pair[no_in_dictionary]]['fragment1_ids']=[mol_frags[0]]
    pair_dict[Fragment_pair[no_in_dictionary]]['fragment2_ids']=[mol_frags[1]]
            
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
    pair_dict[Fragment_pair[no_in_dictionary]]['charges']=charges_list
    pair_dict[Fragment_pair[no_in_dictionary]]['charge0']=charge0
    pair_dict[Fragment_pair[no_in_dictionary]]['charge1']=charge1
    mz=[]
    #if charge0>charge1:
    mz.append(mass0/1)#round(charge0))
    #elif charge0<charge1:
    mz.append(mass1/1)#round(charge1))
    return mz


#Now we loop through all of our fragment pairs, calculate the m/z values and store them in a list
mzs=[]
for i in range(len(Fragment_pair)):
    print(mz(i))
    if type(mz(i)[0])==float:
        vals=mz(i)
        mzs.append(vals[0])
        mzs.append(vals[1])
mzs
mz_peaks=np.unique(mzs)
mz_peaks

print_nested_dict(pair_dict)

#load the verapamil data from the experiment
drug_obs=pd.read_csv('verapamil_exp.csv')
drug_obs_x=drug_obs['M/Z_fragment']
drug_obs_y=np.ones(len(drug_obs_x))
drug_obs_y=drug_obs_y*0.5
drug_obs_x
drug_obs_y

#Make a pillar diagram ###Not so pretty right now.
import matplotlib.pyplot as plt
import numpy as np
#import spectrum_utils.plot as sup
#import spectrum_utils.spectrum as sus
y_vals=np.ones(len(mz_peaks))
#X-axis
plt.xlabel('m/z')
plt.ylabel('Intensity')
plt.title('m/Z values for Verapamil')
plt.xlim(0, 450) #Adjust for other molecules
plt.bar(mz_peaks, y_vals, color='blue')
plt.bar(drug_obs_x, drug_obs_y, color='red')
#save as png
plt.savefig('mz_new.png')


'''#load test file
dir='./test.xyz'
molecule=rdkit.Chem.rdmolfiles.MolFromXYZFile(dir)
rdDetermineBonds.DetermineBonds(molecule,charge=0)
dra=rdkit.Chem.Draw.MolToImage(molecule)
dra.save('molecule.png')''' #Not sure whether this will be important.
key=Fragment_pair[0]

#moving on to finding the correct intensities:
sum_abundance=0
def rel_abund(no_in_dictionary):
    key=Fragment_pair[no_in_dictionary]
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
    if pair_dict[key]['Exclude']==True:
        print(f'{key} is excluded due to error {pair_dict[key]["Error"]}')
        return None
    i=0
    for line in lines:
        if 'point     drms     energy pmode ovlp pmode grad' in line:
            break
        i+=1
    i+=1 #The first line of the intensities
    #now we need to find the last line of the intensities, which is the last line after 'i' that contains more than a " "
    j=i
    while len(lines[j])>1:
        j+=1

    i,j #last line of intensities
    #Now save these lines in a pandas dataframe
    rxnpath_data=[]
    for line in range(i,j):
        rxnpath_data.append((lines[line].split()))
    rxnpath_data = [[float(item) for item in line.split()] for line in lines[i:j]]
    #print(rxnpath_data)
    #Now we have the data in a list of lists, we can convert it to a pandas dataframe
    df=pd.DataFrame(rxnpath_data, columns=['Point', 'drms', 'Energy', 'pmode_ovlp', 'pmode_grad'])
    #pair_dict[key]['xtb_output']=df
    pair_dict[key]['transition_path']=df['Energy']
    pair_dict[key]['transition_energy']=max(pair_dict[key]['transition_path'])-(pair_dict[key]['transition_path'][0])
    #print_nested_dict(pair_dict)
    df.plot(x='Point', y='Energy')
    plt.xlabel('Point')
    plt.ylabel('Energy')
    #plt.xlim(0,path_data['Point'][len(path_data)-1])
    #plt.ylim(min(path_data['Energy']),max(path_data['Energy']))
    plt.title('Energy of the reaction path')
    plt.savefig(f'p{key}/'+f'Energy_path_{key}.png')

    #In order not to redo everything all the time, we can save the pair_dict as a pickle file
    diff=pair_dict[key]['transition_energy']
    diff=diff*4.184 #convert from kcal/mol to J/mol
    k=np.exp(-diff/(8.314*.298)) #rate constant which is proportional to the abundance
    return k


pair1_dict=pair_dict
xtb_error=[]
sum_abundance=0
sum_inv_k=0
'''Needs to be run twice at the moment in order to work?'''
for i in range(len(Fragment_pair)):
    if pair_dict[Fragment_pair[i]]['Exclude']==False or rel_abund(i)!=None:
        sum_abundance+=rel_abund(i)
        rel_abundance=rel_abund(i)
        pair_dict[Fragment_pair[i]]['rate_const']=rel_abundance
        inverse_k=-np.log(pair_dict[Fragment_pair[i]]['rate_const'])
        pair_dict[Fragment_pair[i]]['inverse_k']=inverse_k
        sum_inv_k+=inverse_k
    else:
        xtb_error.append(Fragment_pair[i])
xtb_error
sum_abundance
sum_inv_k

#print all relative abundances
for i in range(len(Fragment_pair)):
    if pair_dict[Fragment_pair[i]]['Exclude']==False:
        probability=pair_dict[Fragment_pair[i]]['rate_const']/sum_abundance
        pair_dict[Fragment_pair[i]]['probability']=probability
        inv_prob=pair_dict[Fragment_pair[i]]['inverse_k']/sum_inv_k #fix with summing the rate constant - then inverse
        pair_dict[Fragment_pair[i]]['inv_probability']=inv_prob
        

print_nested_dict(pair_dict)


pair_dict
plt.plot(pair_dict) 


#Scheme with whatever value you want to show
for key in pair_dict.keys():
    #print(key)
    try:
        transition_energy=pair_dict[key]['transition_energy']
        mass0=pair_dict[key]['mass0']
        mass1=pair_dict[key]['mass1']
        rate_const=pair_dict[key]['rate_const']
        inv_k=pair_dict[key]['inverse_k']
        prob=pair_dict[key]['probability']
        inv_prob=pair_dict[key]['inv_probability']
        print(f'{key} \t {mass0:.3f} \t {mass1:.3f} \t {transition_energy:.3f} \t {rate_const:.5e} \t {inv_k:.3f} \t {prob:.5e} \t {inv_prob:.3f}')
    except:
        pass

print_nested_dict(pair_dict)

#Make a dictionary with the mz_peaks and the corressponding information from pair_dict
mz_dict={}
for mass in mz_peaks:
    mz_dict[mass]={}
    #find_keys_by_value(pair_dict, mass)
    matching_keys = [key for key, value in pair_dict.items() if value.get('mass0') == mass or value.get('mass1') == mass]
    mz_dict[mass]['simulations']=matching_keys
    mz_dict[mass]['number_of_simulations']=len(matching_keys)
    for key in matching_keys:
        if pair_dict[key]['Exclude']==True:
            break
        else:
            mz_dict[mass]['rate_constants']=pair_dict[key]['rate_const']
            mz_dict[mass]['inverse_ks']=pair_dict[key]['inverse_k']
            mz_dict[mass]['probabilities']=pair_dict[key]['probability']
            mz_dict[mass]['inv_probabilities']=pair_dict[key]['inv_probability']
            mz_dict[mass]['transition_energies']=pair_dict[key]['transition_energy']
    mz_dict[mass]['rate_constant']=np.unique(mz_dict[mass]['rate_constants'])
    mz_dict[mass]['inverse_k']=np.unique(mz_dict[mass]['inverse_ks'])
    mz_dict[mass]['probability']=np.unique(mz_dict[mass]['probabilities'])
    mz_dict[mass]['inv_probability']=np.unique(mz_dict[mass]['inv_probabilities'])
    mz_dict[mass]['transition_energy']=np.unique(mz_dict[mass]['transition_energies'])
    mz_dict[mass]['sum_inv_probs']=np.sum(mz_dict[mass]['inv_probabilities'])


print_nested_dict(mz_dict)
print_nested_dict(pair_dict)


#Make a pillar diagram ###Not so pretty right now.
drug_obs=pd.read_csv('verapamil_exp.csv')
drug_obs_x=drug_obs['M/Z_fragment']
drug_obs_y=drug_obs['CE']
drug_obs_y=1/drug_obs_y*0.3
k_obs=np.exp(-drug_obs_y/(8.314*.298))
k_obs_inv=-np.log(k_obs)
k_obs_inv
drug_obs_x
drug_obs_y

x_mz=[]
y_int=[]
for mass in mz_peaks:
    x_mz.append(mass)
    y_int.append(mz_dict[mass]['sum_inv_probs'])

x_mz
drug_obs_x
y_int
drug_obs_y
k_obs_inv

#X-axis
# Clear the current figure
plt.clf()
plt.bar(x_mz, y_int, color='blue', alpha=0.7, label='Sample Data')
plt.bar(drug_obs_x, k_obs_inv, color='red', alpha=0.7, label='Drug Observations')
plt.xlabel('m/z')
plt.ylabel('Intensity')
plt.title('m/Z values for Verapamil')
plt.xlim(0, 450) #Adjust for other molecules
#save as png
plt.savefig('mz_int_k_inv.png')



