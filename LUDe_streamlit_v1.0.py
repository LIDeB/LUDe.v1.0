# -*- coding: utf-8 -*-
"""
Created on Mon May 10 14:27:11 2021

@author: Lucas
"""

##### LUDe #####

# Needed packages
from pathlib import Path
import streamlit as st
import pandas as pd
import base64
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs # para calcular Tanimoto similarity
from rdkit.Chem.Scaffolds import MurckoScaffold # para poder calcular los frameworks
from rdkit.Chem import Descriptors, Lipinski # para poder calcular descriptores
from rdkit.Chem import MCS # para calcular la MCS
from sklearn.metrics import  roc_curve
import random
from time import process_time
import os
import time
import sys
from io import StringIO
import numpy as np
    
from openbabel import openbabel
from molvs import Standardizer

start = time.time()

#---------------------------------#
# Page layout
## Page expands to full width
st.set_page_config(page_title='LIDEB Tools - Decoys',
    layout='wide')

######
# Function to put a picture as header   
def img_to_bytes(img_path):
    img_bytes = Path(img_path).read_bytes()
    encoded = base64.b64encode(img_bytes).decode()
    return encoded

from PIL import Image
image = Image.open('cropped-header-lude.png')
st.image(image)


st.write("[![Website](https://img.shields.io/badge/website-LIDeB-blue)](https://lideb.biol.unlp.edu.ar) [![Twitter Follow](https://img.shields.io/twitter/follow/LIDeB_UNLP?style=social)](https://twitter.com/intent/follow?screen_name=LIDeB_UNLP)")
st.subheader(":pushpin: About Us")
st.markdown("We are a drug discovery team with an interest in the development of publicly available open-source customizable cheminformatics tools to be used in computer-assisted drug discovery. We belong to the Laboratory of Bioactive Research and Development (LIDeB) of the National University of La Plata (UNLP), Argentina. Our research group is focused on computer-guided drug repurposing and rational discovery of new drug candidates to treat epilepsy and neglected tropical diseases.")


# Start the stopwatch / counter 
t1_start = process_time() 


#---------------------------------#
st.write("""
# LIDeB Tools - LUDe

LUDe (LIDEB’s Useful Decoys) is a WebApp that generates, from a set of active compounds, decoys (putative inactive compounds)
which can be used to retrospectively validate virtual screening tools/protocols. Decoys are molecules that have not been tested
against a molecular target of interest but due to their structural features are presumably not prone to bind the target with
high affinity. LUDe finds decoys in a curated ChEMBL27 database; these decoys are paired with the known active compounds in
relation to certain general physicochemical properties (e.g., molecular weight, log P, and others) but are topologically different
from the query compounds. LUDe is conceptually similar to the Directory of Useful Decoys enhanced, but additional filters have been
serially implemented to assure the topological dissimilarity between the decoys and the active compounds used as input.

In this Web App, decoys are obtained through four sequential steps:
    

- Searching molecules with similar physicochemical properties of the input active molecules in a curated ChEMBL database.
- Filtering the selected molecules by dissimilarity against each individual input molecule.
- Randomly selecting a desired number of decoys for each individual input molecule.
- Filtering the selected molecules by the dissimilarity against all the input molecules (set of active compounds used as query). The decoys will have low Tanimoto similarity with the input compounds, and also low Maximum Common Substructure (MCS) ratio and distinctive molecular frameworks (Murcko scaffold). All in all, these three serial filters assure that the decoys will have a different topology in relation to the active compound used as a query. Furthermore, dissimilarity to the remaining active compounds is also checked. Finally, you can download a file with your decoys.

The next workflow summarizes the steps performed by this method:

""")

text = '''
---

'''


image = Image.open('workflow_lude.png')
st.image(image, caption='LUDe Workflow')

st.markdown(text)

st.subheader(":rocket:" "**Fast Tutorial** " "[LUDe](https://www.youtube.com/watch?v=brsC0CPS9U0&ab_channel=LIDeBUNLP)")
st.markdown(" ")


# Name of the dataset
smiles_dataset = "example_molecules"

directory_chembl = str(Path("Base_ChEMBL_final"))


# OPTIONS

#---------------------------------#
# Sidebar - Collects user input features into dataframe
st.sidebar.header('Upload your SMILES')

uploaded_file_1 = st.sidebar.file_uploader("Upload a TXT file with one SMILES per line", type=["txt"])

st.sidebar.markdown("""
[Example TXT input file](https://raw.githubusercontent.com/Capigol/LUDe_v1/main/example_molecules.txt)
""")


lude_setting = st.sidebar.checkbox('Check to change the default configuration', value=False)
if lude_setting == True:
    # Sidebar - Physicochemical features
    st.sidebar.header('Physicochemical features limits')
    lim_MW = st.sidebar.slider('Molecular weight', 0, 50, 20, 1)
    lim_logP = st.sidebar.slider('LogP', 0.0, 1.0, 0.5, 0.1)
    lim_rb = st.sidebar.slider('Rotable bonds', 0, 3, 1, 1)
    lim_Hba = st.sidebar.slider('Num of H Acceptors', 0, 3, 1, 1)
    lim_Hbd = st.sidebar.slider('Num of H Donors', 0, 3, 1, 1)
    lim_charge = st.sidebar.slider('Charge', 0, 3, 1, 1)
    
    # Sidebar - Dissimilarty conditions
    st.sidebar.subheader('Dissimilarty conditions')
    
    fingerprint_radio = st.sidebar.slider('Set the fingerprint radio', 1, 3, 2, 1)
    fingerprint_lenght = st.sidebar.slider('Set the fingerprint lenght', 512, 2048, 1024, 512)
    similarity_metric = st.sidebar.selectbox("Select the similarity metric", ("TanimotoSimilarity", "DiceSimilarity", "CosineSimilarity", "SokalSimilarity", "RusselSimilarity", "KulczynskiSimilarity", "McConnaugheySimilarity"),0)
    
    max_similarity_limit = st.sidebar.slider('Maximum allowed similarity value', 0.0, 1.0, 0.2, 0.05)
    
    # Limit of the fraction of the Maximum Common Substructure
    lim_fraction_MCS = st.sidebar.slider('Limit of the fraction of the Maximum Common Substructure', 0.0, 1.0, 0.5, 0.05)
    
    # Decoys with different framework. Options "Yes" or "No"
    framework_filter = st.sidebar.checkbox('Framework filter',value=True)
    
    # Maximum similarity allowed between decoys and any of the actives
    max_similarity_limit_all = st.sidebar.slider('Maximum similarity allowed between decoys and any of the actives', 0.0, 0.7, 0.2, 0.05)
    
    # Sidebar - Max number of decoys
    st.sidebar.subheader('Maximum number of decoys')
    
    # Max number of decoys by active
    max_decoys = st.sidebar.slider('Maximum number of decoys by active compound', 50, 1000, 50, 50)

else:
    lim_MW = 5
    lim_logP = 0.5
    lim_rb = 1
    lim_Hba = 1
    lim_Hbd = 1
    lim_charge = 1
    fingerprint_radio = 2
    fingerprint_lenght = 1024
    similarity_metric ="TanimotoSimilarity"
    max_similarity_limit = 0.2
    lim_fraction_MCS = 0.5
    framework_filter = True
    max_similarity_limit_all = 0.2
    max_decoys = 50

    
    
st.sidebar.title(":speech_balloon: Contact Us")
st.sidebar.info(
"""
If you are looking to contact us, please
[:e-mail:](mailto:lideb@biol.unlp.edu.ar) or [Twitter](https://twitter.com/LIDeB_UNLP)
""")

#---------------------------------#
# Main panel

# Displays the dataset


st.markdown("""
         ** :point_left: On the left you can find the parameters to set**
         """)


#%%

def decoy_fase1(loaded_smiles, verbose = False):
        
    #import gc
    my_molecules = loaded_smiles[0].tolist()
    
    my_bar = st.progress(0)
    tamanio_total = len(my_molecules)
    t = st.empty()

    active_standarized_smiles = []
    smiles_seleccionados = []
    fp_query = []
    df_analysis = pd.DataFrame()   
    s = Standardizer()
    df_decoy_complete = pd.DataFrame()
    
    for i, molecules in enumerate(my_molecules,start = 1):
       
        t.markdown("Progress: " + str(i) +"/" + str(tamanio_total))
        my_bar.progress(i + 1)

        # Update progress bar
        # if verbose:
        #     print(f"\rProgress: {str(i)} / {str(len(my_molecules))}" )
        
        nombre = "Query_" + str(i)
        
        filtro_physicochemical = [] # molecules passing physicochemical filters
        filtro_tanimoto=[]          # molecules passing similarity filter        
        filtro_MCS=[]               # molecules passing MCS filter
        filtro_fw=[]                # molecules passing framework filter
        decoy_by_active = []

        molecula_ok = molecules.strip()
        
        # Standardization
        mol_standarizado = Standardization(molecula_ok, i, s)
        active_standarized_smiles.append(Chem.MolToSmiles(mol_standarizado))
        # Framework
        core = MurckoScaffold.GetScaffoldForMol(mol_standarizado)
        framework_query = Chem.MolToSmiles(core)
    
        # DESCRIPTORES
        molec_descrip, size_molec = cal_descriptors_DB(mol_standarizado)
        
        # FINGERPRINT
        fp_1 = AllChem.GetMorganFingerprintAsBitVect(mol_standarizado,fingerprint_radio,nBits = fingerprint_lenght,useFeatures=False)
        fp_query.append(fp_1)
        
        # BASE DE DATOS
        # random seed a partir del i de Query, por eso cambia la primera database para cada Query
        # siempre la primer Query va a tener seed = 1, por lo que seleciona df_decoys_physicochem mismo orden
        random.seed(i) 
        total_databases = os.listdir(directory_chembl)
        while total_databases:  
            round_physicochemical = 1
            df_decoys_physicochem = pd.DataFrame()
            # tomo una base de datos random del total, asi no empiezo siempre con la primera
            index = random.randrange(len(total_databases))
            database_random = total_databases[index]
            del total_databases[index]
            df_database = pd.read_csv(directory_chembl + '/' + database_random, sep="\t",index_col=False, header='infer') #abro la base de datos para comparar
            
            # instead of workig with the descriptors value, we worked with the difference 
            # between the decoys molecules and the Query molecule 
            df_database['MW'] = abs(df_database['MW'] - molec_descrip['MW'])
            df_database['LogP'] = abs(df_database['LogP'] - molec_descrip['LogP'])
            df_database['rb'] = abs(df_database['rb'] - molec_descrip['rb'])
            df_database['Hba'] = abs(df_database['Hba'] - molec_descrip['Hba'])
            df_database['Hbd'] = abs(df_database['Hbd'] - molec_descrip['Hbd'])
            df_database['charge'] = abs(df_database['charge'] - molec_descrip['charge'])
            
            lim_MW_DB, lim_logP_DB, lim_rb_DB, lim_Hba_DB, lim_Hbd_DB, lim_charge_DB = lim_MW, lim_logP, lim_rb, lim_Hba, lim_Hbd, lim_charge
            
            while len(df_decoys_physicochem) < 400 or round_physicochemical == 5:
                round_physicochemical = round_physicochemical + 1
                # Applying Physicochemical features limits
                df_decoys_physicochem = df_database[
                    (df_database['MW'] <= lim_MW_DB) &
                    (df_database['LogP'] <= lim_logP_DB) &
                    (df_database['rb'] <= lim_rb_DB) &
                    (df_database['Hba']<= lim_Hba_DB) &
                    (df_database['Hbd']<= lim_Hbd_DB) &
                    (df_database['charge']<= lim_charge_DB)].copy()

                # if there are more than 400 molecules than pass the Physicochemical filter 
                # if there are less than 400 molecules the Physicochemical limits are widen
                if len(df_decoys_physicochem) < 400:
                    # print('The limites were extended')
                    lim_MW_DB = lim_MW_DB * 1.5 
                    lim_logP_DB = lim_logP_DB * 1.5 
                    lim_rb_DB = lim_rb_DB * 1.5 
                    lim_Hba_DB = lim_Hba_DB * 1.5 
                    lim_Hbd_DB = lim_Hbd_DB * 1.5 
                    lim_charge_DB = lim_charge_DB * 1.5
                else:
                    # Calculates the PSS Score for the decoys 
                    df_decoys_PSS = Calculate_PSS(df_decoys_physicochem, molec_descrip) 
                    # Only keeps the top 200 compounds with higher PSS
                    df_decoys_PSS.sort_values(by = 'PSS', ascending = False, inplace = True)
                    df_decoys_PSS = df_decoys_PSS.head(200)
                    break

            filtro_physicochemical.append(len(df_decoys_PSS))
            
            # Calculates MolFromSmiles from selected decoys with Physicochemical
            df_decoys_PSS['mol'] = df_decoys_PSS['SMILE'].apply(lambda x: Chem.MolFromSmiles(x))  ['mol'] = df_decoys_physicochem['SMILE'].apply(lambda x: Chem.MolFromSmiles(x))  
            
            # Dissimilarity filter for each decoy with molecule "Query_i"
            filtro_fw, filtro_MCS, filtro_tanimoto = Dissimilarity_filter(df_decoys_PSS, fp_1, 
                                            mol_standarizado, framework_query, size_molec, smiles_seleccionados, 
                                            decoy_by_active, filtro_fw, filtro_MCS, filtro_tanimoto, nombre, database_random)
            
            # st.write(database_random)
            # st.write('Decoys per active that pass physico and similarity', len(decoy_by_active))
            # instead of having too many decoys we only keep 500 
            if len(decoy_by_active) > 500:
                break

        # after iteration of all the ChEMBL databases or if after the filster 300 decoys per active are reached
        serie_OK = pd.Series({"Query": nombre, "Selected by physicochemical properties": sum(filtro_physicochemical),
                              "Pass the Tc filter per active": len(filtro_tanimoto),"Pass the fMCS filter": len(filtro_MCS), 
                              "Pass the Framework filter": len(filtro_fw), 'Non duplicated Decoys per active:' : len(decoy_by_active)})

        df_decoy_complete = pd.concat([df_decoy_complete, pd.DataFrame(decoy_by_active)], axis = 0, ignore_index= True)
        
        df_analysis = pd.concat([df_analysis, serie_OK], axis = 1, ignore_index=True)
    
    df_analysis = df_analysis.transpose()
    df_analysis.set_index('Query', drop=True, inplace = True)
    
    st.write('')
    st.write('Decoys have been successfully obtained for each loaded molecule!!')
    st.write("---------------------------------------------------------------")
    st.write("Molecules that passed the physicochemical properties filters: " + str(df_analysis['Selected by physicochemical properties'].sum()))
    st.write("Molecules that passed the dissimilarity (structural) filters: " + str(len(smiles_seleccionados)))
    st.write("Of which: " + str(df_decoy_complete.shape[0]) + " are different")

    return (df_decoy_complete, fp_query, df_analysis, active_standarized_smiles)


#%%
def Standardization(molecula_ok, i, s):

    try:
        mol = Chem.MolFromSmiles(molecula_ok) # convierte las moleculas a mol
        mol_standarizado = s.fragment_parent(mol) #Return the fragment parent of a given molecule, the largest organic covalent unit in the molecule
        # mol_standarizado = s.stereo_parent(mol_standarizado, skip_standardize= True) #Return The stereo parentof a given molecule, has all stereochemistry information removed from tetrahedral centers and double bonds.
        mol_standarizado = s.charge_parent(mol_standarizado, skip_standardize= True) #Return the charge parent of a given molecule,  the uncharged version of the fragment parent
        mol_standarizado = s.isotope_parent(mol_standarizado, skip_standardize= True) #Return the isotope parent of a given molecule, has all atoms replaced with the most abundant isotope for that element.
       
        smile_standarizado = Chem.MolToSmiles(mol_standarizado)
        ionized_smile = charges_ph(smile_standarizado)
        smile_checked = smile_obabel_corrector(ionized_smile)
        
        mol_checked = Chem.MolFromSmiles(smile_checked)
    except:
        st.write("**Oh no! There is a problem with standarization of one SMILES.**")
        st.write("**Please check your molecule: **" + str(i))
        st.write("**That is the SMILES: **" + str(molecula_ok))
        sys.exit()
    
    return mol_checked

#%%
def charges_ph(molecule):
    
    # obConversion it's neccesary for saving the objects
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("smi", "smi")
    
    # create the OBMol object and read the SMILE
    mol = openbabel.OBMol()
    obConversion.ReadString(mol, molecule)
    
    # Add H, correct pH and add H again, it's the only way it works
    mol.AddHydrogens()
    mol.CorrectForPH(7.4)
    mol.AddHydrogens()
    
    # transforms the OBMOl objecto to string (SMILES)
    optimized = obConversion.WriteString(mol)
    
    return optimized

#%%

def smile_obabel_corrector(smiles_ionized):
    mol1 = Chem.MolFromSmiles(smiles_ionized, sanitize = False)
    
    # checks if the ether group is wrongly protonated
    pattern1 = Chem.MolFromSmarts('[#6]-[#8-]-[#6]')
    if mol1.HasSubstructMatch(pattern1):
        # gets the atom number for the O wrongly charged
        at_matches = mol1.GetSubstructMatches(pattern1)
        at_matches_list = [y[1] for y in at_matches]
        # changes the charged for each O atom
        for at_idx in at_matches_list:
            atom = mol1.GetAtomWithIdx(at_idx)
            atom.SetFormalCharge(0)
            atom.UpdatePropertyCache()

    pattern12 = Chem.MolFromSmarts('[#6]-[#8-]-[#16]')
    if mol1.HasSubstructMatch(pattern12):
        # gets the atom number for the O wrongly charged
        at_matches = mol1.GetSubstructMatches(pattern12)
        at_matches_list = [y[1] for y in at_matches]
        # changes the charged for each O atom
        for at_idx in at_matches_list:
            atom = mol1.GetAtomWithIdx(at_idx)
            atom.SetFormalCharge(0)
            atom.UpdatePropertyCache()
            
    # checks if the nitro group is wrongly protonated in the oxygen
    pattern2 = Chem.MolFromSmarts('[#6][O-]=[N+](=O)[O-]')
    if mol1.HasSubstructMatch(pattern2):
        # print('NO 20')
        patt = Chem.MolFromSmiles('[O-]=[N+](=O)[O-]', sanitize = False)
        repl = Chem.MolFromSmiles('O[N+]([O-])=O')
        rms = AllChem.ReplaceSubstructs(mol1,patt,repl,replaceAll=True)
        mol1 = rms[0]

    # checks if the nitro group is wrongly protonated in the oxygen
    pattern21 = Chem.MolFromSmarts('[#6]-[O-][N+](=O)=[O-]')
    if mol1.HasSubstructMatch(pattern21):
        # print('NO 21')
        patt = Chem.MolFromSmiles('[O-][N+](=O)=[O-]', sanitize = False)
        repl = Chem.MolFromSmiles('[O][N+](=O)-[O-]')
        rms = AllChem.ReplaceSubstructs(mol1,patt,repl,replaceAll=True)
        mol1 = rms[0]
        
    # checks if the nitro group is wrongly protonated, different disposition of atoms
    pattern22 = Chem.MolFromSmarts('[#8-][N+](=[#6])=[O-]')
    if mol1.HasSubstructMatch(pattern22):
        # print('NO 22')
        patt = Chem.MolFromSmiles('[N+]([O-])=[O-]', sanitize = False)
        repl = Chem.MolFromSmiles('[N+]([O-])-[O-]')
        rms = AllChem.ReplaceSubstructs(mol1,patt,repl,replaceAll=True)
        mol1 = rms[0]

    # checks if the nitro group is wrongly protonated, different disposition of atoms
    pattern23 = Chem.MolFromSmarts('[#6][N+]([#6])([#8-])=[O-]')
    if mol1.HasSubstructMatch(pattern23):
        # print('NO 23')
        patt = Chem.MolFromSmiles('[N+]([O-])=[O-]', sanitize = False)
        repl = Chem.MolFromSmiles('[N+]([O-])[O-]')
        rms = AllChem.ReplaceSubstructs(mol1,patt,repl,replaceAll=True)
        mol1 = rms[0]

    # checks if the nitro group is wrongly protonated, different disposition of atoms
    pattern24 = Chem.MolFromSmarts('[#6]-[#8][N+](=O)=[O-]')
    if mol1.HasSubstructMatch(pattern24):
        # print('NO 24')
        patt = Chem.MolFromSmiles('[O][N+](=O)=[O-]', sanitize = False)
        repl = Chem.MolFromSmiles('[O][N+](=O)[O-]')
        rms = AllChem.ReplaceSubstructs(mol1,patt,repl,replaceAll=True)
        mol1 = rms[0]

    # checks if the 1H-tetrazole group is wrongly protonated
    pattern3 = Chem.MolFromSmarts('[#7]-1-[#6]=[#7-]-[#7]=[#7]-1')
    if mol1.HasSubstructMatch(pattern3):
        # gets the atom number for the N wrongly charged
        at_matches = mol1.GetSubstructMatches(pattern3)
        at_matches_list = [y[2] for y in at_matches]
        # changes the charged for each N atom
        for at_idx in at_matches_list:
            atom = mol1.GetAtomWithIdx(at_idx)
            atom.SetFormalCharge(0)
            atom.UpdatePropertyCache()

    # checks if the 2H-tetrazole group is wrongly protonated
    pattern4 = Chem.MolFromSmarts('[#7]-1-[#7]=[#6]-[#7-]=[#7]-1')
    if mol1.HasSubstructMatch(pattern4):
        # gets the atom number for the N wrongly charged
        at_matches = mol1.GetSubstructMatches(pattern4)
        at_matches_list = [y[3] for y in at_matches]
        # changes the charged for each N atom
        for at_idx in at_matches_list:
            atom = mol1.GetAtomWithIdx(at_idx)
            atom.SetFormalCharge(0)
            atom.UpdatePropertyCache()
        
    # checks if the 2H-tetrazole group is wrongly protonated, different disposition of atoms
    pattern5 = Chem.MolFromSmarts('[#7]-1-[#7]=[#7]-[#6]=[#7-]-1')
    if mol1.HasSubstructMatch(pattern5):
        # gets the atom number for the N wrongly charged
        at_matches = mol1.GetSubstructMatches(pattern4)
        at_matches_list = [y[4] for y in at_matches]
        # changes the charged for each N atom
        for at_idx in at_matches_list:
            atom = mol1.GetAtomWithIdx(at_idx)
            atom.SetFormalCharge(0)
            atom.UpdatePropertyCache()
    smile_checked = Chem.MolToSmiles(mol1)
   
    return smile_checked


#%%

def cal_descriptors_DB(mol_standarizado):

    # Calculate properties and store in dict
    prop_dict = {}
    # molweight
    prop_dict.update({'MW': Descriptors.MolWt(mol_standarizado)})
    # logP
    prop_dict.update({'LogP': Chem.Crippen.MolLogP(mol_standarizado)})
    # HBA
    prop_dict.update({'Hba': Chem.rdMolDescriptors.CalcNumLipinskiHBA(mol_standarizado)})
    # HBD
    prop_dict.update({'Hbd': Chem.rdMolDescriptors.CalcNumLipinskiHBD(mol_standarizado)})
    # rotatable bonds
    prop_dict.update({'rb': Chem.rdMolDescriptors.CalcNumRotatableBonds(mol_standarizado)})
    # Formal (net) charge
    prop_dict.update({'charge': Chem.rdmolops.GetFormalCharge(mol_standarizado)})

    size_molec = Descriptors.HeavyAtomCount(mol_standarizado)
    
    # molec_descrip = pd.Series({'MW' : MolWt, 'LogP': MolLogP, 'Num_Rotatable_Bonds': NumRotatableBonds,'Num_H_Acceptors': NumHAcceptors, 'Num_H_Donors': NumHDonors})

    return prop_dict, size_molec

#%%

def Calculate_PSS(df_decoys_physicochem, molec_descrip):
    ''' 
    Calculates the physicochemical similarity score (PSS)
    for each decoy according to their query molecule
    '''
    # for each Active row with their physicochem properties
    df_final = pd.DataFrame()
    relative_distance= []
    
    df_only_property = df_decoys_physicochem.copy()    
    df_only_property.drop(labels = ['SMILE', 'Framework'], axis = 1, inplace = True)
    
    df_only_framework = df_decoys_physicochem.copy()    
    df_only_framework.drop(labels = ['MW', 'LogP', 'rb', 'Hba', 'Hbd', 'charge'], axis = 1, inplace = True)
    
    # Calculates the max distance between the query_i and all the decoys for the Physicochem prop
    distances = Physicochem_distances(df_only_property, molec_descrip)
    
    # Iterate for each property and distance of the active
    for prop, value in molec_descrip.items():
        relative_distance_prop = []
        f_prop = distances[f'f_{prop}']
        decoy_prop_physicochem = df_only_property[prop]
        for x in decoy_prop_physicochem:
            if f_prop == 0:
                relative_distance_prop.append(1)
            else:
                normalization = 1 - (abs(value - x) / f_prop)
                relative_distance_prop.append(normalization)
        relative_distance.append(relative_distance_prop)
        
    df_dist = pd.DataFrame(relative_distance)
    df_final = pd.concat([df_final, df_dist], axis = 1, ignore_index= True)
   
    df_final = df_final.transpose()
    df_final = df_final.set_index(df_only_framework.index)
    df_only_framework['PSS'] = df_final.mean(axis = 1)
    
    return df_only_framework

#%%
def Physicochem_distances(df_only_property, molec_descrip):
    '''
    Calculates the maximum distances for the five Physicochem descriptors
    between the query molecule and all the selected decoy from the Physicochem df
    '''
    distances = []
    
    physicochem_min = pd.Series([df_only_property[column].min() for column in df_only_property],
                                index = ['MW', 'LogP', 'rb', 'Hba', 'Hbd', 'charge'])
    physicochem_max = pd.Series([df_only_property[column].max() for column in df_only_property],  
                                index = ['MW', 'LogP', 'rb', 'Hba', 'Hbd', 'charge'])
    for idx, value in physicochem_min.items():
        if abs(molec_descrip[idx] - physicochem_min[idx]) < abs(molec_descrip[idx] - physicochem_max[idx]):
            distances.append(abs(molec_descrip[idx] - physicochem_max[idx]))
        else:
            distances.append(abs(molec_descrip[idx] - physicochem_min[idx]))
    
    distances_serie = pd.Series(distances)
    distances_serie.rename(index = {0:'f_MW' , 1:'f_LogP', 2:'f_rb',3:'f_Hba', 4:'f_Hbd', 5:'f_charge'}, inplace=True)
    
    return distances_serie


#%%

def Dissimilarity_filter(df_decoys_PSS, fp_1, mol_standarizado, framework_query, size_molec, smiles_seleccionados, decoy_by_active, filtro_fw, filtro_MCS, filtro_tanimoto, nombre, database_random):
    ''' 
    Dissimilarity conditions for each decoy with molecule "Query_i", which involves:
    filter by Maximum allowed similarity value, calculating Tanimoto similarity
    filter by limit of the fraction of the Maximum Common Substructure
    filter by framework if the framework of the decoys are the same as the Query
    '''
    # st.write(f'Calculating Maximum allowed similarity and MCS for {database_random}')
    
    # filter by Maximum allowed similarity value, calculating Tanimoto similarity
            
    for index, decoy_in_df in df_decoys_PSS.iterrows():
    
        smiles_DB = decoy_in_df['SMILE']
        fp_2 = AllChem.GetMorganFingerprintAsBitVect(decoy_in_df["mol"],fingerprint_radio, 
                                                     nBits = fingerprint_lenght, useFeatures=False)
        similarity_metric_ok = getattr(DataStructs, similarity_metric)
        tan_sim = similarity_metric_ok(fp_1, fp_2)
        
        # Tanimoto similarity Maximum allowed similarity value
        if tan_sim <= float(max_similarity_limit):
            mols = [mol_standarizado,decoy_in_df["mol"]]
            filtro_tanimoto.append(smiles_DB)
            res = MCS.FindMCS(mols, timeout=0.1)
            tamanio_MCS = res.numAtoms
            
            # Limit of the fraction of the Maximum Common Substructure
            if tamanio_MCS / size_molec < lim_fraction_MCS:
                filtro_MCS.append(smiles_DB)
                if framework_filter == True:
                    framework_DB = decoy_in_df['Framework']
                    if framework_query != framework_DB:
                        if not smiles_DB in smiles_seleccionados:
                            decoy_by_active.append(pd.Series({'SMILE': smiles_DB, 'fp':fp_2, 'Query': nombre}))
                        filtro_fw.append(smiles_DB)
                        smiles_seleccionados.append(smiles_DB)
                    else:
                        pass
                else:
                    if not smiles_DB in smiles_seleccionados:
                        decoy_by_active.append(pd.Series({'SMILE': smiles_DB, 'fp':fp_2, 'Query': nombre}))
                    filtro_fw.append(smiles_DB)
                    smiles_seleccionados.append(smiles_DB)
                 
    
     
    return filtro_fw, filtro_MCS, filtro_tanimoto


#%% To take XXX random decoys from each list

def duplicates_filter(df_decoy_complete, fp_query, df_analysis):
    final_decoys = pd.DataFrame()

    # Comparing actives and decoys by tanimoto
    df_Tanimoto_final = Tanimoto_actives_decoys(df_decoy_complete, fp_query, max_similarity_limit_all)
    number_decoys_TC = df_Tanimoto_final.value_counts(subset = ['Query'], sort = False)
    number_decoys_TC.rename('Pass the TC filter total actives:', inplace = True)
    df_analysis_final = pd.merge(df_analysis, number_decoys_TC,how = 'left', left_on = 'Query',right_on = 'Query')
    
    # keeping only max_decoys per Query_j, random selection
    for j, fp_ in enumerate(fp_query, start = 1):
        decoys_from_query = df_Tanimoto_final[df_Tanimoto_final['Query'] == f'Query_{j}'] 
        if decoys_from_query.shape[0] >= max_decoys:
            df_decoys_max = decoys_from_query.sample(n = max_decoys)
        else:
            df_decoys_max = decoys_from_query

        final_decoys = pd.concat([final_decoys, df_decoys_max], axis = 0, ignore_index= True)
        
    st.write("Finally, " + str(df_Tanimoto_final.shape[0]) + " passed the Tanimoto filter by comparing all loaded molecules")
    st.write('')
    st.write('Congratulations, you have obtained ' + str(final_decoys.shape[0]) + " decoys!!!")
    st.write("---------------------------------------------------------------")
    final_decoys.drop(labels = ['fp'], axis = 1, inplace = True)
    return final_decoys, df_analysis_final
   
#%%

def Tanimoto_actives_decoys(df_decoy_complete, fp_query, max_similarity_limit_all):
    ''' 
    From the df with decoys that pass the filters,
    only the rows that are dissimilar to all the actives are kept
    '''
    Tanimoto_all_decoy = []
    
    for _, decoy_in_df in df_decoy_complete.iterrows():
        coef_tan= []
        for fp_activo in fp_query:
            tan_sim=DataStructs.TanimotoSimilarity(decoy_in_df['fp'], fp_activo)
            coef_tan.append(tan_sim)
        if max(coef_tan) < max_similarity_limit_all:
            Tanimoto_all_decoy.append(decoy_in_df)
    df_Tanimoto_final = pd.DataFrame(Tanimoto_all_decoy)    
    
    return df_Tanimoto_final

#%%

def Doppelganger_score_calculation(decoys_fp, actives_fp):
    ''' 
    Following Vogel_2011 DEKOIS
    From the decoys Series with their FP compares with the actives Series with FP
    and calculates the similarity matrix, selecting the max value 
    
    En DeepCoy dice: 
    For each decoy molecule, its doppelganger score is the maximum similarity across all actives
    '''
    Tanimoto_all_decoy = []
    # Creates the similarity matrix for Query and Decoys
    for fp_decoy in decoys_fp:
        coef_tan = []
        for fp_activo in actives_fp:
            tan_sim=DataStructs.TanimotoSimilarity(fp_decoy, fp_activo)
            coef_tan.append(tan_sim)
        Tanimoto_all_decoy.append(coef_tan) 
    df_Tanimoto_final = pd.DataFrame(Tanimoto_all_decoy)
    
    # Determine the maximum TC for each Query atom to all Decoys
    max_TC_decoys = df_Tanimoto_final.max(axis = 0)
    
    max_doppelganger_score = max_TC_decoys.max()
    doppelganger_score = max_TC_decoys.mean()
    
    return df_Tanimoto_final, doppelganger_score, max_doppelganger_score

#%%

def cal_descriptors(mol_standarizado):
    ''' 
    Calculates the same 6 molecular descriptors as in Dud-E
    las query de Alexander calculaba HBD y HBA de dos formas con valores diferentes :
        Chem.Lipinski.NumHAcceptors()
        rdMolDescriptors.CalcNumLipinskiHBA()
    Según mails con Alexander, Lipinski.NumHAcceptors da mejores resultados
    para generar los LUDe empleamos esa funcion y por lo tanto aca tmb
    '''
    
    MolWt = Descriptors.MolWt(mol_standarizado) 
    MolLogP =  Chem.Crippen.MolLogP(mol_standarizado)
    NumRotatableBonds = Lipinski.NumRotatableBonds(mol_standarizado) 
    NumHAcceptors =  Lipinski.NumHAcceptors(mol_standarizado) 
    NumHDonors =  Lipinski.NumHDonors(mol_standarizado) 
    Formalcharge = Chem.rdmolops.GetFormalCharge(mol_standarizado)
    
    molec_descrip = [MolWt, MolLogP, NumRotatableBonds, NumHAcceptors,
                     NumHDonors, Formalcharge]
    return molec_descrip


def dataset_preparation(dataset_smiles, name: str):
    ''' 
    Given the Query or the decoys dataset
    standardize the molecules, calculates the FP and the molecular descriptors
    '''
    dataset_mol = [Chem.MolFromSmiles(smile) for smile in dataset_smiles]       
    dataset_desc_np = np.array([cal_descriptors(mol_active) for mol_active in dataset_mol])
    
    return dataset_desc_np


def normalize_descriptors_percentile(all_desc):
    ''' 
    Following Vogel_2011 DEKOIS
    for each descriptor, calculates the 95th and the 5th percentile of all the molecules
    then do the rest between the 95th and the 5th percentiles and
    normalizes the descriptors values with the aforementioned difference
    '''
    normalized_all_desc = pd.DataFrame()
    percentile = all_desc.quantile(0.95) - all_desc.quantile(0.05)
        
    for column in all_desc:
        normalized_all_desc[column] = all_desc[column]/percentile[column]
        
    return normalized_all_desc


#%%

def doe_score(actives, decoys):
    ''' 
    Following Imrie_2021 DeepCoy
    Calculates the DOE score with the script from the paper
    for every active measure the distance with every molecule
    calculates the AUC ROC 
    '''
    all_feat = list(actives) + list(decoys)
    up_p = np.percentile(all_feat, 95, axis=0)
    low_p = np.percentile(all_feat, 5, axis=0)
    norms = up_p - low_p
    for i in range(len(norms)):
        if norms[i] == 0:
            norms[i] = 1.

    active_norm = [act/norms for act in actives]
    decoy_norm = [dec/norms for dec in decoys]
    all_norm = active_norm + decoy_norm

    active_embed = []
    labels = [1] * (len(active_norm)-1) + [0] * len(decoy_norm)
    for i, act in enumerate(active_norm):
        comp = list(all_norm)
        del comp[i]
        dists = [100 - np.linalg.norm(c-act) for c in comp] # arbitrary large number to get scores in reverse order
        fpr, tpr, _ = roc_curve(labels, dists)
        fpr = fpr[::]
        tpr = tpr[::]
        a_score = 0
        for i in range(len(fpr)-1):
            a_score += (abs(0.5*( (tpr[i+1]+tpr[i])*(fpr[i+1]-fpr[i]) - (fpr[i+1]+fpr[i])*(fpr[i+1]-fpr[i]) )))
        active_embed.append(a_score)

    #print(np.average(active_embed))
    return np.average(active_embed)




#%%
# Funcion para exportar el archivo

def filedownload(df):
    csv = df.to_csv(index=True,header=True)
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64,{b64}" download="generated_decoys.csv">Download CSV File with your decoys</a>'
    return href

def filedownload1(df):
    csv = df.to_csv(index=True,header=True)
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64,{b64}" download="decoys_analysis.csv">Download CSV File with the table</a>'
    return href

def filedownload2(df):
    csv = df.to_csv(index=False,header=False)
    b64 = base64.b64encode(csv.encode()).decode()  # strings <-> bytes conversions
    href = f'<a href="data:file/csv;base64,{b64}" download="decoys_settings.csv">Download CSV File with your settings</a>'
    return href



#%% Settings
def setting_info():
    from datetime import date
    today = date.today()
    fecha = today.strftime("%d/%m/%Y")
    settings = []
    settings.append(["Decoys generated at: " , fecha])
    settings.append(["Physicochemical features limits:",""])
    settings.append(["MW" , "+/- " +  str(lim_MW)])
    settings.append(["logP" , "+/- " + str(lim_logP)])
    settings.append(["Num_rotable_bonds" , "+/- " + str(lim_rb)])
    settings.append(["Num_H_acceptors" , "+/- " + str(lim_Hba)])
    settings.append(["Num_H_donors" , "+/- " + str(lim_Hbd)])

    settings.append(["Topological features ---> Dissimilarty conditions",""])
    settings.append(["Morgan Fingerprints",""])
    settings.append(["Fingerprint radio:" , str(fingerprint_radio)])
    settings.append(["Fingerprint lenght:" , str(fingerprint_lenght)])
    if lim_fraction_MCS == 1.0:
        settings.append(["Not set a limit for the fraction of the Maximum Common Substructure",""])
    else:
        settings.append(["Limit of fraction of the Maximum Common Substructure: " , str(lim_fraction_MCS)])
    settings.append(["Decoys with different framework:" , str(framework_filter)])
    settings.append(["Max Tc similarity between decoys and other actives:", str(max_similarity_limit)])
    settings.append(["Max number of decoys by loaded molecule:", str(max_decoys)])
    settings.append(["\n"])
    settings.append(["Decoy validation", " "])
    settings.append(["DOE_score", str(DOE_score_np)])
    settings.append(["Mean Doppelganger_score", str(doppelganger_score)])
    settings.append(["Max Doppelganger_score", str(max_doppelganger_score)])

    settings_df = pd.DataFrame(settings)
    return settings_df




#%%
# ---------------------------------#

if uploaded_file_1 is not None:
    st.subheader(':point_down: Click RUN to generate decoys')
    run = st.button("RUN")
    if run == True:
        time_start = time.time()
        loaded_smiles = pd.read_csv(uploaded_file_1,sep="\t",header=None)
        lista_resultados = decoy_fase1(loaded_smiles)
        st.write("Now, we are comparing all decoys vs all input SMILES, please wait a moment...")
        
        df_final_decoys, df_analysis = duplicates_filter(lista_resultados[0], lista_resultados[1], lista_resultados[2]) 
        
        
        # DOE and Doppelganger calculation
        df_Tanimoto_final, doppelganger_score, max_doppelganger_score = Doppelganger_score_calculation(lista_resultados[0]["fp"], lista_resultados[1])
        actives_desc = dataset_preparation(lista_resultados[3], 'Query')
        decoys_desc = dataset_preparation(lista_resultados[3], 'Decoy')
        DOE_score_np = doe_score(actives_desc, decoys_desc)

        # st.markdown(" **Here you can dowload the generated decoys**", unsafe_allow_html=True)

        st.write("Decoy validation metrics:")
                # df = duplicates_filter(lista_resultados[0])   
        st.write("Decoy validation metrics:")
        st.write(f"{'DOE_score: '}{round(DOE_score_np, 3)}{'    (best possible score = 0, worst possible score = 0.5)'}")
        st.write(f"{'Doppelganger score: '}{round(doppelganger_score, 3)}{'    (best possible score = 0, worst possible score = 1)'}")
        st.write(f"{'Max Doppelganger score: '}{round(max_doppelganger_score, 3)}{'    (best possible score = 0, worst possible score = 1)'}")

        
        st.markdown(":point_down: **Here you can dowload the generated decoys**", unsafe_allow_html=True)
        st.markdown(filedownload(df_final_decoys), unsafe_allow_html=True)
        st.markdown(text)
        st.markdown(":point_down: **Here you can see a little analysis of the process**", unsafe_allow_html=True)
        st.write(df_analysis)
        
        st.markdown(":point_down: **Here you can dowload this table in a csv file**", unsafe_allow_html=True)
        st.markdown(filedownload1(df_analysis), unsafe_allow_html=True)
        st.markdown(text)
        settings_df = setting_info()
        st.markdown(":point_down: **Here you can download your settings**", unsafe_allow_html=True)
        st.markdown(filedownload2(settings_df), unsafe_allow_html=True)


        st.balloons()

else:
    if st.button('Press to use Example SMILES'):
        st.write("Five SMILES have been loaded as example")
        loaded_smiles = pd.read_csv("example_molecules.txt",sep="\t",header=None)
        
        lista_resultados = decoy_fase1(loaded_smiles)
        st.write("Now, we are comparing all decoys vs all input SMILES, please wait a moment...")
        
        df_final_decoys, df_analysis = duplicates_filter(lista_resultados[0], lista_resultados[1], lista_resultados[2]) 
        
        # DOE and Doppelganger calculation
        df_Tanimoto_final, doppelganger_score, max_doppelganger_score = Doppelganger_score_calculation(lista_resultados[0]["fp"], lista_resultados[1])
        actives_desc = dataset_preparation(lista_resultados[3], 'Query')
        decoys_desc = dataset_preparation(lista_resultados[3], 'Decoy')
        DOE_score_np = doe_score(actives_desc, decoys_desc)

        st.write("Decoy validation metrics:")
        st.write(f"{'DOE_score: '}{round(DOE_score_np, 3)}{'    (best possible score = 0, worst possible score = 0.5)'}")
        st.write(f"{'Doppelganger score: '}{round(doppelganger_score, 3)}{'    (best possible score = 0, worst possible score = 1)'}")
        st.write(f"{'Max Doppelganger score: '}{round(max_doppelganger_score, 3)}{'    (best possible score = 0, worst possible score = 1)'}")


        st.markdown(":point_down: **Here you can dowload the generated decoys**", unsafe_allow_html=True)
        st.markdown(filedownload(df_final_decoys), unsafe_allow_html=True)
        st.markdown(text)
        st.markdown(":point_down: **Here you can see a little analysis of the process**", unsafe_allow_html=True)
        st.write(df_analysis)
        
        st.markdown(":point_down: **Here you can dowload this table in a csv file**", unsafe_allow_html=True)
        st.markdown(filedownload1(df_analysis), unsafe_allow_html=True)
        st.markdown(text)
        settings_df = setting_info()
        st.markdown(":point_down: **Here you can download your settings**", unsafe_allow_html=True)
        st.markdown(filedownload2(settings_df), unsafe_allow_html=True)

        st.balloons()
    else:
        st.info('Awaiting for TXT file to be uploaded.')

      
st.write("""
### Cite us:
*LIDeB Tools: A Latin American resource of freely available, open-source cheminformatics apps*

Denis N. Prada Gori, Lucas N. Alberca, Santiago Rodriguez, Juan I.Alice, Manuel A.Llanos, Carolina L. Bellera, Alan Talevi.

Artificial Intelligence in the Life Sciences

[DOI: 10.1016/j.ailsci.2022.100049](https://www.sciencedirect.com/science/article/pii/S2667318522000198)
    
""")        

#Footer edit

footer="""<style>
a:link , a:visited{
color: blue;
background-color: transparent;
text-decoration: underline;
}

a:hover,  a:active {
color: red;
background-color: transparent;
text-decoration: underline;
}

.footer {
position: fixed;
left: 0;
bottom: 0;
width: 100%;
background-color: white;
color: black;
text-align: center;
}

</style>
<div class="footer">
<p>Made in  🐍 and <img style='display: ; ' href="https://streamlit.io" src="https://i.imgur.com/iIOA6kU.png" target="_blank"></img> Developed with ❤️ by <a style='display:; text-align: center;' href="https://lideb.biol.unlp.edu.ar/" target="_blank">LIDeB</a></p>
</div>
"""
st.markdown(footer,unsafe_allow_html=True)



