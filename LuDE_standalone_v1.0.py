# -*- coding: utf-8 -*-
"""
Input: un archivo .TXT con un SMILE por linea de las moleculas que se desean buscar decoys
una carpeta llamada "Databases_ChEMBL_29" con los archivos .TXT donde se van a buscar los decoys, con las columnas:
    SMILES, Framework, MW, LogP, Num Rotatable Bonds, Num H Acceptors, Num H Donor
Ouput: Decoys_analysis un resumen de los decoys generados para cada molecula
Decoys_settings el setting que se uso para generar los decoys
Generated_decoys los decoys generados 

FUNCIONA CON CHEMBL 30 
"""

##### LUDe #####

# Needed packages
from pathlib import Path
import pandas as pd
import os
import sys
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Descriptors, rdFMCS
from rdkit.Chem.Scaffolds import MurckoScaffold 
from molvs import Standardizer
import random
from openbabel import openbabel
import time
start = time.time()

def read_paired_file(filename):
    '''Reads .smi file '''
    '''Returns array containing smiles strings of molecules'''
    smiles, names = [], []

    with open(filename, 'r') as f:
        for line in f:
            if line:
                smiles.append(line.strip().split(' ')[0:2])
    smiles_final = pd.Series([d[0] for d in smiles])
    return smiles_final

# Name of the dataset
smiles_dataset = "SMILES_Activos_hCAA"

# Folder with files
directory = str(Path(r"C:\Users\adm\Desktop\Pruebas"))
# loaded_smiles = pd.read_csv(directory + "\\" + smiles_dataset + '.txt', sep="\t",header=None)
loaded_smiles = read_paired_file(f'{directory}\\{smiles_dataset}.txt')
directory_chembl = str(Path(r'C:\Users\adm\OneDrive-Denis\OneDrive - biol.unlp.edu.ar\Postdoc\Lude\Construction of Database\Bases_Cehmbl_listas'))
# OPTIONS

# Physicochemical features limits
lim_MW = 20         # Molecular Weight
lim_logP = 0.5      # logP
lim_rb = 1   # Rotable bonds
lim_Hba = 1    # Num of H Acceptors
lim_Hbd = 1      # Num of H Donors
lim_charge = 1      # Formal charge

# Dissimilarty conditions
fingerprint_radio =  2          # 1, 2, 3
fingerprint_lenght =  1024      # 512, 1024, 2048
similarity_metric = "TanimotoSimilarity" # "TanimotoSimilarity", "DiceSimilarity", "CosineSimilarity", "SokalSimilarity", "RusselSimilarity", "KulczynskiSimilarity", "McConnaugheySimilarity"
max_similarity_limit = 0.2      # Maximum allowed similarity value: from 0.0 to 1.0

# Limit of the fraction of the Maximum Common Substructure
lim_fraction_MCS = 0.5 # from 0.0 to 1.0

# Decoys with different framework. Options "True" or "False"
framework_filter = True

# Maximum similarity allowed between decoys and any of the actives
max_similarity_limit_all = 0.2 # from 0.0 to 0.7

# Max number of decoys by active compound
max_decoys = 50 

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
        print("**Oh no! There is a problem with standarization of one SMILES.**")
        print("**Please check your molecule: **" + str(i))
        print("**That is the SMILES: **" + str(molecula_ok))
        sys.exit()
    
    return mol_checked

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

def Dissimilarity_filter(df_decoys_PSS, fp_1, mol_standarizado, framework_query, size_molec, smiles_seleccionados, decoy_by_active, filtro_fw, filtro_MCS, filtro_tanimoto, nombre, database_random):
    ''' 
    Dissimilarity conditions for each decoy with molecule "Query_i", which involves:
    filter by Maximum allowed similarity value, calculating Tanimoto similarity
    filter by limit of the fraction of the Maximum Common Substructure
    filter by framework if the framework of the decoys are the same as the Query
    '''
    print(f'Calculating Maximum allowed similarity and MCS for {database_random}')
    
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
            res = rdFMCS.FindMCS(mols)
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


#%%

def decoy_fase1(loaded_smiles, verbose = False):
        
    # my_molecules = loaded_smiles[0].tolist()
    my_molecules = loaded_smiles.tolist()
    
    smiles_seleccionados = []
    fp_query = []
    df_analysis = pd.DataFrame()   
    s = Standardizer()
    df_decoy_complete = pd.DataFrame()
    
    for i, molecules in enumerate(my_molecules,start = 1):
       
        # Update progress bar
        if verbose:
            print(f"\rProgress: {str(i)} / {str(len(my_molecules))}" )
        
        nombre = "Query_" + str(i)
        
        filtro_physicochemical = [] # molecules passing physicochemical filters
        filtro_tanimoto=[]          # molecules passing similarity filter        
        filtro_MCS=[]               # molecules passing MCS filter
        filtro_fw=[]                # molecules passing framework filter
        decoy_by_active = []

        molecula_ok = molecules.strip()
        
        # Standardization
        mol_standarizado = Standardization(molecula_ok, i, s)
        
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
            df_database = pd.read_csv(directory_chembl + '\\' + database_random, sep="\t",index_col=False, header='infer') #abro la base de datos para comparar
            
            # instead of workig with the descriptors value, we worked with the difference 
            # between the decoys molecules and the Query molecule 
            df_database['MW'] = abs(df_database['MW'] - molec_descrip['MW'])
            df_database['LogP'] = abs(df_database['LogP'] - molec_descrip['LogP'])
            df_database['rb'] = abs(df_database['rb'] - molec_descrip['rb'])
            df_database['Hba'] = abs(df_database['Hba'] - molec_descrip['Hba'])
            df_database['Hbd'] = abs(df_database['Hbd'] - molec_descrip['Hbd'])
            df_database['charge'] = abs(df_database['charge'] - molec_descrip['charge'])
            
            lim_MW_DB, lim_logP_DB, lim_rb_DB, lim_Hba_DB, lim_Hbd_DB, lim_charge_DB = lim_MW, lim_logP, lim_rb, lim_Hba, lim_Hbd, lim_charge
            
            while len(df_decoys_physicochem) < 400 or round_physicochemical < 5:
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
                    # print(round_physicochemical, lim_MW_DB, lim_logP_DB, lim_rb_DB, lim_Hba_DB, lim_Hbd_DB, lim_charge_DB)
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
            
            print(database_random)
            print('Decoys per active that pass physico and similarity', len(decoy_by_active))
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
    
    print('')
    print('Decoys have been successfully obtained for each loaded molecule!!')
    print("---------------------------------------------------------------")
    print("Molecules that passed the physicochemical properties filters: " + str(df_analysis['Selected by physicochemical properties'].sum()))
    print("Molecules that passed the dissimilarity (structural) filters: " + str(len(smiles_seleccionados)))
    print("Of which: " + str(df_decoy_complete.shape[0]) + " are different")

    return (df_decoy_complete, fp_query, df_analysis)


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
        
    print("Finally, " + str(df_Tanimoto_final.shape[0]) + " passed the Tanimoto filter by comparing all loaded molecules")
    print('')
    print('Congratulations, you have obtained ' + str(final_decoys.shape[0]) + " decoys!!!")
    print("---------------------------------------------------------------")
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

#%% Settings
def setting_info():
    from datetime import date
    today = date.today()
    fecha = today.strftime("%d/%m/%Y")
    settings = []
    settings.append(["Decoys generated at: " , fecha])
    settings.append(["\n"])
    settings.append(["Physicochemical features limits:",""])
    settings.append(["MW" , "+/- " +  str(lim_MW)])
    settings.append(["logP" , "+/- " + str(lim_logP)])
    settings.append(["Num_rotable_bonds" , "+/- " + str(lim_rb)])
    settings.append(["Num_H_acceptors" , "+/- " + str(lim_Hba)])
    settings.append(["Num_H_donors" , "+/- " + str(lim_Hbd)])
    settings.append(["\n"])
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
    settings_df = pd.DataFrame(settings)
    return settings_df


#%%

lista_resultados = decoy_fase1(loaded_smiles, verbose = True)
print("Now, we are comparing all decoys vs all input SMILES, please wait a moment...")
df_final_decoys, df_analysis = duplicates_filter(lista_resultados[0], lista_resultados[1], lista_resultados[2]) 
settings_df = setting_info()  

# Here you can dowload the generated decoys
df_final_decoys.to_csv(f'{directory}\\Generated_decoys_{smiles_dataset}.csv',index=False,header=True)

# Here you can see a little analysis of the process
df_analysis.to_csv(f'{directory}\\Decoys_analysis_{smiles_dataset}.csv',index=True,header=True)

# Here you can download your settings
settings_df.to_csv(f'{directory}\\Decoys_settings_{smiles_dataset}.csv',index=False,header=False)

end = time.time()
hours, rem = divmod(end-start, 3600)
minutes, seconds = divmod(rem, 60)
print("{:0>2}:{:0>2}:{:05.2f}".format(int(hours),int(minutes),seconds)) 


