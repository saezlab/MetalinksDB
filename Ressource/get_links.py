# Load standard packages
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Load needed files 
taf_gtp = pd.read_csv('~/Documents/Database/targets_and_families.csv', skiprows=1) # Guide to Pharmacology targets (GtP)
info = pd.read_csv('~/Documents/Python/global_ID_mapping_NA.tsv', dtype = str, index_col = 0, sep = '\t') # Kenos mapping table
hlinks = pd.read_csv('~/Documents/Database/protein_chemical_links_human.tsv', sep = '\t') # Stitch human chemical protein links

print('loading done')

# Cut out all catalytic, g-coupled and nuclear receptors from GtP
receptors = [x in ['catalytic_receptor', 'gpcr', 'nhr'] for x in taf_gtp['Type']]
receptor_df = taf_gtp.iloc[receptors, [0,2,11,12,15,16,18]]
receptor_df.index = receptor_df['HGNC symbol']

# Map symbols to ENSP IDs (here done via the web interface, but available via API)

# Load ID mapping and clean ENSP IDs
uniprot_to_ENSP = pd.read_csv('~/Documents/Startover/uniprot_to_ENSP.tsv', sep = '\t', index_col=0)
uniprot_to_ENSP['clean'] = [x[:15] for x in uniprot_to_ENSP['To']]

# Load Stitch-Chebi mappings (also done via uniprot web-interface)
chebi1 = pd.read_csv('~/Documents/Database/IDs1.txt', sep = '\t', names = ['CID', 'chebi'])
chebi2 = pd.read_csv('~/Documents/Database/IDs2.txt', sep = '\t', names = ['CID', 'chebi'])
chebi = pd.concat([chebi1, chebi2])

# Prepare a dictionary and map the IDs
chebi_dict = dict(zip(chebi.iloc[:,0], chebi.iloc[:,1]))
CIDs = hlinks.iloc[:,0]
CIDs_ints = [int(x[4:]) for x in CIDs]
CIDs_chebi = [chebi_dict.get(x) for x in CIDs_ints]
hlinks['chebi_ID'] = CIDs_chebi
hlinks = hlinks.dropna() #Remove links with NAs
hlinks['ENSP'] = [x[5:] for x in hlinks.iloc[:,1]]

print('mapping done')

# Reorder columns
hlinks.index = hlinks.iloc[:,1]
hlinks = hlinks.iloc[:,2:9]

# List ENSP IDs (there is shurely a more efficient way)
hlinks_ensp = list(hlinks['ENSP']) 
receptor_ensp = list(uniprot_to_ENSP['clean'])

# Exclude all connection that aren't connected to a receptor
receptor_bool = [x in receptor_ensp for x in hlinks['ENSP']]
hlinks_receptor = hlinks.iloc[receptor_bool,:]

# Exclude connections that are solely based on text mining
textmin = list(hlinks_receptor.iloc[:,0:3].sum(axis=1) > 0)
hlinks_receptor_notextmin = hlinks_receptor.iloc[textmin,:]

# Create dictionaries to add metadata (also inefficient)
chebi_dict = dict(zip(info.iloc[:,1], info.index))
uniprot_dict = dict(zip(uniprot_to_ENSP.iloc[:,1], uniprot_to_ENSP.index))
symbol_dict = dict(zip(receptor_df.iloc[:,5], receptor_df.iloc[:,2]))
name_dict = dict(zip(receptor_df.iloc[:,5], receptor_df.iloc[:,3]))
type_dict = dict(zip(receptor_df.iloc[:,5], receptor_df.iloc[:,0]))

# Add metadata (inefficient)
hlinks_receptor_notextmin['name_keno'] = [chebi_dict.get(x) for x in hlinks_receptor_notextmin['chebi_ID']]
hlinks_receptor_notextmin['uniprot'] = [uniprot_dict.get(x) for x in hlinks_receptor_notextmin['ENSP']]
hlinks_receptor_notextmin['target_symbol'] = [symbol_dict.get(x) for x in hlinks_receptor_notextmin['uniprot']]
hlinks_receptor_notextmin['target_name'] = [name_dict.get(x) for x in hlinks_receptor_notextmin['uniprot']]
hlinks_receptor_notextmin['target_type'] = [type_dict.get(x) for x in hlinks_receptor_notextmin['uniprot']]

# Restrict connections to metabolites contained in kenos table
subs_bool = [x is not None for x in hlinks_receptor_notextmin['name_keno']]
subs = hlinks_receptor_notextmin.iloc[subs_bool,:]

subs.to_csv('test.csv')

print('DONE!')


