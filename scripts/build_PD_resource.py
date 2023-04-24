import sys
import os
import argparse
import configparser
import pandas as pd
import scipy.io as sio

this_script = os.path.dirname(sys.argv[0])
script_path = os.path.abspath(this_script)

config = configparser.ConfigParser()
config.read(f'{script_path}/../config.ini')


### Set arguments
parser = argparse.ArgumentParser(
    description='Get Prodution Degration (PD) resource for CCC inference')

# required arguments
requiredArguments = parser.add_argument_group('required arguments')


# optional arguments
parser.add_argument('--recon_dir',
                    help='Directory of recon dataset',
                    required=False,
                    default=config['DefaultPaths']['recon_dir'])

parser.add_argument('--met_map_dir',
                    help='',
                    required=False,
                    default=config['DefaultPaths']['met_map_dir'])

parser.add_argument('--met_map_keno_dir',
                    help='',
                    required=False,
                    default=config['DefaultPaths']['met_map_keno_dir'])

parser.add_argument('--met_map_hmdb_dir',
                    help='',
                    required=False,
                    default=config['DefaultPaths']['met_map_hmdb_dir'])

parser.add_argument('--recon_mapping_dir',
                    help='',
                    required=False,
                    default=config['DefaultPaths']['recon_mapping_dir'])

parser.add_argument('--recon_symbols_dir',
                    help='',
                    required=False,
                    default=config['DefaultPaths']['recon_symbols_dir'])

parser.add_argument('--aux_dir',
                    help='',
                    required=False,
                    default=config['DefaultPaths']['aux_dir'])

parser.add_argument('--out_dir',
                    help='PD resource for CCC inference',
                    required=False,
                    default=config['DefaultPaths']['out_dir'])




args = parser.parse_args()
args_dict = vars(args)

sys.path.append(args_dict['aux_dir'])
from aux import *


recon = sio.loadmat(args_dict['recon_dir'])
symbols = pd.read_csv(args_dict['recon_symbols_dir']   , sep=';')
recon = recon['Recon3D']

print(f"Successfully loaded Recon3D model from {args_dict['recon_dir']}")
data = recon


rxn_gene_df = pd.DataFrame(data['rxnGeneMat'][0][0])
reaction_ids = data['rxns'][0][0].flatten()
reaction_ids = [x[0] for x in reaction_ids]
mets = data['mets'][0][0].flatten()
mets = [x[0] for x in mets]
rxn_gene_df.columns = symbols['symbols']
rxn_gene_df.index = reaction_ids
S = pd.DataFrame(data['S'][0][0].toarray(), index=mets, columns=reaction_ids)
lb_ub = pd.DataFrame(data['lb'][0][0], index=reaction_ids, columns=['lb'])
lb_ub['ub'] = data['ub'][0][0]
lb_ub['rev'] = lb_ub.apply(lambda x: 'reversible' if x['lb'] < 0 and x['ub'] > 0 else 'irreversible', axis=1)
lb_ub['direction'] = lb_ub.apply(lambda x: 'forward' if x['ub'] > 0 else 'backward', axis=1)

print(f'successfully extracted data from model')

reaction_to_genes = get_gene_symbols(rxn_gene_df)

reaction_to_metabolites_prod = get_metabolites(S, d = 1)
reaction_to_metabolites_deg = get_metabolites(S, d = -1)

metabolite_to_gene = get_metabolite_to_gene(reaction_to_metabolites_prod, reaction_to_metabolites_deg, reaction_to_genes, lb_ub)

print(f'collapsed metabolites to genes, now have {len(metabolite_to_gene)} metabolite to gene links')

metmap1 = pd.read_csv(args_dict['met_map_keno_dir'], sep='\t', dtype=object)
metmap2 = pd.read_csv(args_dict['met_map_dir'], sep='\t', dtype=str)
metmap3 = pd.read_csv(args_dict['met_map_hmdb_dir'], sep=',', dtype=str)
df =      pd.read_csv(args_dict['recon_mapping_dir'], sep=',', dtype=str)

print(f'loaded metabolite mapping files')

dfs = preprocess_metmaps(df, metmap1, metmap2, metmap3)

test = fill_missing_values(df1=dfs[0], df2=dfs[3], df3=dfs[1])
test = fill_missing_values(test[0], test[1], test[2], 'kegg_id', 'chebi_id', 'hmdb_id', 'pubchem_id')
test = fill_missing_values(test[0], test[1], test[2], 'hmdb_id', 'chebi_id', 'kegg_id', 'pubchem_id')
test = fill_missing_values(test[0], test[1], test[2], 'pubchem_id', 'chebi_id', 'kegg_id', 'hmdb_id')

print(f'filled missing values in metabolite mapping files')

met_dict = dict(zip(mets, test[0]['hmdb_id']))
metabolite_to_gene['hmdb_id'] = metabolite_to_gene['metabolite_id'].apply(lambda x: met_dict[x])

# drop metabolite and reaction ids
metabolite_to_gene.drop(['metabolite_id', 'reaction_id'], axis=1, inplace=True)
metabolite_to_gene.drop_duplicates(inplace=True)
metabolite_to_gene.dropna(subset=['hmdb_id'], inplace=True)


metabolite_to_gene.to_csv(f'{args_dict["out_dir"]}/metabolite_to_gene_test.csv', index=False)

print(metabolite_to_gene.head())








































