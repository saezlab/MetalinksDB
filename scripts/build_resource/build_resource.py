import sys
import os
import argparse
import configparser
import pandas as pd

this_script = os.path.dirname(sys.argv[0])
script_path = os.path.abspath(this_script)

config = configparser.ConfigParser()
config.read(f'{script_path}/config.ini')


### Set arguments
parser = argparse.ArgumentParser(
    description='Infer signalling link between single cells or bulk samples')

# required arguments
requiredArguments = parser.add_argument_group('required arguments')


# optional arguments
parser.add_argument('--STITCH_actions_dir',
                    help='Directory of Stitch dataset',
                    required=False,
                    default=config['DefaultPaths']['STITCH_actions_dir'])


parser.add_argument('--STITCH_detailed_dir',
                    help='blabla',
                    required=False,
                    default=config['DefaultPaths']['STITCH_detailed_dir'])

parser.add_argument('--GtP_dir',
                    help='guide to pharma resource for CCC inference',
                    required=False,
                    default=config['DefaultPaths']['GtP_dir'])

parser.add_argument('--met_map_dir',
                    help='',
                    required=False,
                    default=config['DefaultPaths']['met_map_dir'])

parser.add_argument('--met_map_keno_dir',
                    help='',
                    required=False,
                    default=config['DefaultPaths']['met_map_keno_dir'])

parser.add_argument('--mode_args',
                    help='Modes of action to include in the resource',
                    required=False,
                    default = ['inhibition', 'activation'], 
                    nargs='+',
                    choices = ['activation', 'binding', 'catalysis', 'expression', 'inhibition', 'pred_bind', 'reaction'])

parser.add_argument('--out_dir',
                    help='Metabolite-protein link resource for CCC inference',
                    required=False,
                    default=config['DefaultPaths']['out_dir'])

parser.add_argument('--confidence_cutoffs',
                    help='Confidence cutoffs for thresholding interactions',
                    required=False,
                    default=[150, 150, 700],
                    nargs=3,
                    type=int)

parser.add_argument('--aux_dir',
                    help='Directory of auxillary scripts',
                    required=False,
                    default=config['DefaultPaths']['aux_dir'])

parser.add_argument('--receptor_args',
                    help='Receptors types to include in the resource',
                    required=False,
                    default = ['catalytic receptor', 'gpcr', 'nhr', 'other_protein'],
                    nargs='+',
                    choices = ['catalytic receptor', 'enzyme', 'gpcr', 'lgic', 'nhr', 'other_ic','other_protein', 'transporter', 'vgic'])

args = parser.parse_args()
args_dict = vars(args)

sys.path.append(args_dict['aux_dir'])
from aux import *

# load stitch data
actions = pd.read_csv(args_dict['STITCH_actions_dir'], sep='\t')
details = pd.read_csv(args_dict['STITCH_detailed_dir'], sep='\t')

print(f"Loaded {actions.shape[0]} actions and {details.shape[0]} details")

# cut down dataframe actions to include only the chosen modes of action
actions = actions[actions['mode'].isin(args_dict['mode_args'])]

# apply function to actions dataframe
actions['item_id_a'], actions['item_id_b'] = zip(*actions.apply(flip_item_id, axis=1))

# rename item_id_a and item_id_b to chemical and protein
actions = actions.rename(columns={'item_id_a': 'chemical', 'item_id_b': 'protein'})

# reduce details dataframe to only include tuples of chemical and protein that are in actions dataframe
details = details[details['chemical'].isin(actions['chemical']) & details['protein'].isin(actions['protein'])]

print(f"After mode filtering, {actions.shape[0]} actions and {details.shape[0]} detailed interactions remain")

del actions

# create column in details clipping the first four characters from the chemical column
details['pubchem_id'] = details['chemical'].str[4:]

metmap1 = pd.read_csv(args_dict['met_map_keno_dir'], sep='\t')
metmap2 = pd.read_csv(args_dict['met_map_dir'], sep=';')

details['pubchem_id'] = float_to_string(details['pubchem_id'])

metmap1['pubchem_id'] = object_to_string(metmap1['pubchem_id'])
metmap2['CID'] = object_to_string(metmap2['CID'])

metmap2 = metmap2.rename(columns={'CID': 'pubchem_id'})

merged_table = get_hmdb_ids(details, metmap1, metmap2)

merged_table = drop_nan(merged_table, 'hmdb_id', 'HMDB')

# cut down merged table to only include columns that are in details and hmdb_id
merged_table = merged_table[['hmdb_id','protein', 'database', 'experimental', 'textmining','combined_score', 'pubchem_id', 'Name']].drop_duplicates()

details = merged_table

print(f"After merging with HMDB, {details.shape[0]} interactions remain")


details['ensp'] = details.apply(clip_ensembl, axis=1)


# convert ensp ids to gene symbols and add to details dataframe
details['symbol'] = ensp_to_genesymbol(details['ensp'])

# remove rows with gene symbols that are 'NA'
details = details[details['symbol'] != 'NA']

print(f"After gene symbol conversion, {details.shape[0]} interactions remain")


gtp = pd.read_csv(args_dict['GtP_dir'], sep=',', skiprows=1)

poi = args_dict['receptor_args']

gtp_cut = gtp[gtp['Type'].isin(poi)]
receptors = gtp_cut['HGNC symbol']

details = details[details['symbol'].isin(receptors)]

print(f"After filtering for receptors, {details.shape[0]} interactions remain")

# further cut down details dataframe to only include rows with a database score higher than 150, experimental score higher than 150, and a textmining score higher than 700
details = details[(details['database'] > args_dict['confidence_cutoffs'][0]) | (details['experimental'] > args_dict['confidence_cutoffs'][1]) | (details['textmining'] > args_dict['confidence_cutoffs'][2])]

print(f"After score filtering, {details.shape[0]} interations remain")


details.to_csv(args_dict['out_dir'], sep='\t', index=False)


print('done')





















