import sys
import os
import argparse
import configparser
from tqdm import tqdm
import numpy as np
import pandas as pd


this_script = os.path.dirname(sys.argv[0])
script_path = os.path.abspath(this_script)

config = configparser.ConfigParser()
config.read(f'{script_path}/../config.ini')


### Set arguments
parser = argparse.ArgumentParser(
    description='Infer signalling link between single cells or bulk samples')

# required arguments
requiredArguments = parser.add_argument_group('required arguments')

# requiredArguments.add_argument('--format',
#                                help='Single cell or bulk RNAseq',
#                                required=True,
#                                default='bulk')


# optional arguments
parser.add_argument('--data_dir',
                    help='Directory of bulk dataset',
                    required=False,
                    default=config['DefaultPaths']['bulk_dir'])

parser.add_argument('--module_dir',
                    help='Directory of module metalinks',
                    required=False,
                    default=config['DefaultPaths']['module_dir'])

parser.add_argument('--pd_dir',
                    help='Production degradation ressource to use for metabolite abundance estimation',
                    required=False,
                    default=config['DefaultPaths']['pd_dir'])

parser.add_argument('--trans_dir',
                    help='Metabolite transportation ressource to use for metabolite abundance estimation',
                    required=False,
                    default=config['DefaultPaths']['trans_dir'])

parser.add_argument('--links_dir',
                    help='Metabolite-protein link resource for CCC inference',
                    required=False,
                    default=config['DefaultPaths']['links_dir'])

parser.add_argument('--out_dir',
                    help='Metabolite-protein link resource for CCC inference',
                    required=False,
                    default=config['DefaultPaths']['out_dir'])


args = parser.parse_args()
args_dict = vars(args)

sys.path.append(args_dict['module_dir'])
import metalinks as ml


# load bulk or scData
transcriptome = ml.ld.read_data(args_dict['data_dir'])


# estimate metabolite abundance
metabolome = ml.mi.infer_metab_abund(transcriptome, args_dict['pd_dir'], args_dict['trans_dir'])


# infer cell-cell communication
enriched_interactions = ml.ccc.infer_CCC(transcriptome, metabolome, args_dict['links_dir'])


#produce summary in plots
ml.vis.prod_sum_plots(enriched_interactions, args_dict['out_dir'])         





















