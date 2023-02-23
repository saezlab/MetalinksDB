import sys
import os
import argparse
import configparser
import pandas as pd
import scipy.io as sio

this_script = os.path.dirname(sys.argv[0])
script_path = os.path.abspath(this_script)

config = configparser.ConfigParser()
config.read(f'{script_path}/config.ini')


### Set arguments
parser = argparse.ArgumentParser(
    description='Get Prodution Degration (PD) resource for CCC inference')

# required arguments
requiredArguments = parser.add_argument_group('required arguments')


# optional arguments
parser.add_argument('--recon_dir',
                    help='Directory of recon dataset',
                    required=False,
                    default=config['PD_paths']['recon_dir'])



args = parser.parse_args()
args_dict = vars(args)


recon = sio.loadmat(args_dict['recon_dir'])

print(recon)