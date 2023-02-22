from pypath.utils import mapping
import numpy as np
import pandas as pd


# write function to flip item_id_a and item_id_b value if item_id_a starts with 9606
def flip_item_id(row):
    if row['item_id_a'].startswith('9606'):
        return row['item_id_b'], row['item_id_a']
    else:
        return row['item_id_a'], row['item_id_b']

# create columns in details dataframe, clipping the '9606.' prefix from the protein column
def clip_ensembl(row):
    return row['protein'].replace('9606.', '')

# writing function that takes a list of ensp ids and returns a list of gene symbols
def ensp_to_genesymbol(ensp_list):
    gene_symbol_list = []
    for element in ensp_list:
        symbol = mapping.map_name(element, 'ensp_biomart', 'genesymbol')
        if symbol != set():
            gene_symbol_list.append(symbol.pop())
        else:
            symbol = mapping.map_name(element, 'ensp', 'genesymbol')
            if symbol != set():
                gene_symbol_list.append(symbol.pop())
            else:
                gene_symbol_list.append('NA')
    return gene_symbol_list



#write function that takes a series of floats and converts them first to integers and then to strings
def float_to_string(series):
    series = series.astype(int)
    series = series.astype(str)
    return series

# convert a series of objects to strings and cut off the last two characters
def object_to_string(series):
    series = series.astype(str)
    series = series.str[:-2]
    return series

def get_hmdb_ids(df, metmap1, metmap2, metmap3):
    df = df.merge(metmap1, on='pubchem_id', how='left')
    df = df.merge(metmap2, on='pubchem_id', how='left')
    df = df.merge(metmap3, on='pubchem_id', how='left')
    return df


def drop_nan(df, col1, col2, col3):
    df = df.copy()
    df = df.dropna(subset=[col1, col2, col3], how='all')
    df[col1] = df[col1].fillna(df[col2])
    df[col1] = df[col1].fillna(df[col3])
    df = df.drop(col2, axis=1)
    df = df.drop(col3, axis=1)
    # if length is 9 insert two zeros after the first 4 characters
    df[col1] = np.where(df[col1].str.len() == 9, df[col1].str[:4] + '00' + df[col1].str[4:], df[col1])
    return df