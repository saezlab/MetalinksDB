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
    df.drop_duplicates(inplace=True)
    df = df.merge(metmap2, on='pubchem_id', how='left')
    df.drop_duplicates(inplace=True)
    df = df.merge(metmap3, on='pubchem_id', how='left')
    df.drop_duplicates(inplace=True)
    return df



def get_hmdb_ids_s(df, metmap3):
    df = df.merge(metmap3, on='pubchem_id', how='left')
    df.drop_duplicates(inplace=True)
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


def preprocess_metmaps(df, metmap1, metmap2, metmap3):
    metmap2.rename(columns={'CID': 'pubchem_id', 'KEGG' : 'kegg_id', 'HMDB' : 'hmdb_id', 'ChEBI' : 'chebi_id'}, inplace=True)
    metmap2['chebi_id'] = 'CHEBI:' + metmap2['chebi_id']
    df.rename(columns={'metCHEBIID': 'chebi_id', 'metKEGGID': 'kegg_id', 'metHMDBID': 'hmdb_id', 'metPubChemID': 'pubchem_id'}, inplace=True)
    metmap1['pubchem_id'] = metmap1['pubchem_id'].apply(lambda x: str(int(x)) if not pd.isnull(x) else x)
    df['chebi_id'] = df['chebi_id'].apply(lambda x: 'CHEBI:' + x if not pd.isnull(x) and not x.startswith('CHEBI:') else x)

    # remove spaces from the beginning and end of the strings in df['kegg_id'] by extracting a substring starting with a 'C' and having 6 characters
    df['kegg_id'] = df['kegg_id'].apply(lambda x: x[x.find('C'):x.find('C')+6] if not pd.isnull(x) else x)
    # do the same for metmap2['kegg_id']
    metmap2['kegg_id'] = metmap2['kegg_id'].apply(lambda x: x[x.find('C'):x.find('C')+6] if not pd.isnull(x) else x)
    # change dtype of hmdb_id to str
    metmap2['hmdb_id'] = metmap2['hmdb_id'].astype(str)
    # if length of metmap2[str3 + '_id'] is lower that 11 add zeros after the first 4 characters to top up to 11 characters only if the string is not nan
    metmap2['hmdb_id'] = metmap2['hmdb_id'].apply(lambda x: x if len(x) == 11 else x[:4] + '0'*(11-len(x)) + x[4:])
    # set value in hmd_id to nan if it starts with 'nan
    metmap2['hmdb_id'] = metmap2['hmdb_id'].apply(lambda x: np.nan if x.startswith('nan') else x)
    # count length of strings in hmdb_id
    metmap2['hmdb_id'].str.len().value_counts()

    # do the same for df['hmdb_id']
    df['hmdb_id'] = df['hmdb_id'].astype(str)
    df['hmdb_id'] = df['hmdb_id'].apply(lambda x: x if len(x) == 11 else x[:4] + '0'*(11-len(x)) + x[4:])
    df['hmdb_id'] = df['hmdb_id'].apply(lambda x: np.nan if x.startswith('nan') else x)
    df['hmdb_id'].str.len().value_counts()

    # rename metmap3['accession'] to 'hmdb_id'
    metmap3.rename(columns={'accession': 'hmdb_id'}, inplace=True)
    # remove protein accesion column from metmap3
    metmap3.drop(columns=['protein_accession'], inplace=True)
    # add a 'CHEBI:' to the beginning of the strings in metmap3['chebi_id'] if they are not nan
    metmap3['chebi_id'] = metmap3['chebi_id'].apply(lambda x: 'CHEBI:' + x if not pd.isnull(x) else x)

    df1 = df[['chebi_id', 'kegg_id', 'hmdb_id', 'pubchem_id']]
    df2 = metmap1[['chebi_id', 'kegg_id', 'hmdb_id', 'pubchem_id']]
    df3 = metmap2[['chebi_id', 'kegg_id', 'hmdb_id', 'pubchem_id']]
    df4 = metmap3[['chebi_id', 'kegg_id', 'hmdb_id', 'pubchem_id']]

    return df1, df2, df3, df4


def fillna_with_map(df, str1, str2, str3, str4, dict1, dict2, dict3):
    df = df.copy()
    dictlist = [dict1, dict2, dict3]
    counter = 0
    for col in [str2, str3, str4]:
        df[col] = df[col].fillna(df[str1].map(dictlist[counter]))
        counter += 1
    return df


def check_fill(df):
    counter = []
    df = df.copy()
    # extract name of columns bi extracting the characters of the first column names until the first '_'
    str1 = df.columns[0][:df.columns[0].find('_')]
    for i in range(0, len(df)):
        row = df.iloc[i]
        # if row has one na value, check if other two rows have the same value, if yes, fill in the na value
        if row.isna().sum() <= 1:
            # if all three values are the same continue to next row
            if row[0] == row[1] and row[0] == row[2]:
                continue
            elif row[0] == row[1]:
                df.iloc[i,2] = row[0]
                continue
            elif row[0] == row[2]:
                df.iloc[i,1] = row[0]
                continue
            elif row[1] == row[2]:
                df.iloc[i,0] = row[1]
                continue
            else:
            # entry that is nan will be set to on of the other two values
                df.iloc[i, df.iloc[i].isna()] = row[~row.isna()][0]
                counter.append(row)
        # if row has two na values set both to the third value
        elif row.isna().sum() == 2:
            df.iloc[i, df.iloc[i].isna()] = row[~row.isna()]
        # if row has no na values, do nothing
    print(f' found {len(counter)} entries that had a conflict in type: {str1}')
    return df



# write function that takes in a dataframe and a list of columns and fills creates a dictionary with the values of the columns
# keys are the first columns and values are the values of the other columns, if the value is nan, it is not added to the dictionary
def create_dict(df, list_of_columns):
    df1 = df.dropna(subset=[list_of_columns[1]])
    keys = df1[list_of_columns[0]]
    values = df1[list_of_columns[1]]
    out_dict = dict(zip(keys, values))
    return out_dict

def fill_missing_values(df1, df2, df3, str1 = 'chebi_id', str2 = 'kegg_id', str3 = 'hmdb_id', str4 = 'pubchem_id'):
    df1 = df1.copy()
    # store unique values in each column of df1 in one list
    before = []
    for i in range(0, len(df1.columns)):
        before.append(df1.iloc[:,i].unique())

    dfx = df1.dropna(subset=[str1])
    # order the columns in dfx
    dfx = dfx[[str1, str2, str3, str4]]
    test = pd.merge(dfx, df2, how='left', on= str1).drop_duplicates()
    test = pd.merge(test, df3, how='left', on= str1).drop_duplicates()

    test.iloc[:,[1,4,7]] = check_fill(test.iloc[:,[1,4,7]])
    test.iloc[:,[2,5,8]] = check_fill(test.iloc[:,[2,5,8]])
    test.iloc[:,[3,6,9]] = check_fill(test.iloc[:,[3,6,9]])
    
    dict1 = create_dict(test, [str1, str2 + '_x'])
    dict2 = create_dict(test, [str1, str3 + '_x'])
    dict3 = create_dict(test, [str1, str4 + '_x'])

    df1 = fillna_with_map(df1, str1, str2, str3, str4, dict1, dict2, dict3) # still used the dictionaries but in the loop
    df2 = fillna_with_map(df2, str1, str2, str3, str4, dict1, dict2, dict3)
    df3 = fillna_with_map(df3, str1, str2, str3, str4, dict1, dict2, dict3)

    after = []
    for i in range(0, len(df1.columns)):
        after.append(df1.iloc[:,i].unique())

    for i in range(0, len(before)):
        print(f'before: {len(before[i])}, after: {len(after[i])}, additional: {len(set(after[i]) - set(before[i]))}')

    return df1, df2, df3


# write function that creates a dataframe with the reaction ids and associated metabolites
def get_metabolites(S):
    S = S.copy()
    S[S != 0] = 1
    S = S.stack().reset_index()
    S.columns = ['metabolite_id', 'reaction_id', 'value']
    S = S[S['value'] == 1]
    S = S.drop('value', axis=1)
    S.drop_duplicates(inplace=True)
    return S

# write function that creates a dataframe with the reaction ids and associated gene symbols
def get_gene_symbols(rxn_gene_df):
    row_sums = rxn_gene_df.sum(axis=1)
    rxn_genes = rxn_gene_df.index[row_sums > 0]
    rxn_gene_df = rxn_gene_df.loc[rxn_genes]
    rxn_gene_df = rxn_gene_df.stack().reset_index()
    rxn_gene_df.columns = ['reaction_id', 'gene_id', 'value']
    rxn_gene_df = rxn_gene_df[rxn_gene_df['value'] == 1]
    rxn_gene_df = rxn_gene_df.drop('value', axis=1)
    rxn_gene_df.drop_duplicates(inplace=True)
    return rxn_gene_df


def get_metabolites(S, d = 1):
    S = S.copy()
    S[S != d] = 0
    S[S == d] = 1
    S = S.stack().reset_index()
    S.columns = ['metabolite_id', 'reaction_id', 'value']
    S = S[S['value'] == 1]
    S = S.drop('value', axis=1)
    S.drop_duplicates(inplace=True)
    return S

# write function that creates a dataframe matching the metabolite_ids of reaction_to_metabolites_prod and reaction_to_metabolites_deg
#  to gene names in reaction_to_genes
# in this dataframe create a column that indicated whether the metabolite is produced or degraded

def get_metabolite_to_gene(reaction_to_metabolites_prod, reaction_to_metabolites_deg, reaction_to_genes, lb_ub):
    metabolite_to_gene = pd.merge(reaction_to_metabolites_prod, reaction_to_genes, on='reaction_id')
    metabolite_to_gene_deg = pd.merge(reaction_to_metabolites_deg, reaction_to_genes, on='reaction_id')
    metabolite_to_gene_deg['direction'] = 'degrading'
    metabolite_to_gene = pd.concat([metabolite_to_gene, metabolite_to_gene_deg])
    # name the entries in the direction column that are not degrading as producing
    metabolite_to_gene['direction'] = metabolite_to_gene['direction'].apply(lambda x: 'producing' if x != 'degrading' else x)
    # for reactions that are reversible we will add the other direction to the metabolite_to_gene dataframe
    reversible_reactions = lb_ub[lb_ub['rev'] == 'reversible'].index
    rev_df = metabolite_to_gene[metabolite_to_gene['reaction_id'].isin(reversible_reactions)]
    rev_df['direction'] = rev_df['direction'].apply(lambda x: 'degrading' if x == 'producing' else 'producing')
    metabolite_to_gene = pd.concat([metabolite_to_gene, rev_df])
    return metabolite_to_gene



































