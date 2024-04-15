import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import gridspec
import matplotlib.colors as mcolors
from scipy.stats import ttest_ind

def load_prepro_metalinks(anno_path = '/Users/ef6/Documents/Saez/metalinks/Data/Intermediate/HMDB/hmdb_metabolites_explained.csv', 
                          MR_path = '/Users/ef6/Documents/GitHub/metalinks_analysis/metalinksDB/MR_500500900.csv', 
                          PD_path = '/Users/ef6/Documents/GitHub/metalinks_analysis/metalinksDB/PD.csv'):
    df = pd.read_csv(anno_path)
    MR_original = pd.read_csv(MR_path, sep=',')
    # remove " from HMDB and symbol column
    MR_original['HMDB'] = MR_original['HMDB'].str.replace('"', '')
    MR_original['Symbol'] = MR_original['Symbol'].str.replace('"', '')

    PD_original = pd.read_csv(PD_path)

    out = df[(df['kingdom'] == 'Inorganic compounds') & (df['accession'].isin(PD_original['HMDB']))]['accession']
    PD = PD_original[PD_original['HMDB'].isin(out) == False]

    return PD_original, MR_original, df, PD

def load_prepro_other_dbs(MR_original, PD_original,
                          cellphone_path = '/Users/ef6/Documents/GitHub/metalinks/data/CellphoneDB/Cellphone_suptab4_curated.xlsx',
                          neuronchat_path = '/Users/ef6/Documents/Saez/metalinks/Data/Source/Other_DBs/NeuronChatDB_human.csv',
                          mebocost_PD_path = '/Users/ef6/Documents/Saez/metalinks/Data/Source/Other_DBs/MebocostPD.tsv',
                          mebocost_MR_path = '/Users/ef6/Documents/Saez/metalinks/Data/Source/Other_DBs/MebocostDB.tsv',
                          neuronchat_table_path = '/Users/ef6/Documents/Saez/metalinks/Data/Intermediate/Mapping/Neuronchat_table.csv'):
    cpdb = pd.read_excel(cellphone_path)
    metabolites = cpdb.iloc[:,2]
    sensors = cpdb.iloc[:,3]
    sensors = [x.split('_')[0] for x in sensors]
    unique_metabolites_cpdb = np.unique(metabolites)
    unique_sensors_cpdb = np.unique(sensors)

    ncdb = pd.read_csv(neuronchat_table_path, sep=',')
    unique_metabolites_ncdb = np.unique(ncdb['HMDB'])

    ncdb_cut = pd.read_csv(neuronchat_path, sep=',')
    ncdb_cut['Sensor'] = ncdb_cut['interaction_name'].str.split('_').str[1]
    unique_sensors_ncdb = np.unique(ncdb_cut['Sensor'])

    ncdb_dict = dict(zip(ncdb['Query'], ncdb['HMDB']))
    ncdb_cut['Query'] = ncdb_cut['interaction_name'].str.split('_').str[0]
    ncdb_cut['HMDB'] = ncdb_cut['Query'].map(ncdb_dict)

    ncdb_cut['gene'] = ncdb_cut['interaction_name'].str.split('_').str[1]

    mdb = pd.read_csv(mebocost_PD_path, sep = '\t')
    unique_metabolites_mdb = mdb['HMDB_ID'].unique()
    l = []
    for i in unique_metabolites_mdb:
        # append i to l if the direction column of mdb subset by i contains 'substrate' and 'product'
        if 'substrate' in mdb[mdb['HMDB_ID'] == i]['direction'].values and 'product' in mdb[mdb['HMDB_ID'] == i]['direction'].values:
            l.append(i)
    unique_metabolites_mdb = l

    unique_metabolites_mldb = np.unique(PD_original['HMDB'])

    # MR  
    mdb = pd.read_csv(mebocost_MR_path, sep = '\t')
    unique_metabolites_mdb_MR = mdb['HMDB_ID'].unique()
    unique_sensors_mdb = np.unique([x.split('[')[0] for x in  mdb.iloc[:,3]])

    unique_metabolites_mldb_MR = np.unique(MR_original['HMDB'])
    unique_sensors_mldb = np.unique(MR_original['Symbol'])

    int_ncdb = np.unique([x + '_' + y for x,y in zip(ncdb_cut['HMDB'], ncdb_cut['gene'])])
    int_cpdb = np.unique([x + '_' + y for x,y in zip(metabolites, sensors)])
    int_mdb = np.unique([x + '_' + y for x,y in zip(mdb['HMDB_ID'], [x.split('[')[0] for x in  mdb.iloc[:,3]])])
    int_mldb = np.unique([x + '_' + y for x,y in zip(MR_original['HMDB'], MR_original['Symbol'])])

    data1 = [unique_metabolites_ncdb, unique_metabolites_cpdb, unique_metabolites_mdb, unique_metabolites_mldb] #PD
    data2 = [unique_metabolites_ncdb, unique_metabolites_cpdb, unique_metabolites_mdb_MR, unique_metabolites_mldb_MR] #MR met
    data3 = [unique_sensors_ncdb, unique_sensors_cpdb, unique_sensors_mdb, unique_sensors_mldb] #MR gene
    data4 = [int_ncdb, int_cpdb, int_mdb, int_mldb]

    return data1, data2, data3, data4


def preprocess_data_for_barplot(data, database_names, index_names):
    df = pd.DataFrame({database_names[0]: data[0], 
                       database_names[1]: data[1],
                       database_names[2]: data[2],
                       database_names[3]: data[3],
                       database_names[4]: data[4],
                       database_names[5]: data[5]}, index=index_names)
    df = df.reset_index()
    df = pd.melt(df, id_vars=['index'], value_vars=database_names)
    df.columns = ['index', 'Database', 'Count']
    df = df.sort_values(by=['index'])
    df['Database'] = pd.Categorical(df['Database'], categories=database_names)
    df = df.sort_values(by=['Database'])
    df['index'] = pd.Categorical(df['index'], categories=index_names)
    df = df.sort_values(by=['index'])
    df['Count'] = df['Count'].astype(int)
    return df

def mr_barplot(df):
    p = (ggplot(df, aes(x='index', y='Count', fill='Database')) +
        geom_bar(stat='identity', position='dodge') +
        theme_classic() +
        theme(axis_text_x=element_text(rotation=0)) +
        labs(x='', y='Count', fill='Database') +
        # make a logarithmic scale
        scale_y_log10()
    )
    return p

def mr_barplot2(df):
    plt.figure(figsize=(10, 6))
    sns.barplot(data=df, x='index', y='Count', hue='Database')
    plt.xlabel('')
    plt.ylabel('Count')
    plt.legend(title='Database')
    plt.yscale('log')
    plt.xticks(rotation=0)
    sns.despine()
    plt.show()


def prepare_fractions(data, df):
    nc_df = df[df['accession'].isin(data[0])]
    cp_df = df[df['accession'].isin(data[1])]
    cl_df = df[df['accession'].isin(data[2])]
    me_df = df[df['accession'].isin(data[3])]
    sc_df = df[df['accession'].isin(data[4])]
    ml_df = df[df['accession'].isin(data[5])]

    classes = df.columns[6:10]
    res = []
    for met_class in classes:
        fractions = pd.DataFrame()
        fractions['MetalinksDB'] = ml_df[met_class].value_counts(normalize=True)
        fractions['scConnect'] = sc_df[met_class].value_counts(normalize=True)
        fractions['NeuronChat'] = nc_df[met_class].value_counts(normalize=True)
        fractions['CellphoneDB'] = cp_df[met_class].value_counts(normalize=True)
        fractions['Cellinker'] = cl_df[met_class].value_counts(normalize=True)
        fractions['MebocostDB'] = me_df[met_class].value_counts(normalize=True)
        s = fractions.sum(axis=1)
        fractions = fractions[s > 0.09]
        colsums = fractions.sum(axis=0)
        others = pd.DataFrame({'MetalinksDB': 1 - colsums['MetalinksDB'], 
                                'scConnect': 1 - colsums['scConnect'],
                               'NeuronChat': 1 - colsums['NeuronChat'],
                               'CellphoneDB': 1 - colsums['CellphoneDB'], 
                                'Cellinker': 1 - colsums['Cellinker'],                               
                               'MebocostDB': 1 - colsums['MebocostDB']
                               }, index=['Others'])
        fractions = pd.concat([fractions, others])
        fractions = fractions.fillna(0)
        res.append(fractions)

    return res


def figure_2(df, hm, values, matrix1, matrix2, subject='Enzyme sets'):
    labels = ['NeuronChat', 'CellphoneDB', 'MebocostDB', 'MetalinksDB']
    row_labels = ['NCDB', 'CPDB', 'MDB', 'MLDB']
    col_labels = ['NCDB', 'CPDB', 'MDB', 'MLDB']
    data1 = np.array(matrix1)
    data2 = np.array(matrix2)
    
    fig = plt.figure(figsize=(12, 8))
    gs = gridspec.GridSpec(2, 1, height_ratios=[2,2])
    gs_upper = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[0, 0], wspace=0.1, width_ratios=[1, 1, 1])
    gs_lower = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1, 0], wspace=0.2, width_ratios=[1, 2, 1])

    colors = ['#B2C9AB', '#92B6B1', '#788AA3', '#932A61']
    cmap_custom = mcolors.LinearSegmentedColormap.from_list('custom_cmap', ['#512D55', '#B2C9AB'])

    
    ax1 = plt.subplot(gs_upper[ 0])
    sns.barplot(data=df, x='index', y='Count', hue='Database', ax=ax1, palette=colors)
    ax1.set_yscale('log')
    ax1.set_xlabel('')
    ax1.set_ylabel('')
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax1.spines['left'].set_linewidth(2)
    ax1.spines['bottom'].set_linewidth(2)
    ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # remove box around the legend
    ax1.get_legend().get_frame().set_linewidth(0.0)

    ax6 = plt.subplot(gs_upper[1])
    # make a white plain for ax6
    ax6.axis('off') 

    ax2 = plt.subplot(gs_upper[0, 2])
    sns.heatmap(hm, annot=True, cmap='cividis', fmt='.2f', vmin=0, vmax=1, ax=ax2 )
    ax2.set_xticklabels(ax2.get_xticklabels(), rotation=0)
    cbar = ax2.collections[0].colorbar
    cbar.set_label('Fraction')

    ax3 = plt.subplot(gs_lower[0])
    sns.barplot(y=labels, x=values, orient='h', ax=ax3, palette=colors)
    ax3.grid(False)
    ax3.set_facecolor('white')
    for i, value in enumerate(values):
        ax3.text(value + 200, i, str(value), va='center', color='black')
    ax3.spines['left'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_linewidth(2)
    ax3.spines['bottom'].set_linewidth(2)
    ax3.spines['bottom'].set_bounds(0, max(values))
    ax3.invert_xaxis()
    ax3.set_yticklabels([])
    ax3.set_yticks([])
    ax3.set_ylabel(subject)

    ax4 = plt.subplot(gs_lower[1])
    sns.heatmap(data1, annot=True, cmap='cividis', fmt='.2f', xticklabels=col_labels, yticklabels=row_labels, ax=ax4)
    ax4.set_xticklabels(ax4.get_xticklabels(), rotation=0)
    ax4.set_yticklabels(ax4.get_yticklabels(), rotation=0)
    cbar1 = ax4.collections[0].colorbar
    cbar1.set_label('Jaccard Index', size=12)
    ax4.spines['right'].set_visible(False)
    ax4.spines['top'].set_visible(False)
    ax4.spines['left'].set_linewidth(2)
    ax4.spines['bottom'].set_linewidth(2)

    ax5 = plt.subplot(gs_lower[2])
    sns.heatmap(data2, annot=True, cmap='cividis', fmt='.2f', xticklabels=['shared', 'unique'], yticklabels=row_labels, vmin=0, vmax=1, ax=ax5)
    ax5.set_xticklabels(ax5.get_xticklabels(), rotation=0)
    ax5.set_yticklabels([])
    ax5.set_yticks([])
    cbar2 = ax5.collections[0].colorbar
    cbar2.set_label('Fraction', size=12)
    ax5.spines['right'].set_visible(False)
    ax5.spines['top'].set_visible(False)
    ax5.spines['left'].set_linewidth(2)
    ax5.spines['bottom'].set_linewidth(2)

    ax1.text(-0.07, 1.08, 'a', transform=ax1.transAxes, fontsize=20, fontweight='bold', va='top', ha='right')
    ax2.text(-1.00, 1.08, 'b', transform=ax2.transAxes, fontsize=20, fontweight='bold', va='top', ha='right')
    ax3.text(-0.02, 1.08, 'c', transform=ax3.transAxes, fontsize=20, fontweight='bold', va='top', ha='right')
    ax4.text(-0.02, 1.08, 'd', transform=ax4.transAxes, fontsize=20, fontweight='bold', va='top', ha='right')
    ax5.text(-0.07, 1.08, 'e', transform=ax5.transAxes, fontsize=20, fontweight='bold', va='top', ha='right')


    plt.tight_layout()
    plt.show()

def prepro_hm(hm):

    hm = hm.fillna(0).iloc[:, ::-1]
    hm = hm[['NeuronChat', 'CellphoneDB','Cellinker', 'MebocostDB', 'scConnect', 'MetalinksDB']]
    hm.columns = ['NeuronChat', 'CellphoneDB','Cellinker', 'MebocostDB', 'scConnect', 'MetalinksDB']
    # hm.index = ['Organoheterocyclic compounds', 'Benzenoids',
    #     'Lipids and lipid-like molecules', 'Organic acids and derivatives',
    #     'Organic oxygen compounds', 'Nucleosides, nucleotides',
    #     'Organic nitrogen compounds', 'Homogeneous non-metal compounds',
    #     'Others']
    return hm

def count_shared_and_unique_metabolites(hmdb_ids):
    matrix_size = len(hmdb_ids)
    count_matrix = np.zeros((matrix_size, 2))

    for i in range(matrix_size):
        shared_metabolites = set(hmdb_ids[i])
        unique_metabolites = set(hmdb_ids[i])

        for j in range(matrix_size):
            if i != j:
                shared_metabolites = shared_metabolites.intersection(hmdb_ids[j])
                unique_metabolites = unique_metabolites.difference(hmdb_ids[j])

        count_matrix[i, 0] = len(shared_metabolites)
        count_matrix[i, 1] = len(unique_metabolites)

        #convert count matrix to integers
        count_matrix = count_matrix / len([item for sublist in hmdb_ids for item in sublist])

    return count_matrix

def get_jaccard(hmdb_ids):
    matrix_size = len(hmdb_ids)
    jaccard_matrix = np.zeros((matrix_size, matrix_size))

    for i in range(matrix_size):
        for j in range(i, matrix_size):
            set1 = set(hmdb_ids[i])
            set2 = set(hmdb_ids[j])
            intersection = len(set1.intersection(set2))
            union = len(set1.union(set2))
            jaccard_index = intersection / union
            jaccard_matrix[i, j] = jaccard_index
            jaccard_matrix[j, i] = jaccard_index

    # Print the Jaccard matrix
    return(jaccard_matrix)


def cosmos_DE_before(MR_path = '/Users/ef6/Documents/GitHub/metalinks_analysis/metalinksDB/MR_500500900_Kidney-pred.csv'):
    metabolites = pd.read_csv('/Users/ef6/Documents/Saez/metalinks/Data/Intermediate/COSMOS/metab_ttop_tumour_vs_healthy.csv', index_col=0)
    RNA = pd.read_csv('/Users/ef6/Documents/Saez/metalinks/Data/Intermediate/COSMOS/RNA_ttop_tumorvshealthy.csv', index_col=0)
    COSMOS_mapping = pd.read_csv('/Users/ef6/Documents/Saez/metalinks/Data/Intermediate/Mapping/ocean_mets.csv', index_col=0)
    MR = pd.read_csv(MR_path)
    MR['HMDB'] = MR['HMDB'].str.replace('"', '')
    MR['Symbol'] = MR['Symbol'].str.replace('"', '')

    RNA_IDS = pd.read_csv('/Users/ef6/Documents/Saez/metalinks/Data/Intermediate/Mapping/COSMOS_RNA.csv', index_col=0)
    RNA_IDS = RNA_IDS[RNA_IDS['SYMBOL'].isin(MR['Symbol'])]
    metabolites = metabolites.join(COSMOS_mapping, how='left')
    metabolites = metabolites[metabolites['HMDB'].isin(MR['HMDB'])]
    RNA = RNA[RNA.index.isin(RNA_IDS['ENTREZID'])]
    RNA_dict = dict(zip(RNA_IDS['ENTREZID'], RNA_IDS['SYMBOL']))
    RNA.index = [RNA_dict.get(x) for x in RNA.index]
    # cut MR to only include hmdb id that are in metabolites['HMDB'] and RNA.ndex
    MR_cut = MR[(MR['HMDB'].isin(metabolites['HMDB'])) & (MR['Symbol'].isin(RNA.index))]
    RNA['Symbol'] = np.array(RNA.index)
    # rename HMDB column to hmdb_id
    metabolites.index = metabolites['HMDB']
    # left join RNA and metabolites to MR
    df = MR_cut.join(RNA, on= 'Symbol', how='left', lsuffix='_MR', rsuffix='_RNA')
    # join metabolites and df on the hmdb_id columns
    df.drop(columns=['Symbol_RNA'], inplace=True)
    df.rename(columns={'Symbol_MR': 'symbol'}, inplace=True)

    df = df.join(metabolites, on='HMDB',  lsuffix='_RNA', rsuffix='_metabolites')

    df['metalinks'] = df['t_RNA'] + df['t_metabolites']
    # order after metalinks
    df.sort_values(by='metalinks', ascending=False, inplace=True)

    return df




def metalinks_DE():

    # loading and preprocessing
    COSMOS_RNA = pd.read_csv('/Users/ef6/Documents/GitHub/metalinks_benchmark/data/RNA_counts_vsn.tsv',  sep = '\t', index_col=0).T
    COSMOS_MET = pd.read_csv('/Users/ef6/Documents/GitHub/metalinks_benchmark/data/raw_metabolomic_vsn.csv', sep = ',', index_col=0).T
    COSMOS_mapping = pd.read_csv('/Users/ef6/Documents/Saez/metalinks/Data/Intermediate/Mapping/ocean_mets.csv', index_col=0)
    COSMOS_MET.columns = COSMOS_mapping['HMDB']

    #decide on how to treat zero values
    # COSMOS_MET.fillna(COSMOS_MET.mean(axis=1).mean(), inplace=True)

    COSMOS_MET.dropna(axis=1, inplace=True)
    COSMOS_MET.index = COSMOS_MET.index.str.replace('KI', 'H')
    COSMOS_MET.index = COSMOS_MET.index.str.replace('TU', 'T')
    COSMOS_MET = COSMOS_MET.sort_index()
    COSMOS_RNA = COSMOS_RNA[COSMOS_RNA.index.isin(COSMOS_MET.index)].sort_index()
    COSMOS_MET = COSMOS_MET.T
    COSMOS_RNA = COSMOS_RNA.T
    MR_path = '/Users/ef6/Documents/BC/TestDBs/Kidney_after_DB_F.csv'
    MR = pd.read_csv(MR_path)
    MR['HMDB'] = MR['HMDB'].str.replace('"', '')
    MR['Symbol'] = MR['Symbol'].str.replace('"', '')
    COSMOS_RNA = COSMOS_RNA[COSMOS_RNA.index.isin(MR['Symbol'])]

    # scale values in matrix per column.sum(axis=1)
    COSMOS_RNA = COSMOS_RNA.apply(lambda x: (x - np.mean(x)) / np.std(x), axis=1)
    COSMOS_MET = COSMOS_MET[COSMOS_MET.index.isin(MR['HMDB'])]
    COSMOS_MET = COSMOS_MET.apply(lambda x: (x - np.mean(x)) / np.std(x), axis=1)

    # cut MR to only include hmdb id that are in metabolites['HMDB'] and RNA.ndex
    MR_cut = MR[(MR['HMDB'].isin(COSMOS_MET.index)) & (MR['Symbol'].isin(COSMOS_RNA.index))]
    df = MR_cut.join(COSMOS_RNA, on= 'Symbol', how='left', lsuffix='_MR', rsuffix='_RNA')
    df = df.join(COSMOS_MET, on='HMDB',  lsuffix='_RNA', rsuffix='_metabolites')
    df.index = df['Symbol']

    # calculate metalinks scores by taking an average of both values
    test = (np.array(df.iloc[:,14:30]) + np.array(df.iloc[:,30:]))/2
    meta = pd.DataFrame( test) 
    meta.index = df.index
    meta = pd.concat([df.iloc[:,:3], meta], axis=1)
    meta.columns = ['HMDB', 'Symbol', 'MetName'] + list(COSMOS_RNA.columns) 
    meta.index = meta['MetName'] + '_' + meta['Symbol']

    healthy = meta[meta.columns[3:][meta.columns[3:].str.endswith('H')]]
    tumor = meta[meta.columns[3:][meta.columns[3:].str.endswith('T')]]

    ttest = pd.DataFrame(index=meta.index)
    ttest['ttest'] = ttest_ind(healthy, tumor, axis=1)[0]
    ttest['pvalue'] = ttest_ind(healthy, tumor, axis=1)[1]
    ttest['logFC'] = np.log2(np.mean(tumor, axis=1) / np.mean(healthy, axis=1))
    res = pd.concat([meta, ttest], axis=1)
    df.index = res.index # quick and dirty
    res['Pathways'] = df['m.pathways']
    res['Diseases'] = df['m.diseases']
    res.sort_values(by='pvalue', ascending=True, inplace=True)
    res['Symbol'] = res['Symbol'].str.replace('"', '')
    res['LR'] = res['Symbol'] + '_' + res['MetName']

    return res




# # Assuming you have a DataFrame called 'data' with columns 'LR', 'ttest', and 'pvalue'
# data = res[['LR', 'ttest', 'pvalue']]

# data = data[data['pvalue'] < 0.005]

# # Sort the DataFrame by absolute t-statistic in descending order
# data_sorted = data.assign(abs_tstatistic=data['ttest'].abs()).sort_values(by='abs_tstatistic', ascending=True)

# # Set the plot size and create subplots
# plt.figure(figsize=(8, 6))
# ax = plt.subplot()

# # Plot the lollipops
# col = '#932A61'
# ax.hlines(y=data_sorted['LR'], xmin=0, xmax=data_sorted['ttest'], color=col, alpha=0.7, linewidth=2)
# ax.scatter(data_sorted['ttest'], data_sorted['LR'], color=col, alpha=0.7, s=100 * (1 - data_sorted['pvalue']))

# # Customize the plot
# ax.set_xlabel('t-statistic')
# ax.axvline(x=0, color='gray', linestyle='--', linewidth=1)
# ax.grid(True, axis='y', linestyle='--')

# # Show the plot
# plt.tight_layout()
# plt.show()