library(R.matlab)
library(stringr)
library(readr)
library(ocean)
library(org.Hs.eg.db)
library(annotate)
library(plyr)
library(dplyr)


setwd('~/Documents/Database_old/recon3D_full/')
recon3D <- readMat("~/Documents/Database_old/recon3D_full/Recon3D_301.mat")
recon3D <- recon3D$Recon3D

# get the dataframe for all reactions and the respective genes
rxn_gene_df = as.data.frame(recon3D[[22]])
reaction_ids = unlist(recon3D[[5]])
genes = unlist(recon3D[[10]])
mets = unlist(recon3D[[2]])
rownames(rxn_gene_df) = reaction_ids
entrez = str_split(genes, '\\.', simplify = T)[,1]
symbols = getSYMBOL(entrez, data='org.Hs.eg')
colnames(rxn_gene_df) = symbols

# create a dataframe listing all genes needed in a reaction
row_sums_rnx_genes = rowSums(rxn_gene_df) 
rows = sum(rxn_gene_df) + sum(row_sums_rnx_genes==0) # make a row for every gene in a reaction and every reaction that has no genes
#reaction_to_genes = data.frame(matrix(nrow = rows, ncol = 2))
reaction_to_genes = list()
counter = 1
for(reaction in rownames(rxn_gene_df)){
  print(counter)
  if(row_sums_rnx_genes[reaction] == 0){
    #reaction_to_genes[counter,] = c(reaction, 'Orphan')
    reaction_to_genes[[counter]] = c(reaction, 'Orphan')
    counter = counter + 1
  }else{
    #reaction_to_genes[counter:(counter + row_sums_rnx_genes[reaction] - 1), ] = cbind(rep(reaction, row_sums_rnx_genes[reaction]), colnames(rxn_gene_df)[rxn_gene_df[reaction,] > 0])
    reaction_to_genes[[counter]] = cbind(rep(reaction, row_sums_rnx_genes[reaction]), colnames(rxn_gene_df)[rxn_gene_df[reaction,] > 0])
    #counter = counter + row_sums_rnx_genes[reaction]
    counter = counter + 1
  }
}

reaction_to_genes = ldply(reaction_to_genes, rbind)

# make rows unique 
reaction_to_genes <- reaction_to_genes %>% distinct(.keep_all = TRUE)

#get the stochio matrix
S <- as.matrix(recon3D[[1]])

#get reversible reactions
lbs <- as.data.frame(cbind(recon3D[[6]],recon3D[[7]]))
lbs$direction <- ifelse((recon3D[[7]] + recon3D[[6]]) >= 0,"forward","backward")
lbs$rev = ifelse((lbs$V1 == -1000), 'reversible', 'irreversible')
lbs$rxnid = reaction_ids

colnames(S) = reaction_ids
rownames(S) = mets
S = S[,row_sums_rnx_genes !=0] # exclude orphan reactions

reactions_df = list()
reactions_df_rev = list()
counter = 1
counter2 = 1
for(reaction in colnames(S)){
  print(counter)
  reactands = rownames(S)[S[,reaction] < 0]
  products = rownames(S)[S[,reaction] > 0]
  genes = reaction_to_genes[reaction_to_genes[,1] == reaction,2]
  a = cbind(rep(reactands, length(genes)), sort(rep(genes, length(reactands))))
  b = cbind(sort(rep(genes, length(products))), rep(products, length(genes)))
  c = rbind(a,b)
  if(lbs$rev[counter] == 'reversible'){
    a_rec = c()
    b_rev = c()
    tryCatch({
      a_rev = a[,c(2,1)]
      b_rev = b[,c(2,1)]
    }, error = function(counter){paste('reaction:', counter, 'has either no products or reactands')}
    )
    c_rev = rbind(a_rev, b_rev)
    reactions_df_rev[[counter2]] = as.data.frame(c_rev)
    counter2 = counter2 + 1
  }
  reactions_df[[counter]] = as.data.frame(c)
  counter = counter + 1
}
reactions_df <- ldply(reactions_df, rbind)
reactions_df_rev <- ldply(reactions_df_rev, rbind)
reactions_full = rbind(reactions_df, reactions_df_rev)
write_csv2(reactions_full, '~/Documents/Database_old/recon3D_full/reactions_df_full.csv')

reaction_ends = str_extract(colnames(S), '[:lower:]{1,}')
types = unique(reaction_ends)

#Data cleanup
reactions_df$compartment  = paste0(str_split(reactions_df$V1, '\\[', simplify = T)[,2], str_split(reactions_df$V2, '\\[', simplify = T)[,2])
reactions_df$V1 = str_split(reactions_df$V1, '\\[', simplify = T)[,1]
reactions_df$V2 = str_split(reactions_df$V2, '\\[', simplify = T)[,1]
reactions_df$compartment = str_split(reactions_df$compartment, '\\]', simplify = T)[,1]

HMDB_list = recon3D[[16]]
HMDB_list[sapply(recon3D[[16]], isEmpty)] = 'NA'
HMDB = unlist(HMDB_list)
chebi_list = recon3D[[30]]
chebi_list[sapply(recon3D[[30]], isEmpty)] = 'NA'
chebi = unlist(chebi_list)
chebi = paste0(str_split(chebi, 'CHEBI:', simplify = T)[,1], str_split(chebi, 'CHEBI:', simplify = T)[,2])
chebi[chebi == 'NA'] = NA
chebi[is.na(chebi) == F] = paste0('CHEBI:', chebi[is.na(chebi) == F])
met_mapping = as.data.frame(cbind(mets, HMDB, chebi))

write_csv2(met_mapping, '~/Documents/Database_old/recon3D_full/met_mapping.csv')
write_csv2(reactions_df, '~/Documents/Database_old/recon3D_full/reactions_df_clean.csv')
















