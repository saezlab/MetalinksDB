##################### general stuff ############################################
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
lbs_cut = lbs[row_sums_rnx_genes !=0,]
###################### create met/gene links for all reactions #################
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
  if(lbs_cut$rev[counter] == 'reversible'){
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
reactions_full = reactions_full %>% distinct()
#write_csv2(reactions_full, '~/Documents/Database_old/recon3D_full/reactions_df_full.csv')

#################### create met/gene links for transport reactions #############
#get all transport reactions
reaction_names = as.data.frame(unlist(recon3D[[24]])[row_sums_rnx_genes !=0])
subsystem = as.data.frame(unlist(recon3D[[29]])[row_sums_rnx_genes !=0])
colnames(subsystem) = 'ss'
colnames(reaction_names) = 'desc'
reaction_ends = str_extract(colnames(S), '[:lower:]{1,}')
reaction_ends[is.na(reaction_ends)] = 'no_ending'
types = unique(reaction_ends)
transport = colnames(S)[reaction_ends %in% c('t', 'tmi', 'te', )]
transport1 = rowSums(sapply(c('influx', 'efflux', 'Efflux', 'transport', 'Transport', 'uptake', 'Uptake', 'Symport', 'symport', 'Antiport', 
                              'antiport', 'uniport', 'Uniport', 'excretion', 'release'), grepl, x = reaction_names$desc, fixed = T))
transport2 = sapply(c('t'), grepl, x = reaction_ends, fixed = T)
transport3 = rowSums(sapply(c('transport', 'Transport'), grepl, x = subsystem$ss, fixed = T))

#get outtake reactions
ST = S[,transport3>0] 
lbsT = lbs_cut[transport3>0,]

treactions = list()
reactands = list()
products = list()
counter = 1
for(reaction in colnames(ST)){
  reactands[[counter]] = rownames(ST)[ST[,reaction] < 0]
  products[[counter]] = rownames(ST)[ST[,reaction] > 0]
  genes = reaction_to_genes[reaction_to_genes[,1] == reaction,2]
  counter = counter +1
}
reactands = ldply(reactands, rbind)
products <- ldply(products, rbind)
in_out = cbind(reactands,products)

comp = list()
met = list()
for(i in 1:9){
  a = str_split(in_out[,i], '\\[', simplify = T)[,2]
  b = str_split(a, '\\]', simplify = T)[,1]
  c = str_split(in_out[,i], '\\[', simplify = T)[,1]
  comp[[i]] = b
  met[[i]]  = c
}
comp = t(ldply(comp, rbind))
met = t(ldply(met, rbind))
in_out_df = as.data.frame(cbind(met, comp))
rownames(in_out_df) = colnames(ST)

in_out_df$rev = lbsT$rev
in_out_df$label = 'intra'
in_out_df$transport_in = 'no transport_in'
in_out_df$transport_out = 'no transport_out'
in_out_df$score = 0
in_out_m = as.matrix(in_out_df)
in_mets = list()
out_mets = list()
counter = 1
for(rxn in rownames(in_out_df)){
  line= in_out_m[rxn,]
  abundances = as.data.frame(table(line[1:9]))
  candidates = abundances$Var1[abundances$Freq ==2]
  in_score = 0
  out_score = 0
  if(length(candidates) > 0){
    in_out_df[rxn,'label'] = 'outer_membrane_transport'
    pos_met = which(line %in% candidates)
    pos_comp = pos_met + 9
    pos_in = pos_comp[1:length(pos_met)/2]
    pos_out = pos_comp[length(pos_met)+1/2:length(pos_met)]
    ins = line[pos_in]
    outs= line[pos_out]
    tin = list()
    tout = list()
    for (j in 1:length(candidates)){
      if(ins[j] == 'e' & outs[j] == 'c'){
        in_score = in_score +1
        tin[[j]] = candidates[j]
      }else if(ins[j] == 'c' & outs[j] == 'e'){
        out_score = out_score +1
        tout[[j]] = candidates[j]
      }
    }
    in_out_df[rxn, 'transport_in'] =  paste(unlist(tin), collapse = ', ')
    in_out_df[rxn, 'transport_out'] =  paste(unlist(tout), collapse = ', ')
    if(lbsT$rev[counter] == 'reversible'){
      fused = as.vector(unlist(c(tin, tout)))
      tin = fused
      tout = fused
    }
    in_mets[[counter]] = as.vector(unlist(tin))
    out_mets[[counter]] = as.vector(unlist(tout))
    counter = counter + 1
  }
  
  in_out_df[rxn,'score'] = in_score + out_score
}
in_mets = unlist(in_mets)
out_mets = unlist(out_mets)
in_met_uni = unique(in_mets)
out_met_uni = unique(out_mets)

rev_df = in_out_df[in_out_df$rev == 'reversible',]
#rownames(rev_df) = paste0(rownames(rev_df), '_rev')
rev_df = rev_df[,c(5:9,1:4,14:18,10:13,19:20, 22,21,23)]
colnames(rev_df) = colnames(in_out_df)


#################### get intake df #############################################
intake = rbind(in_out_df[in_out_df$transport_in != '',], rev_df[rev_df$transport_in != '',] )
intake = intake[intake$label != 'intra',]
  
STI = ST[,rownames(intake)]


reactions_df = list()
counter = 1
for(reaction in colnames(STI)){
  print(counter)
  if(lbsT$rev[lbsT$rxnid == reaction] == 'irreversible'){
    reactands = rownames(STI)[STI[,reaction] < 0]
    genes = reaction_to_genes[reaction_to_genes[,1] == reaction,2]
    a = cbind(rep(reactands, length(genes)), sort(rep(genes, length(reactands))))
    reactions_df[[counter]] = a
    counter = counter + 1
  }else if(lbsT$rev[lbsT$rxnid == reaction] == 'reversible'){
    products = rownames(S)[S[,reaction] > 0]
    genes = reaction_to_genes[reaction_to_genes[,1] == reaction,2]
    b = cbind(rep(products, length(genes)), sort(rep(genes, length(products))))
    reactions_df[[counter]] = b
    counter = counter + 1
  }
}

reactions_df <- ldply(reactions_df, rbind) # fuse small dfs
reactions_df <- reactions_df %>% distinct(.keep_all = TRUE) #erase dublicated columns
clean_mets = str_split(reactions_df$`1`, '\\[', simplify = T)[,1]
reactions_df = reactions_df[clean_mets %in% in_met_uni, ] # clean out antiporter molecules
#write_csv2(reactions_df, '~/Documents/GitHub/metabolicCCC/Ressource/reactions_tin.csv')


#################### get outtake df #############################################
outtake = rbind(in_out_df[in_out_df$transport_out != '',], rev_df[rev_df$transport_out != '',] )
outtake = outtake[outtake$label != 'intra',]

STO = ST[,rownames(outtake)]

reactions_df = list()
counter = 1
for(reaction in colnames(STO)){
  print(counter)
  if(lbsT$rev[lbsT$rxnid == reaction] == 'irreversible'){
    reactands = rownames(STO)[STO[,reaction] < 0]
    genes = reaction_to_genes[reaction_to_genes[,1] == reaction,2]
    a = cbind(rep(reactands, length(genes)), sort(rep(genes, length(reactands))))
    reactions_df[[counter]] = a
    counter = counter + 1
  }else if(lbsT$rev[lbsT$rxnid == reaction] == 'reversible'){
    products = rownames(S)[S[,reaction] > 0]
    genes = reaction_to_genes[reaction_to_genes[,1] == reaction,2]
    b = cbind(rep(products, length(genes)), sort(rep(genes, length(products))))
    reactions_df[[counter]] = b
    counter = counter + 1
  }
}

reactions_df <- ldply(reactions_df, rbind) # fuse small dfs
reactions_df <- reactions_df %>% distinct(.keep_all = TRUE) #erase dublicated columns
clean_mets = str_split(reactions_df$`1`, '\\[', simplify = T)[,1]
reactions_df = reactions_df[clean_mets %in% out_met_uni, ] # clean out antiporter molecules
#write_csv2(reactions_df, '~/Documents/GitHub/metabolicCCC/Ressource/reactions_tout.csv')

###################### create met/gene links for prod/deg reactions #################
SPD = S[,transport3 == 0]

reactions_df = list()
reactions_df_rev = list()
counter = 1
counter2 = 1
for(reaction in colnames(SPD)){
  print(counter)
  reactands = rownames(SPD)[SPD[,reaction] < 0]
  products = rownames(SPD)[SPD[,reaction] > 0]
  genes = reaction_to_genes[reaction_to_genes[,1] == reaction,2]
  a = cbind(rep(reactands, length(genes)), sort(rep(genes, length(reactands))))
  b = cbind(sort(rep(genes, length(products))), rep(products, length(genes)))
  c = rbind(a,b)
  if(lbs_cut$rev[counter] == 'reversible'){
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
reactions_full = reactions_full %>% distinct()

# filter out SLC and ABC proteins manually

transport_slc1 = sapply(c('SLC', 'ABC'), grepl, x = reactions_full$V1, fixed = T)
transport_slc2 = sapply(c('SLC', 'ABC'), grepl, x = reactions_full$V2, fixed = T)
transport_slc = rowSums(cbind(transport_slc1, transport_slc2))
reactions_full = reactions_full[transport_slc==0,]

#write_csv2(reactions_full, '~/Documents/Database_old/recon3D_full/reactions_df_pd.csv')
#################### create metabolite mapping #################################
#Data cleanup
reactions_df$compartment  = paste0(str_split(reactions_df$V1, '\\[', simplify = T)[,2], str_split(reactions_df$V2, '\\[', simplify = T)[,2])
reactions_df$V1 = str_split(reactions_df$V1, '\\[', simplify = T)[,1]
reactions_df$V2 = str_split(reactions_df$V2, '\\[', simplify = T)[,1]
reactions_df$compartment = str_split(reactions_df$compartment, '\\]', simplify = T)[,1]

mets = unlist(recon3D[[2]])
HMDB_list = recon3D[[16]]
HMDB_list[sapply(recon3D[[16]], isEmpty)] = 'NA'
HMDB = unlist(HMDB_list)
for(i in 1:length(HMDB)){
  if(nchar(HMDB[i]) == 9){
    a = str_split(HMDB[i], 'MDB', simplify = T)[2]
    HMDB[i] = paste0('HMDB00', a)
  } else if(nchar(HMDB[i]) == 10){
    a = str_split(HMDB[i], 'MDB', simplify = T)[2]
    HMDB[i] = paste0('HMDB0', a)
  }
}




chebi_list = recon3D[[30]]
chebi_list[sapply(recon3D[[30]], isEmpty)] = 'NA'
chebi = unlist(chebi_list)
chebi = paste0(str_split(chebi, 'CHEBI:', simplify = T)[,1], str_split(chebi, 'CHEBI:', simplify = T)[,2])
chebi[chebi == 'NA'] = NA
chebi[is.na(chebi) == F] = paste0('CHEBI:', chebi[is.na(chebi) == F])
met_mapping = as.data.frame(cbind(mets, HMDB, chebi))

write_csv2(met_mapping, '~/Documents/Database_old/recon3D_full/met_mapping2.csv')
write_csv2(reactions_df, '~/Documents/Database_old/recon3D_full/reactions_df_clean.csv')
















