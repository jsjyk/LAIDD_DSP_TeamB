

library(dplyr)

getwd()
dirv = '/home/mentee1/data/filtered'

list.files(dirv)

gene.exp = read.csv(file=paste0(dirv,'/CCLE_TPMLogp1_filtered.csv'))
gene.mut = read.csv(file=paste0(dirv,'/Genetic_feature_filtered.csv'))

drug.res = read.csv(paste0(dirv,'/sanger-dose-response.csv'))

cell.info = read.csv(file=paste0(dirv,'/Cell_info_filtered.csv'))
drug.info = read.csv(file=paste0(dirv,'/screened_compounds_rel_8.4.csv'))

head(gene.exp[,1:4])
head(gene.mut[,1:4])



# get common cell line IDs
cell.ids = intersect(gene.exp$COSMIC.ID, gene.mut$COSMIC.ID)
length(cell.ids)
length(unique(gene.exp$COSMIC.ID))
length(unique(gene.mut$COSMIC.ID))



# GDSC drug response data
head(drug.res)

dim(drug.res)
subset(drug.res, COSMIC_ID %in% cell.ids) %>% dim
subset(drug.res, COSMIC_ID %in% cell.ids)$DRUG_ID %>% unique %>% length()

drug.res = subset(drug.res, COSMIC_ID %in% cell.ids)
drug.res$response = ifelse(drug.res$MAX_CONC > drug.res$IC50_PUBLISHED,1,0)
table(drug.res$response)

hist(log10(drug.res$IC50_PUBLISHED))

getwd()
dir.create('/home/mentee1/hanbi/processed/')
write.csv(drug.res, file='/home/mentee1/hanbi/processed/drug_response_labeled.csv', 
                    row.names=F)

# cell line - tissue info summary
head(cell.info)
cell.info$GDSC.Desc2 %>% table %>% sort # cancer type

# drug MoA summary
drug.info$TARGET_PATHWAY %>% table %>% sort


