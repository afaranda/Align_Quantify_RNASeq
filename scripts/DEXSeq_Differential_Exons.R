##############################################################################
#                                                                            #
#  File: DEXSeq_Differential_Exons.R                                         #
#  Author: Adam Faranda                                                      #
#  Created: April 5, 2021                                                    #
#  Purpose:  Use the DEXSeq method to identify genes with significant        #
#            a significant difference in exon composition.                   #
#                                                                            #
##############################################################################
print(date())
library(DEXSeq)
library(dplyr)
options(echo=T)

# Enter Working Directory and Load Raw Data
setwd('/work/abf/23Mar2021_DNA_Link_LTS_Ensembl101')
wd<-getwd()
results<-paste(wd,'DEXSeq_Output',sep='/')
data_dir<-paste(wd,'Counts',sep='/')


## Build sample table
## The "group" column must be a factor that assigns samples to experimental
## conditions. The first level of this factor will be the refference in the
## design matrix used to model exon abundance.
sample_files <- list.files(data_dir, pattern="rf_DEXSeq_Counts")
samples<-data.frame(
  countFile=paste(
        data_dir, "/",
        list.files(data_dir, pattern="rf_DEXSeq_Counts"),
        sep=""
  ),
  group=factor(x=c(), levels=c()) # EDIT THIS FACTOR TO MATCH STUDY
) %>% group_by(group) %>%
  mutate(sample = paste(group, row_number(), sep="_")) %>%
  select(sample, group, countFile) %>% as.data.frame()

## Build DEXSeq Analysis object from Counts
dxd <- DEXSeqDataSetFromHTSeq(
  countfiles = samples$countFile,
  sampleData = samples,
  design = ~ sample + exon + group:exon,
  flattenedfile = "/work/abf/MouseEnsembl101/Ens_Mm101_DEXSeq_Annot.gff"
)

#gid<-unique(geneIDs(dxd))[1:10]
#dxd<-dxd[geneIDs(dxd) %in% gid,]

dxd<-estimateSizeFactors(dxd)
dxd<-estimateDispersions(dxd)
# plotDispEsts(dxd)
dxd <- testForDEU(dxd)
dxd <- estimateExonFoldChanges( dxd, fitExpToVar="group")
dxr <- DEXSeqResults(dxd)

save(dxd, dxr, file="Counts/DEXSeq_Results_Object.Rdata")
print(date())

