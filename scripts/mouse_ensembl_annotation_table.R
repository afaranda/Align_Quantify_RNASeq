# Setup Workspace
library(dplyr)
library(synapser)
wd <- "/work/abf/MouseEnsembl104"

## Import Table with gene length's estimated from the GTF file
length_table <- read.table(
    "EnsMm_104_Length_GC.tsv", header=F, quote="",
    sep="\t"
)

names(length_table) <- c(
  "gene_id", "tx_id", "eu_length", "eu_gc",
  "pl_length","pl_gc","is_principal"
)
print(paste("Rows in length table:", nrow(length_table)))


## Import Table with gene metadata retreived from Ensembl via the
## Biomart web interface May 24, 2021
meta_table <- read.table(
    "mouse_gene_metadata_EnsMm_104.txt", header=T, quote="",
    sep="\t"
)

names(meta_table) <- c(
  "gene_id", "gene_id_version", "SYMBOL",
  "SEQNAME", "DESCRIPTION", "GENEBIOTYPE"
)

meta_table <- meta_table %>%
    rowwise() %>%
    mutate(
        gene_id_version=gsub(paste0(gene_id,"."),"",gene_id_version)
    )
    
print(paste("Rows in meta table:", nrow(meta_table)))

gene_annotations <- inner_join(
  length_table, meta_table, by="gene_id"
)

print(paste("Rows in annotation table:", nrow(gene_annotations)))

print(head(gene_annotations))
write.csv(
  gene_annotations,
  file="Mouse_Gene_Annotations.csv"
)

## Synapse Login 
synapser::synLogin()

## Specify Project and code synapse IDs
syn_project <- "syn25714110"
syn_code_dir <- "syn25714112"
syn_gene_meta <- "syn25809619"
length_file <- "/work/abf/MouseEnsembl104/Mouse_Gene_Annotations.csv"

# Add this script to the code dir
synapse_push <- File(
  path="/work/abf/MouseEnsembl104/mouse_ensembl_annotation_table.R",
  parent=syn_code_dir
)

synapse_push <- synStore(
  synapse_push
)

# Upload files in the counts folder
length_file <- File(
    path=length_file,
    parent=syn_gene_meta
  )

length_file <- synStore(
    length_file
)

syn_act <- Activity(
  name="upload_gene_annotations",
  description="upload gene length table to Synapse Project"
)

syn_act$executed(synapse_push)
synSetProvenance(
  length_file,
  syn_act
)
