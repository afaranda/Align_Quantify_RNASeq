library(DegNorm)
args <- commandArgs(trailingOnly=TRUE)
wd <- getwd()
gtf_file <- "/work/abf/MouseEnsembl104/Mus_musculus.GRCm39.104.gtf"
bam_path <- paste(wd, "Alignments", sep="/")

# Use a commandline argument to specify ribosomal or all-genomic alignments
if(args[1] == "geno"){
    bam_files <-paste(
        bam_path, 
        list.files(bam_path, pattern="sorted_alignment\\.bam$"),
        sep="/"
    )
} else if(args[1] == "ribo"){
    bam_files <-paste(
        bam_path, 
        list.files(bam_path, pattern="sorted_rf_alignment\\.bam$"),
        sep="/"
    )

}
print(bam_files)

# Calculate gene level coverage estimates
gene_coverage_batch <-read_coverage_batch(bam_files, gtf_file, cores=8)
fn<-paste0(wd,"/RSeQC_Results/DegNorm_Coverage_Estimates_",args[1],".Rdata")
save(gene_coverage_batch, file=fn)


# Estimate Gene level Degradation Indices
degnorm_results <- degnorm(
    read_coverage = gene_coverage_batch[[1]],
    counts = gene_coverage_batch[[2]],
    down_sampling = 1,
    grid_size = 10,
    loop = 100,
    cores = 8
)

fn<-paste0(wd,"/RSeQC_Results/DegNorm_Analysis_Results_",args[1],".Rdata")
print(fn)
save(degnorm_results, file=fn)
