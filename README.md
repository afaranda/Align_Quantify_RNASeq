# Align_Quantify_RNASeq

This repository contains a collection of shell scripts intended to
automate the processing of RNA Sequencing data.

### Basic Workflow
1. Clone This Reopsitory
2. Build A Genome Index
  * Download any Genomic data not automatically fetched by build script
  * Run appropriate Build Script
  * Download Ribosomal sequence and build index using hisat2-build
  * Update PE_ribsosomal_content_analysis.sh with appropriate length
  * Run appropriate \*_make_rseqc_bedfile.sh script (update filepaths)
  * Generate additional RSeQC bedfiles for targeted analysis eg. crystallins.
  * Estimate Gene Level length and GC content using Extract_Length_GC.sh
3. Set Up Analysis Script
  * Edit header for appropriate file paths / index labels
4. Run Main Analysis
  * Push main analysis script and cross fingers
5. Run extended analyses

### Scripts

### Dependencies
Package | URL
------- | ------
Hisat2 | http://daehwankimlab.github.io/hisat2/
bedtools | https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html
bedops | https://bedops.readthedocs.io/en/latest/
htseq-count | https://htseq.readthedocs.io/en/master/count.html
FastQC | https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
Trim Galore! | https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
RSeQC | http://rseqc.sourceforge.net/
MultiQC | https://multiqc.info/

### To Do
  * Fix call to dexseq_counts.py to point to packaged installation
  * setup per-sample workflows to define environment variables independently
  * dockerize pipeline

### Changelog

**April 19, 2021 -- Added DEXSeq Analysis**
  * Added DEXSeq analysis to paired and single end pipelines
  * Added generic R script for running differential Exon Analysis
  * Copied dexseq_counts.py from dexseq package

**March 1, 2021 -- updated Paired End Pipeline**
  * revised approach to ribosomal filtering
  * added kallisto abundance estimates
  * added DEGNorm analysis
