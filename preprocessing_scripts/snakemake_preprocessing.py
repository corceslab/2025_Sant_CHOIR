#!/usr/bin/env python
# coding: utf-8

# ---------------------------------------------------------------------------
# Snakemake CHOIR benchmarking pipeline (Step 1: Pre-processing)
# ---------------------------------------------------------------------------
import os
import sys
import numpy as np
import pandas as pd
import snakemake.io
import glob

# ---------------------------------------------------------------------------
# Set up
# ---------------------------------------------------------------------------

# Define config variables
in_dir = config['in_dir']
out_dir = config['out_dir']
temp_dir = config['temp_dir']
scripts_dir = config['scripts']
metadata_file = config['metadata']

# Check if metadata file exists
if os.path.isfile(metadata_file):
    metadata = pd.read_csv(metadata_file)
    print("Metadata file exists.")
else:
    error = "Designated metadata file does not exist! exiting..."
    sys.exit(error)
    
# Import values from metadata
name = metadata.loc[metadata['variable'] == 'name']['value'].tolist()[0]
sample_ids_char = metadata.loc[metadata['variable'] == 'sample_ids']['value'].tolist()[0]
sample_ids = sample_ids_char.split(sep = " ") # Create list of sample IDs, uses space as separator
strategy = metadata.loc[metadata['variable'] == 'strategy']['value'].tolist()[0]
tech = metadata.loc[metadata['variable'] == 'tech']['value'].tolist()[0]
cells_nuclei = metadata.loc[metadata['variable'] == 'cells_nuclei']['value'].tolist()[0]
file_type = metadata.loc[metadata['variable'] == 'file_type']['value'].tolist()[0]
cellranger = metadata.loc[metadata['variable'] == 'cellranger']['value'].tolist()[0]
species = metadata.loc[metadata['variable'] == 'species']['value'].tolist()[0]
loaded = metadata.loc[metadata['variable'] == 'loaded']['value'].tolist()[0]
cell_metadata = metadata.loc[metadata['variable'] == 'cell_metadata']['value'].tolist()[0]

print("-----------------------------------")
print("Input directory: ", in_dir)
print("Output directory: ", out_dir)
print("Temporary directory: ", temp_dir)
print("Scripts directory: ", scripts_dir)
print("-----------------------------------")
print("Dataset: ", name)
print("Strategy: ", strategy)
print("Sequencing tech:", tech)
print("Cells or nuclei? ", cells_nuclei)
print("Input file type? ", file_type)
print("Cellranger info? ", cellranger)
print("Species: ", species)
print("# loaded per library: ", loaded)
print("Cell metadata provided? ", cell_metadata)
print("-----------------------------------")

# Define expected input and output files
aligned_read_files = []
preprocessed_output_files = []

if tech == 'RNA' or tech == "RNA_ATAC":
    aligned_read_files.extend(["barcodes.tsv","features.tsv"])
    preprocessed_output_files.append(expand("{input_dir}/quality_checks/RNA_preprocessed_file_list.csv",
                                     input_dir = in_dir))
if tech == 'ATAC' or tech == "RNA_ATAC":
    aligned_read_files.extend(["atac_fragments.tsv.gz", "atac_fragments.tsv.gz.tbi"])
    preprocessed_output_files.append(expand("{input_dir}/quality_checks/ATAC_preprocessed_file_list_ArchR.txt",
                                     input_dir = in_dir))

if file_type == "mtx":
    rna_qc_input = ['barcodes.tsv', 'features.tsv']
    atac_input = ["atac_fragments.tsv.gz", "atac_fragments.tsv.gz.tbi", "peaks.bed"]
else:
    rna_qc_input = aligned_read_files
    atac_input = aligned_read_files


# ---------------------------------------------------------------------------
# Snakemake rules
# ---------------------------------------------------------------------------

localrules: all
    
rule all:
    input:
        cell_metadata = expand("{input_dir}/intermediate_files/preprocessed/cell_metadata.csv",
                   input_dir = in_dir),
    output:
        final="PipelineCompletion.txt"
    threads: 1
    resources:
        mem_qsub = "2G",
        job_time = "00:15:00",
    shell:
        "echo Pipeline created by Cathrine Petersen.; "
        "echo Completed Successfully at `date` | tee {output.final}; "
        
# Read alignment ---------------------------------------

# For bam files, first convert back to fastq using cellranger bamtofastq
if file_type == 'bam':
    rule ConvertBAM:
        input:
            expand("{input_dir}/bam_files/{sample}.bam", input_dir = in_dir, sample = sample_ids),
        output:
            file_list = in_dir + "/quality_checks/{sample}_fastq_file_list.txt"
        threads: 4
        resources:
            mem_qsub = "24G",
            job_time = "08:00:00",
        shell:
            "bash {scripts_dir}/ConvertBAM.sh -i {in_dir} -s {wildcards.sample}"

# For provided/generated fastq files, align reads using cellranger count
if file_type == 'bam' or file_type == "fastq":
    # If bam files were provided, make this rule dependent on previous rule
    AlignReads_input = list()
    if file_type == 'bam':
        AlignReads_input.append(expand("{input_dir}/quality_checks/{sample}_fastq_file_list.txt", 
                                           input_dir = in_dir, sample = sample_ids))
    rule AlignReads:
        input:
            AlignReads_input,
        output:
            expand("{input_dir}/aligned_reads/{sample}/{files}",
                   input_dir = in_dir,
                   files = aligned_read_files,
                   allow_missing = True),
            expand("{input_dir}/quality_checks/{sample}_web_summary.html",
                   input_dir = in_dir,
                   allow_missing = True)
        threads: 8
        resources:
            mem_qsub = "8G",
            job_time = "24:00:00",
        shell:
            "bash {scripts_dir}/AlignReads.sh -i {in_dir} -a {species} -t {tech} -c {cells_nuclei} -s {wildcards.sample}"

# Quality control, normalization, & integration ---------------------------------------

# RNA
if tech == 'RNA' or tech == 'RNA_ATAC':
    rule QualityControlRNA:
        input:
            expand("{input_dir}/aligned_reads/{sample}/{files}",
                   input_dir = in_dir,
                   files = rna_qc_input,
                   allow_missing = True),
        output:
            expand("{input_dir}/intermediate_files/qc/{sample}/qc_cell_metadata.rds", 
                   input_dir = in_dir, 
                   allow_missing = True),
            expand("{input_dir}/intermediate_files/qc/{sample}/qc_matrix.rds", 
                   input_dir = in_dir, 
                   allow_missing = True),
        threads: 4
        resources:
            mem_qsub = "16G",
            job_time = "12:00:00",
        shell: 
            "Rscript {scripts_dir}/QualityControlRNA.R {in_dir} {wildcards.sample} {strategy} {tech} {cells_nuclei} {cellranger} {species} {loaded} {cell_metadata}"
            
    rule FilterIntegrateNormRNA:
        input:
            expand("{input_dir}/intermediate_files/qc/{sample}/qc_cell_metadata.rds", 
                   input_dir = in_dir, 
                   sample = sample_ids),
            expand("{input_dir}/intermediate_files/qc/{sample}/qc_matrix.rds", 
                   input_dir = in_dir, 
                   sample = sample_ids),
        output:
            preprocessed_list = in_dir + "/quality_checks/RNA_preprocessed_file_list.csv",
        threads: 1
        resources:
            mem_qsub = "96G",
            job_time = "48:00:00",
        shell:
            "Rscript {scripts_dir}/FilterIntegrateNormRNA.R {in_dir} {strategy} {sample_ids_char}"

                      
# ATAC
if tech == 'ATAC' or tech == 'RNA_ATAC':
    rule ArrowFilesATAC:
        input:
            expand("{input_dir}/aligned_reads/{sample}/{files}",
                   input_dir = in_dir,
                   files = atac_input,
                   sample = sample_ids),
        output:
            expand("{input_dir}/intermediate_files/arrow_files/arrow_files.rds",
                   input_dir = in_dir),
        threads: 1
        resources:
            mem_qsub = "96G",
            job_time = "24:00:00",
        shell:
            "Rscript {scripts_dir}/ArrowFilesATAC.R {in_dir} {species} {sample_ids_char}"
            
    rule PreProcessArchR:
        input:
            expand("{input_dir}/intermediate_files/arrow_files/arrow_files.rds",
                   input_dir = in_dir),
        output:
            preprocessed_list = in_dir + "/quality_checks/ATAC_preprocessed_file_list_ArchR.txt",
        threads: 1
        resources:
            mem_qsub = "24G",
            job_time = "12:00:00",
        shell:
            "Rscript {scripts_dir}/PreProcessArchR.R {in_dir} {species} {cell_metadata} {strategy} {sample_ids_char}"

# Compile cell metadata file ---------------------------------------
rule CellMetadata:
    input:
        preprocessed_output_files
    output:
        cell_metadata = expand("{input_dir}/intermediate_files/preprocessed/cell_metadata.csv",
                   input_dir = in_dir),
    threads: 1
    resources:
        mem_qsub = "4G",
        job_time = "02:00:00",
    shell:
        "Rscript {scripts_dir}/CellMetadata.R {in_dir} {tech}"