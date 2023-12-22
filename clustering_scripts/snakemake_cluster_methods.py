#!/usr/bin/env python
# coding: utf-8

# ---------------------------------------------------------------------------
# Snakemake CHOIR benchmarking pipeline (Step 2: Clustering)
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
parameter_file = config['parameters']

# Check if metadata file exists
if os.path.isfile(metadata_file):
    metadata = pd.read_csv(metadata_file)
    print("Metadata file exists.")
else:
    error = "Designated metadata file does not exist! exiting..."
    sys.exit(error)
    
# Import values from metadata
name = metadata.loc[metadata['variable'] == 'name']['value'].tolist()[0]
data_type = metadata.loc[metadata['variable'] == 'type']['value'].tolist()[0]
strategy = metadata.loc[metadata['variable'] == 'strategy']['value'].tolist()[0]
tech = metadata.loc[metadata['variable'] == 'tech']['value'].tolist()[0]
cells_nuclei = metadata.loc[metadata['variable'] == 'cells_nuclei']['value'].tolist()[0]
species = metadata.loc[metadata['variable'] == 'species']['value'].tolist()[0]
cell_metadata = metadata.loc[metadata['variable'] == 'cell_metadata']['value'].tolist()[0]

print("-----------------------------------")
print("Input directory: ", in_dir)
print("Output directory: ", out_dir)
print("Temporary directory: ", temp_dir)
print("Scripts directory: ", scripts_dir)
print("-----------------------------------")
print("Dataset: ", name)
print("Data type: ", data_type)
print("Strategy: ", strategy)
print("Sequencing tech:", tech)
print("Cells or nuclei? ", cells_nuclei)
print("Species: ", species)
print("Cell metadata provided? ", cell_metadata)
print("-----------------------------------")

# Check if parameter list file exists
if os.path.isfile(parameter_file):
    parameter_key = re.search(r'/([^/]+).csv$', parameter_file).group(1)
    cluster_parameters = pd.read_csv(parameter_file)
    print("Parameter list file exists.")
    parameter_indices = np.array(cluster_parameters.index) + 1
else:
    error = "Designated parameter list file does not exist! exiting..."
    sys.exit(error)

print("Running ", np.max(parameter_indices), " cluster method/parameter combinations.")

# Submission parameters depend on data type

if data_type == "small":
    mem_A = "24G"
    time_A = "96:00:00"
    mem_G = "8G"
    time_G = "96:00:00"
if data_type == "large":
    mem_A = "128G"
    time_A = "96:00:00"
    mem_C = "16G"
    time_C = "96:00:00"

# Indices of cluster/parameter combinations to be run with 1 vs. 16 threads
category = np.array(cluster_parameters.loc[:,'snakemake_category'])
parameter_indices_A = parameter_indices[np.where(category == "A")]
parameter_indices_B = parameter_indices[np.where(category == "B")]

# Define files expected after running all cluster method/parameter combinations
compile_clusters_input = []
# Add cluster results
compile_clusters_input.append(expand("{input_dir}/intermediate_files/clusters/{param_key}/clusters_parameter_set_{indices}_A.csv",
                                     input_dir = in_dir,
                                     param_key = parameter_key,
                                     indices = parameter_indices_A))
compile_clusters_input.append(expand("{input_dir}/intermediate_files/clusters/{param_key}/clusters_parameter_set_{indices}_B.csv",
                                     input_dir = in_dir,
                                     param_key = parameter_key,
                                     indices = parameter_indices_B))

# Add timekeeping records
compile_clusters_input.append(expand("{input_dir}/intermediate_files/clusters/{param_key}/time_parameter_set_{indices}_A.csv",
                                     input_dir = in_dir,
                                     param_key = parameter_key,
                                     indices = parameter_indices_A))
compile_clusters_input.append(expand("{input_dir}/intermediate_files/clusters/{param_key}/time_parameter_set_{indices}_B.csv",
                                     input_dir = in_dir,
                                     param_key = parameter_key,
                                     indices = parameter_indices_B))

# ---------------------------------------------------------------------------
# Snakemake rules
# ---------------------------------------------------------------------------

localrules: all
    
rule all:
    input:
        compiled_cluster_results = in_dir + "/intermediate_files/clusters/compiled_clusters_" + parameter_key + ".csv",
        compiled_time_records = in_dir + "/intermediate_files/clusters/compiled_time_" + parameter_key + ".csv",
    output:
        final = "PipelineCompletion.txt"
    threads: 1
    resources:
        mem_qsub = "2G",
        job_time = "00:15:00",
    shell:
        "echo Pipeline created by Cathrine Petersen.; "
        "echo Completed Successfully at `date` | tee {output.final}; "
        
# Clustering

rule A_Cluster:
    input:
    output:
        cluster_ids = in_dir + "/intermediate_files/clusters/" + parameter_key + "/clusters_parameter_set_{parameter_index}_A.csv",
        time = in_dir + "/intermediate_files/clusters/" + parameter_key + "/time_parameter_set_{parameter_index}_A.csv",
    threads: 1
    resources:
        mem_qsub = mem_A,
        job_time = time_A,
    shell:
        "Rscript {scripts_dir}/RunClusterMethod.R {scripts_dir} {in_dir} {temp_dir} {parameter_file} {parameter_key} {wildcards.parameter_index}; "
        '[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID" > {in_dir}/intermediate_files/clusters/{parameter_key}/qstat_parameter_set_{wildcards.parameter_index}_A.txt; '
        'grep "usage" {in_dir}/intermediate_files/clusters/{parameter_key}/qstat_parameter_set_{wildcards.parameter_index}_A.txt > {in_dir}/intermediate_files/clusters/{parameter_key}/usage_parameter_set_{wildcards.parameter_index}_A.txt; '
        'rm {in_dir}/intermediate_files/clusters/{parameter_key}/qstat_parameter_set_{wildcards.parameter_index}_A.txt'
        
rule B_Cluster:
    input:
    output:
        cluster_ids = in_dir + "/intermediate_files/clusters/" + parameter_key + "/clusters_parameter_set_{parameter_index}_B.csv",
        time = in_dir + "/intermediate_files/clusters/" + parameter_key + "/time_parameter_set_{parameter_index}_B.csv",
    threads: 16
    resources:
        mem_qsub = mem_B,
        job_time = time_B,
    shell:
        "Rscript {scripts_dir}/RunClusterMethod.R {scripts_dir} {in_dir} {temp_dir} {parameter_file} {parameter_key} {wildcards.parameter_index}; "
        '[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID" > {in_dir}/intermediate_files/clusters/{parameter_key}/qstat_parameter_set_{wildcards.parameter_index}_B.txt; '
        'grep "usage" {in_dir}/intermediate_files/clusters/{parameter_key}/qstat_parameter_set_{wildcards.parameter_index}_B.txt > {in_dir}/intermediate_files/clusters/{parameter_key}/usage_parameter_set_{wildcards.parameter_index}_B.txt; '
        'rm {in_dir}/intermediate_files/clusters/{parameter_key}/qstat_parameter_set_{wildcards.parameter_index}_B.txt'

# Compile clustering results

rule CompileClusters:
    input:
        compile_clusters_input
    output:
        compiled_cluster_results = in_dir + "/intermediate_files/clusters/compiled_clusters_" + parameter_key + ".csv",
        compiled_time_records = in_dir + "/intermediate_files/clusters/compiled_time_" + parameter_key + ".csv",
    threads: 1
    resources:
        mem_qsub = "12G",
        job_time = "02:00:00",
    shell:
        "Rscript {scripts_dir}/CompileClusters.R {in_dir} {parameter_key}"

