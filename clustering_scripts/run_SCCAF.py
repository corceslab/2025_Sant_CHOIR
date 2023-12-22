#!/usr/bin/env python
# coding: utf-8

# ---------------------------------------------------------------------------
# Run SCCAF
# ---------------------------------------------------------------------------
import random
import time
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
from SCCAF import *
import scanpy.external as sce

def run_SCCAF(input_dir,
              cluster_parameter_index,
              input_data,
              resolution, 
              min_acc,
              set_seed,
              batch_correction_method):

    # SCCAF ---------------------------
    
    random.seed(set_seed)
    
    # Input
    obj = sc.read_csv(input_dir + "/intermediate_files/preprocessed/" + input_data + ".csv", first_column_names = True)
    obj = obj.transpose()
    
    if batch_correction_method == "Harmony":
        # Add cell metadata
        metadata = pd.read_csv(input_dir + "/intermediate_files/preprocessed/cell_metadata.csv", index_col='CellID')
        obj.obs["Batch"] = metadata["Batch"].tolist()
        
    # Start time
    start_time = time.time()
    
    # Prep
    sc.pp.highly_variable_genes(obj)
    if "log" in input_data:
        sc.pp.scale(obj)
    sc.tl.pca(obj, svd_solver='arpack')
    
    if batch_correction_method == "Harmony":
        # Harmony batch correction
        sce.pp.harmony_integrate(obj, 'Batch')
        sc.pp.neighbors(obj, use_rep='X_pca_harmony')
    else:
        sc.pp.neighbors(obj)
    
    # Cluster
    sc.tl.leiden(obj, resolution = resolution, key_added = 'leiden_r1')

    obj.obs['L1_Round0'] = obj.obs['leiden_r1']
    
    try:
        SCCAF_optimize_all(ad = obj,
                           prefix = 'L1',
                           min_acc = min_acc,
                           use = 'pca', 
                           plot = False)
    except KeyError:
        try:
            SCCAF_optimize_all(ad = obj,
                           prefix = 'L1',
                           min_acc = min_acc,
                           use = 'pca', 
                           plot = False)
        except KeyError:
            # Cluster assignments
            cell_ids = obj.obs.index.tolist()
            clusters = ['NA'] * len(cell_ids)

            # Time
            time_elapsed = 'NA'
        else:
            # End time
            end_time = time.time()

            # Output ---------------------------

            # Cluster assignments dataframe
            cell_ids = obj.obs.index.tolist()
            clusters = obj.obs.L1_result.tolist()

            # Time
            time_elapsed = [end_time - start_time]
    else:
        # End time
        end_time = time.time()

        # Output ---------------------------

        # Cluster assignments dataframe
        cell_ids = obj.obs.index.tolist()
        clusters = obj.obs.L1_result.tolist()

        # Time
        time_elapsed = [end_time - start_time]
    
    # List
    output = [cell_ids, clusters, time_elapsed]
    
    # Write csv with original clusters
    cluster_labels = pd.DataFrame(obj.obs['leiden_r1'].values, index=obj.obs.index.tolist(), columns=['Cluster'])
    cluster_labels.to_csv(input_dir + "/intermediate_files/clusters/SCCAF_original_clusters/" + 
                          "cluster_labels_" + cluster_parameter_index + ".csv")
    
    return output

