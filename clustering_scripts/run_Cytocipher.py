#!/usr/bin/env python
# coding: utf-8

# ---------------------------------------------------------------------------
# Run Cytocipher
# ---------------------------------------------------------------------------
import random
import time
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
import cytocipher as cc
import scanpy.external as sce

def run_Cytocipher(input_dir,
                   cluster_parameter_index,
                   input_data,
                   alpha,
                   resolution, 
                   n_markers,
                   use_p_adj,
                   set_seed,
                   n_cores,
                   batch_correction_method):

    # Cytocipher ---------------------------
    
    random.seed(set_seed)
    
    if use_p_adj == "TRUE":
        use_p_adj = True
    else:
        use_p_adj = False
    
    # Input
    obj = sc.read_csv(input_dir + "/intermediate_files/preprocessed/" + input_data + ".csv", first_column_names = True)
    obj = obj.transpose()
    
    if batch_correction_method == "Harmony":
        # Add cell metadata
        metadata = pd.read_csv(input_dir + "/intermediate_files/preprocessed/cell_metadata.csv", index_col='CellID')
        obj.obs["Batch"] = metadata["Batch"].tolist()
        print(obj.obs["Batch"])
    
    # Start time
    start_time = time.time()
    
    # Prep
    sc.pp.highly_variable_genes(obj)
    sc.tl.pca(obj, svd_solver='arpack')
    
    if batch_correction_method == "Harmony":
        # Harmony batch correction
        sce.pp.harmony_integrate(obj, 'Batch')
        print("Harmony done.")
        sc.pp.neighbors(obj, use_rep='X_pca_harmony')
        print("Neighbors done.")
    else:
        sc.pp.neighbors(obj)
    
    # Cluster
    sc.tl.leiden(obj, resolution = resolution, key_added = 'leiden')
    
    print("Leiden done.")
    
    try:
        cc.tl.get_markers(obj, 'leiden', n_top = n_markers)
    except ValueError:
        # Cluster assignments
        cell_ids = obj.obs.index.tolist()
        clusters = ['NA'] * len(cell_ids)
        # Time
        time_elapsed = 'NA'
    else:
        try:
            cc.tl.code_enrich(obj, 'leiden', n_cpus = n_cores)
        except Exception:
            # Cluster assignments
            cell_ids = obj.obs.index.tolist()
            clusters = ['NA'] * len(cell_ids)
            # Time
            time_elapsed = 'NA'
        else:
            try:
                cc.tl.merge_clusters(obj, 'leiden', n_cpus = n_cores, p_cut = alpha, n_top_genes = n_markers, p_adjust = use_p_adj)

            except (Exception, ValueError):
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
                clusters = obj.obs.leiden_merged.tolist()

                # Time
                time_elapsed = [end_time - start_time]
    
    # List
    output = [cell_ids, clusters, time_elapsed]
    
    # Write csv with original clusters
    cluster_labels = pd.DataFrame(obj.obs['leiden'].values, index=obj.obs.index.tolist(), columns=['Cluster'])
    cluster_labels.to_csv(input_dir + "/intermediate_files/clusters/Cytocipher_original_clusters/" + 
                          "cluster_labels_" + cluster_parameter_index + ".csv")
    
    return output
    

