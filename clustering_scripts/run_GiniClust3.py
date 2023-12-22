#!/usr/bin/env python
# coding: utf-8

# ---------------------------------------------------------------------------
# Run GiniClust3
# ---------------------------------------------------------------------------
import random
import time
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
import giniclust3 as gc

def run_GiniClust3(input_dir,
                   cluster_parameter_index,
                   input_data, 
                   selection, 
                   p_value, 
                   min_gini_value, 
                   resolution, 
                   method, 
                   set_seed):
 
    # GiniClust3 ---------------------------
    random.seed(set_seed)

    # Input
    obj = sc.read_csv(input_dir + "/intermediate_files/preprocessed/" + input_data + ".csv", first_column_names = True)
    obj = obj.transpose()
    
    if "raw" in input_data:
        sc.pp.normalize_per_cell(obj, counts_per_cell_after = 1e4)
        
    # Start time
    start_time = time.time()

    # Cluster
    
    # Gini index
    try:
        gc.gini.calGini(obj,
                        selection = selection,
                        p_value = p_value,
                        min_gini_value = min_gini_value)
    except (SystemExit, KeyError):
        # Cluster assignments
        cell_ids = obj.obs.index.tolist()
        clusters = ['NA'] * len(cell_ids)
        # Time
        time_elapsed = 'NA'
    else:
        try:
            obj_gini = gc.gini.clusterGini(obj, 
                                           neighbors = 15, 
                                           resolution = resolution,
                                           method = method)
        except (ValueError, KeyError):
            # Cluster assignments
            cell_ids = obj.obs.index.tolist()
            clusters = ['NA'] * len(cell_ids)
            # Time
            time_elapsed = 'NA'
        else:
            try:
                # Fano factor
                gc.fano.calFano(obj)
            except KeyError:
                # Cluster assignments
                cell_ids = obj.obs.index.tolist()
                clusters = ['NA'] * len(cell_ids)
                # Time
                time_elapsed = 'NA'
            else:
                obj_fano=gc.fano.clusterFano(obj,
                                           neighbors = 15,
                                           resolution = resolution,
                                           method = method)
                # Consensus
                obj_consensus = {}
                obj_consensus['giniCluster'] = np.array(obj.obs['rare'].values.tolist())
                obj_consensus['fanoCluster'] = np.array(obj.obs['fano'].values.tolist())
                gc.consensus.generateMtilde(obj_consensus)
                try:
                    gc.consensus.clusterMtilde(obj_consensus)
                except (ValueError, KeyError):
                    # Cluster assignments
                    cell_ids = obj.obs.index.tolist()
                    clusters = ['NA'] * len(cell_ids)

                    # Time
                    time_elapsed = 'NA'
                else:
                    # End time
                    end_time = time.time()

                    # Output ---------------------------

                    # Cluster assignments
                    cell_ids = obj.obs.index.tolist()
                    clusters = obj_consensus['finalCluster']

                    # Time
                    time_elapsed = [end_time - start_time]
    

    # List
    output = [cell_ids, clusters, time_elapsed]
    
    return output

