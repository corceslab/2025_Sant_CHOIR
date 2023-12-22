#!/usr/bin/env python
# coding: utf-8

# ---------------------------------------------------------------------------
# Run PanoView
# ---------------------------------------------------------------------------
import random
import time
import os
import pandas as pd
from PanoramicView import scPanoView

def run_PanoView(input_dir,
                 cluster_parameter_index,
                 input_data,
                 Zscore,
                 GeneLow,
                 set_seed):
    
    # Set up ---------------------------
    
    # Working directory
    os.chdir(input_dir + "/intermediate_files/preprocessed")
    
    if "MNN" in input_data:
        Normal = False
        Log2 = False
    else:
        Normal = True
        Log2 = True
    
    # PanoView ---------------------------
    
    random.seed(set_seed)
    
    # Input
    obj = scPanoView.PanoView(filename = input_dir + "/intermediate_files/preprocessed/" + input_data)
    
    # Start time
    start_time = time.time()
    
    try:
        obj.RunSearching(GeneLow = GeneLow, 
                         Zscore = Zscore, 
                         Normal = Normal, 
                         Log2 = Log2)
        print("RunSearching ran.")
    except TypeError:
        # Cluster assignments
        cell_ids = obj.raw_exp.columns.tolist()
        clusters = ['NA'] * len(cell_ids)
        # Time
        time_elapsed = 'NA'
    else:
        try:
            obj.OutputResult()
            print("OutputResult ran.")
        except ValueError:
            # Cluster assignments
            cell_ids = obj.raw_exp.columns.tolist()
            clusters = ['NA'] * len(cell_ids)
            # Time
            time_elapsed = 'NA'
        else:
            # End time
            end_time = time.time()

            # Output ---------------------------

            # Cluster assignments
            cell_ids = obj.cell_membership.iloc[:,0].tolist()
            clusters = obj.cell_membership.iloc[:,2].tolist()
            
            print(obj.cell_membership.iloc[1:5,2].tolist())

            # Time
            time_elapsed = [end_time - start_time]
    
    # List
    output = [cell_ids, clusters, time_elapsed]
    
    return output

