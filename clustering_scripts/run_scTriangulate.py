#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import sys
import scanpy as sc
from sctriangulate import *
from sctriangulate.preprocessing import *

input_dir = "/gladstone/mucke/sequencing/Cathrine/cluster_benchmarking/datasets/Srivatsan_2021_mouse_embryo_spatial"

obj = sc.read("/gladstone/mucke/sequencing/Cathrine/cluster_benchmarking/datasets/Srivatsan_2021_mouse_embryo_spatial_adata.h5ad")

sctri = ScTriangulate(dir=input_dir + '/output',adata=obj,query=['CHOIR_1','Seurat_1','Cytocipher_1'])

sctri.compute_metrics(parallel=False)

sctri.serialize('break_point_after_metrics.p')

sctri = ScTriangulate.deserialize('output/break_point_after_metrics.p')

sctri.compute_shapley(parallel=False)

sctri.serialize('break_point_after_shapley.p')

sctri = ScTriangulate.deserialize('output/break_point_after_shapley.p')

sctri.prune_result()

sctri.serialize('break_point_after_prune.p')






