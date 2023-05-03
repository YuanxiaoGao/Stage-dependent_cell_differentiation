#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  4 11:52:22 2018

@author: gao
"""
import numpy as np

num_dupl=100
division_times=np.array([5],dtype=int)
grid_num=41

fp = open("run.sh", "w")
for num_duplicate in range(num_dupl):  
	for n_cluster in division_times:  
		for alpha_cluster in range(0,grid_num):  
			for beta_cluster in range(0,grid_num):  
	
					jobname = "dynamic_heatmap_%d_%d_%d_%d" % (num_duplicate,n_cluster,alpha_cluster,beta_cluster)
					fname = "script/dynamic_heatmap_%d_%d_%d_%d.sh" % (num_duplicate,n_cluster,alpha_cluster,beta_cluster)
					fp.write("sbatch %s\n" % fname)
					bashfp = open(fname, "w")
			
					bashfp.write("#!/bin/sh\n")
					bashfp.write("#SBATCH --time=00-01:30:00\n")
					bashfp.write("#SBATCH --job-name %s\n" % jobname)
					bashfp.write("#SBATCH -o out/%s\n" % jobname)
					bashfp.write("#SBATCH -e err/%s\n" % jobname)
				 
					# record running time				 
					if num_duplicate==num_dupl-1 and alpha_cluster==grid_num-1 and beta_cluster==grid_num-1:
						bashfp.write("#SBATCH --mail-type=ALL\n")
						bashfp.write("#SBATCH --mail-user=gao@evolbio.mpg.de\n")
					 
					 
					bashfp.write("python dynamic_heatmap.py %d %d %d %d\n" % (num_duplicate,n_cluster,alpha_cluster,beta_cluster) )
			
		
		