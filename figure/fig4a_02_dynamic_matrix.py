#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 14:40:40 2018
@author: gao
"""
import sys
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
#-------------------------------------------------------------------------------
# set display width
np.get_printoptions()['linewidth']
np.set_printoptions(linewidth=2000)
#-------------------------------------------------------------------------------
'''Parameters'''
alpha=[-1,2]
beta=[-1,2]
n=5
M,SAM=[1000,100]
num_dup=100
#-------------------------------------------------------------------------------
'''read data '''
grid_num=41
alpha_expo_range=np.array([alpha[0],alpha[1]])                   # exponential distribution of alpha
grid_alpha=np.linspace(alpha_expo_range[0],alpha_expo_range[1],num=grid_num,endpoint=True) # split log(a) into grid_num points
alpha_list=10**grid_alpha

beta_expo_range=np.array([beta[0],beta[1]])
grid_beta=np.linspace(beta_expo_range[0],beta_expo_range[1],num=grid_num,endpoint=True) # split log(a) into grid_num points
beta_list=10**grid_beta

result=[[[] for i in range(grid_num)] for i in range(grid_num)]
all_result=[[[[np.nan] for i in range(grid_num)] for i in range(grid_num)] for dup in range(num_dup)]
np.shape(all_result[0][0][0])

for alpha_cluster in range(0,grid_num):                               # how many figures or cell divisions
	for beta_cluster in range(0,grid_num):                           # how many figures or cell divisions

		max_dup=[[] for i in range(num_dup)]
		for dup in range(0,num_dup):                               # how many figures or cell divisions
			max_dup[dup]=np.loadtxt('/Users/gao/Desktop/dds/SimulationCode/06_analytic/00_dynamic_bc/data_n5_m1000_sam100/%s_%s_%s_%s.txt'%(dup,n,alpha_cluster,beta_cluster))
			all_result[dup][alpha_cluster][beta_cluster]=max_dup[dup]

		"find the maximum over duplicates"
		grate_list=[]
		for dup in range(0,num_dup):
			grate_list.append(max_dup[dup][n,0])

		"save maximum"
		result[alpha_cluster][beta_cluster]=max_dup[np.argmax(grate_list)]

"""classify the duplicates: store the 10 duplicates' data"""
'''categaries: ND -1, RD=0; igd=1; isd=2, id=3  g->gg=1 no differentiation;'''
all_result_classify=[[[0 for i in range(grid_num)] for i in range(grid_num)]for dup in range(num_dup)]   # set all as DD
all_result_grate=[[[0 for i in range(grid_num)] for i in range(grid_num)]  for dup in range(num_dup)]

for dup in range(0,num_dup):                               # how many figures or cell divisions
	for alpha_cluster in range(0,grid_num):                               # how many figures or cell divisions
		for beta_cluster in range(0,grid_num):                           # how many figures or cell divisions

			item =all_result[dup][alpha_cluster][beta_cluster]
			all_result_grate[dup][beta_cluster][alpha_cluster]=item[n][0]

			item1=item[:n,:]
			for i in range(1,n+1):
#			for i in range(1,2):
				# ND
				if all(x==1 for x in item1[:,5]):                          # check whether last column f_g->gg are 1
					all_result_classify[dup][beta_cluster][alpha_cluster]=-1  
				# ISD	
				elif all(x==1 for x in item1[-i:,0]) and item1[-1:,5]!=1: 
					all_result_classify[dup][beta_cluster][alpha_cluster]=2 
				# IGD
				elif all(x==1 for x in item1[-i:,5]) and item1[-1:,0]!=1: 
					all_result_classify[dup][beta_cluster][alpha_cluster]=1 	
				# ID
				elif all(x==1 for x in item1[-i:,0])and all(x==1 for x in item1[-i:,5]) : 
					all_result_classify[dup][beta_cluster][alpha_cluster]=3 

all_classify_narry=np.array([[np.array(i) for i in item] for item in all_result_classify])
all_grate_narry=np.array([[np.array(i) for i in item] for item in all_result_grate])

"percentage of 10"
class_percent_nd=[[[] for i in range(grid_num)] for i in range(grid_num)]
class_percent_rd=[[[] for i in range(grid_num)] for i in range(grid_num)]
class_percent_isd=[[[] for i in range(grid_num)] for i in range(grid_num)]
class_percent_igd=[[[] for i in range(grid_num)] for i in range(grid_num)]
class_percent_id=[[[] for i in range(grid_num)] for i in range(grid_num)]

for alpha_cluster in range(0,grid_num):                               # how many figures or cell divisions
	for beta_cluster in range(0,grid_num):                           # how many figures or cell divisions
		dup_num=np.count_nonzero(~np.isnan(all_classify_narry[:,alpha_cluster,beta_cluster]))
		class_percent_nd[alpha_cluster][beta_cluster]=np.count_nonzero(all_classify_narry[:,alpha_cluster,beta_cluster]==-1)/dup_num #ND fraction
		class_percent_rd[alpha_cluster][beta_cluster]=np.count_nonzero(all_classify_narry[:,alpha_cluster,beta_cluster]==0)/dup_num #dd
		class_percent_isd[alpha_cluster][beta_cluster]=np.count_nonzero(all_classify_narry[:,alpha_cluster,beta_cluster]==2)/dup_num #isd
		class_percent_igd[alpha_cluster][beta_cluster]=np.count_nonzero(all_classify_narry[:,alpha_cluster,beta_cluster]==1)/dup_num #isd
		class_percent_id[alpha_cluster][beta_cluster]=np.count_nonzero(all_classify_narry[:,alpha_cluster,beta_cluster]==3)/dup_num #isd

class_percent_nd_narry=np.array([np.array(i) for i in class_percent_nd])
class_percent_rd_narry=np.array([np.array(i) for i in class_percent_rd])
class_percent_isd_narry=np.array([np.array(i) for i in class_percent_isd])
class_percent_igd_narry=np.array([np.array(i) for i in class_percent_igd])
class_percent_id_narry=np.array([np.array(i) for i in class_percent_id])

with open('dynamic_nd.txt', 'wb') as f:
	np.savetxt(f, class_percent_nd_narry, fmt='%1.8f')                     # demical number
with open('dynamic_rd.txt', 'wb') as f:
	np.savetxt(f, class_percent_rd_narry, fmt='%1.8f')                     # demical number
with open('dynamic_isd.txt', 'wb') as f:
	np.savetxt(f, class_percent_isd_narry, fmt='%1.8f')                     # demical number
with open('dynamic_igd.txt', 'wb') as f:
	np.savetxt(f, class_percent_igd_narry, fmt='%1.8f')                     # demical number
with open('dynamic_igsd.txt', 'wb') as f:
	np.savetxt(f, class_percent_id_narry, fmt='%1.8f')                     # demical number

