#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 13:29:38 2022

@author: gao
"""

"""
f: fraction of germ cells
p1: p_{g->s}
p2: p_{s->g}
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
#-------------------------------------------------------------------------------
# set display width
np.get_printoptions()['linewidth']
np.set_printoptions(linewidth=2000)

#-------------------------------------------------------------------------------
# fraction of germ-role with respect to n
def fra(n,p1,p2):
	return ((1-p1-p2)**n*p1+p2)/(p1+p2)

# growth time 
def time(alpha,beta,b,c,p1,p2,n):	
	t=0
	for i in range(1,n+1):
		f=fra(i-1,p1,p2)
		t_i=(1+c*(f*p1+beta*(1-f)*p2))/(1+b*(1-f)**alpha)
		t=t+t_i
	return t

# growth rate
def grate(alpha,beta,b,c,p1,p2,n):
	f=fra(n,p1,p2)
	dno=np.log((2**n)*f)
	numo=time(alpha,beta,b,c,p1,p2,n)
	return (dno)/(numo)

'''Parameters'''
n=5
alpha=1
beta=1
	   	      
grid_num=41
bound_l=-1
bound_r=2		      
alpha_expo_range=np.array([bound_l,bound_r])                   # exponential distribution of alpha
grid_alpha=np.linspace(alpha_expo_range[0],alpha_expo_range[1],num=grid_num,endpoint=True) # split log(a) into grid_num points
b_list=10**grid_alpha

beta_expo_range=np.array([bound_l,bound_r])                      
grid_beta=np.linspace(beta_expo_range[0],beta_expo_range[1],num=grid_num,endpoint=True) # split log(a) into grid_num points
c_list=10**grid_beta

""" data for optimal strategies"""
sets=20
scale=1/sets

ab_matrix=[[0 for i in range(grid_num)] for j in range(grid_num)]    # 0 ND; 1 RD; 2 ISD
ab_growth=[[np.log(2) for i in range(grid_num)] for j in range(grid_num)]     # optimal growth rate 
ab_p_matrix=[[0 for i in range(grid_num)] for j in range(grid_num)]  # p1
ab_q_matrix=[[-1 for i in range(grid_num)] for j in range(grid_num)]  # p2

for i in range(grid_num):
	b=b_list[i]
	for j in range(grid_num):
		c=c_list[j]
		
		# maximal isd
		y_isd=[]
		for p in range(1,sets):                # don't let p1=1, because p1=1,p2=0 no offspring
			p1=scale*p 	
			y_isd.append(grate(alpha,beta,b,c,p1,0,n))
		y_isd=np.array(y_isd)	
		max_isd=np.nanmax(y_isd)
		index0=np.where(y_isd==max_isd)
		sta_pg_isd=[scale*(index0[0][0]+1)]     # p1 and p2
		
		# maximal rd
		y_rd=[[] for i in range(sets)]
		for p in range(1,sets+1):                    # 0<p_g->s=1 
			p1=scale*p 	
			for q in range(1,sets+1):               # p_s->g>0
				q1=scale*q 	
				
				y_rd[p-1].append(grate(alpha,beta,b,c,p1,q1,n))
		y_rd_array=np.array([np.array(item) for item in y_rd])
		max_rd=np.nanmax(y_rd)
		index1=np.where(y_rd==max_rd)
		sta_pq_rd=[scale*(index1[0][0]+1),scale*(index1[1][0]+1)]   # p1 and p2
		
		if np.max([np.log(2),max_isd,max_rd])==max_isd:
			ab_matrix[j][i]=2  
			ab_growth[j][i]=max_isd
			ab_p_matrix[j][i]=sta_pg_isd[0]
			ab_q_matrix[j][i]=np.nan         
		elif np.max([np.log(2),max_isd,max_rd])==max_rd:
			ab_matrix[j][i]=1   
			ab_growth[j][i]=max_rd
			ab_p_matrix[j][i]=sta_pq_rd[0]
			ab_q_matrix[j][i]=sta_pq_rd[1]        

static_matrix=np.array(ab_matrix)

with open('static_m.txt', 'wb') as f:
	np.savetxt(f, static_matrix, fmt='%1.0f')                     # demical number
