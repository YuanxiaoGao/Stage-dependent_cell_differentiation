#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 13:29:38 2022

@author: gao
"""

'''
f: fraction of germ cells
p1: p_{g->s}
p2: p_{s->g}

Parameters and Variables:
	Para: parameters for the cell dibision rate --- to be determined
		np.array([[benefit parameters],[cost parameters]]) (np array: float elements)
		np.array([[b,alpha,...],[c,beta,...]])
	Comp: [# of Soma cells, # of Germ cells] (np array: int elements)
	T: time to division (float) 	
	DP: six developement probabilities for cells
		[S->SS, S->GS, S->GG, G->SS, G->GS, G->GG] (np array: float elements)
	n: cell division times
	dp0: initial dtrategy
	ddp0: initial dynamic strategy
    SAM: sampling size of neighbourhood dynamic strategies that started with dp0	
	M: initial samplying size of dp0
	RD,IGD,ISD,IGSD : stra_cluster 0,1,2,3
'''
"""!!!!!!!here , we let  the delta be any value, that is the D1 AND D2 don't have ant constraint"""

import sys
import numpy as np
#------------------------------------------------------------------------------------------------------------
''' inital conditions;'''
       #-------- sys.argv
num_duplicate=int(sys.argv[1])                      # number of duplicates
n_cluster=int(sys.argv[2])                          # division times
alpha_cluster=int(sys.argv[3])                      # [0.01,100]
stra_cluster=int(sys.argv[4])                       # [0.01,100] 

       #-------- sys to local
b_points=[0.1,0.1,60,60]
c_points=[0.1,60,0.1,60]

n=n_cluster
alpha=b_points[alpha_cluster]    
beta=c_points[alpha_cluster]

		# ------set para	
"For calculating convenience, alpha is actually b and beta is c in the following setting. alpha=beta=1."
para=np.array([[alpha,1],[beta,1]])                # para=n	np.array([[b,alpha,...],[c,beta,...]])

		# ------set simulation times and dps
portion=[810,90,90,10]                              #total 1000	
M=portion[stra_cluster]*n**2		                  # the number of global initial dp   1000*n**2
SAM=np.int(n**2)                                   #sampling size of neighbourhood dynamic strategies	

#===========NUll part================growth rate formula============================
"""return array of switching probabilities p1 and p2; return the fraction of g and s after each i;
[    [np.array([p1,p2]),...],[[np.array([fg,fs]),....]]      ]"""
def Fration(n,dp):                              # dp is a sp series with number n
	p12=[np.array([0,0])]                       #p_g->s & p_s->g     total n+1 items
	fgs=[np.array([1,0])]                       #total n+1 items
	for i in range(n):
		p1=dp[i][1][0]+0.5*dp[i][1][1]
		p2=dp[i][0][2]+0.5*dp[i][0][1]
		p12.append(np.array([p1,p2]))
	for i in range(n+1):
		matrix=np.array([[1-p12[i][0],p12[i][1]],[p12[i][0],1-p12[i][1]]])
		fgs.append(matrix.dot(fgs[i]))
	del fgs[0]
	
	return [p12,fgs]                     

def Growth_num(n,dp,para):
	p12,fgs=Fration(n,dp)
	t=0
	for i in range(1,n+1):
		cost=1+para[1][0]*(fgs[i-1][0]*p12[i][0]+para[1][1]*fgs[i-1][1]*p12[i][1])
		bene=1+para[0][0]*(fgs[i-1][1])**para[0][1]
		t_i=cost/bene
		t=t+t_i
	growth_rate=np.log(2**n*fgs[n][0])/t
	return growth_rate

'''Input a initial dp for the first cell division; Return a dynmaic strategy (includes n dps). dp0->dp1->dp2->dp3->dp4..'''
'''if the delta is any value, then dynamic means choose from pool'''

def Dynamic_strategy(dp0,n):
	
	dpn=[dp0]	
	neighbour_dps=np.array(Dp_space())
	random_num=np.random.randint(n_dp,size=n-1)
	for i in range(n-1):				
		dpn.append(neighbour_dps[random_num[i]])
		
	return dpn[::-1]

#------------------------------------------------------------------------------------------------------------
'''Return SAM dynamic strategies including n dps started from dp0 (for the first cell division). 13**(n=5)=371293'''

def Dynamic_stra_set(dp0,SAM,n):
	dy_sta_set=[]
		
	while len(dy_sta_set)<SAM:
		one_sta=Dynamic_strategy(dp0,n)		
		test=np.array(one_sta)
		# for IGD exclude the ND
		if stra_cluster>0 and all(item==1 for item in test[:,1,2]): 
			pass
		else:
			one_sta=np.array(one_sta)	
			dy_sta_set.append(one_sta)
			
	return 	dy_sta_set                          # shape: sam * n * 2 *3

#------------------------------------------------------------------------------------------------------------
'''Find the optimal dynamic strategy among a set of dynamic strategies.'''

def Dp_lambda(dp0,para,SAM,n):
	
	growth_table=[]
	dynamic_strategy=Dynamic_stra_set(dp0,SAM,n)            # A randomly chose SAM dynmaic strategies for the following calculation.

	for i in range(len(dynamic_strategy)):                  # dp shape: n * 2 *3
		dp=dynamic_strategy[i]
		lambda0=Growth_num(n,dp,para)
		growth_table.append(np.array([i,lambda0]))         # insert lambda at the end
	
	array_data=np.array(growth_table)
	growth_max=np.nanmax(array_data[:,1])			   	 # maximum lambda
	index_growth_max=np.where(array_data==growth_max)      # index for maximum
	dp_max=dynamic_strategy[index_growth_max[0][0]]        #  dp for maximum  :np.ndarray
	
	return [dp_max,growth_max]

#-------------------global optimal-----------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------
'''All dps' DP=np.array([[s_ss,s_gs,s_gg];[g_ss,g_gs,g_gg]]); totally 11*12/2=66*66=4356'''
# we choose the random initial optimal strategy here

def Dp_space():
	size=np.int(11)          	         		         # how many discrete dp = [size*(size+1)/2]^2 
	q=np.linspace(0.0,1.0,size)                       # take sequence values between 0-1
	for i in range(size):
		q[i]=round(q[i],2)
		
	initial_para=[]                                  
	for i in range(size):
		for j in range(size):
			for k in range(size):
				for l in range(0,size):  
					if round(1-q[i]-q[j],2)>=0 and round(1-q[k]-q[l],2)>=0:
						dpi=np.round([[q[i],1-q[i]-q[j],q[j]],[q[k],1-q[k]-q[l],q[l]]],decimals=1)						
						initial_para.append(np.array(dpi))
	return initial_para                              # return a list

n_dp=4356	

dp_pool=[[] for i in range(4)]			#RD=4225; IGD=65; ISD=65; IGSD=1
for stra in Dp_space():
	if stra[0][0]==stra[1][2]==1:
		dp_pool[3].append(stra)       # IGSD
	elif stra[0][0]==1 and stra[1][2]<1:
		dp_pool[2].append(stra)
	elif stra[0][0]<1 and stra[1][2]==1:
		dp_pool[1].append(stra)
	else:
		dp_pool[0].append(stra)

#------------------------------------------------------------------------------------------------------------
'''Choose M initial DP=np.array([[s_ss,s_gs,s_gg];[g_ss,g_gs,g_gg]]) from the totally 
11*12/2=66*66=4356 dps; Return an array which is the optial initial dynamic strategy 
and its growth rate np.array([ ddp_array, grate ])'''

def Golbal_initial_ds(M,para,n,SAM):
	poo_size=len(dp_pool[stra_cluster])
	initial_num=min(poo_size,M)
	
	M_index=np.random.randint(poo_size,size=initial_num)               # choose M intial numbers
		
	initial_para=dp_pool[stra_cluster]	
	initial_dps=[]
	for i in M_index:
		initial_dps.append(initial_para[i])               # M initial dp0 (2*3)
				
	neighbour_optimal=[]
	for dp0 in initial_dps:
		optimal_dp_grate=Dp_lambda(dp0,para,SAM,n)   # each dp0 find the best one in its neighbour
		neighbour_optimal.append(optimal_dp_grate)		
	neighbour_optimal=np.array(neighbour_optimal)
	max_grate=np.nanmax(neighbour_optimal[:,1])
	op_index=np.where(neighbour_optimal[:,1]==max_grate)
		
	# return one	
	return neighbour_optimal[op_index[0][0]]            # np.ndarray    np.aray([ddp0,grate])
	
##------------------------------------------------------------------------------------------------------------
'''output result'''

def Output_global_ddp(M,para,SAM,n):
	data_ddp_grate=Golbal_initial_ds(M,para,n,SAM) 
	num_global_ddp=len(data_ddp_grate)
	result=[]
	if num_global_ddp==2:
		for i in data_ddp_grate[0]:
			result.append(i.flatten())	
		result.append(np.array([data_ddp_grate[-1],0,0,0,0,0]))	 # add grate at the last row
		result=np.array(result)
	else:
			
		for i in range(num_global_ddp-1):
			ddp_g=data_ddp_grate[i]
			for j in ddp_g:
				result.append(j.flatten())
		result.append(np.array([data_ddp_grate[-1],0,0,0,0,0]))	 # add grate at the last row
		result=np.array(result)
	return np.array(result)                                          # output np.array([ddp,grate]) or np.array([ddp,ddp, grate]) ddp are flat
	
result=Output_global_ddp(M,para,SAM,n)

with open('data/%s_%s_%s_%s.txt'%(num_duplicate,n_cluster,alpha_cluster,stra_cluster), 'wb') as f:
	np.savetxt(f, result, fmt='%1.8f')                     # demical number
	
