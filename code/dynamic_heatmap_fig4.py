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
'''
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import random

#------------------------------------------------------------------------------------------------------------
''' inital conditions;'''
       #-------- sys.argv
num_duplicate=int(sys.argv[1])                      # number of duplicates
n_cluster=int(sys.argv[2])                          # division times
alpha_cluster=int(sys.argv[3])                      # [0.01,100]
beta_cluster=int(sys.argv[4])                       # [0.01,100]

#-------- sys to local	   
n=n_cluster
grid_num=41	   	      
alpha_expo_range=np.array([-1,3])                   # exponential distribution of alpha
grid_alpha=np.linspace(alpha_expo_range[0],alpha_expo_range[1],num=grid_num,endpoint=True) # split log(a) into grid_num points
alpha_list=10**grid_alpha
alpha=alpha_list[alpha_cluster]    

beta_expo_range=np.array([-1,3])                      
grid_beta=np.linspace(beta_expo_range[0],beta_expo_range[1],num=grid_num,endpoint=True) # split log(a) into grid_num points
beta_list=10**grid_beta
beta=beta_list[beta_cluster]

# ------set para	
"For calculating convenience, alpha is actually b and beta is c in the following setting. alpha=beta=1."
para=np.array([[alpha,1],[beta,1]])                       # para=n	np.array([[b,alpha,...],[c,beta,...]])

		# ------set simulation times and dps
M=1000                                               # the number of global initial dp
SAM=np.int(100)                                     #sampling size of neighbourhood dynamic strategies	

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

#===========First part================Neighbour optimal============================
'''for a randomly chose intial dp=np.array([[s_ss,s_gs,s_gg];[g_ss,g_gs,g_gg]]), find all its neighbours including itself
Return a list.'''

def Neighbour_dps(dp):
	delta_dp=np.array([[1,-1,0],[1,0,-1],[-1,1,0],[0,1,-1],[-1,0,1],[0,-1,1]])
	dp_set=[dp]
			
	for i in range(2):                          # for +1 increase switch, 0 no change, -1 decrease switch 
		for j in range(6):                     # for +1 increase switch, 0 no change, -1 decrease switch 			
			if i==0:				
				new0=dp[i]+0.1*delta_dp[j]
				new1=dp[1]
				new=np.array([new0,new1])
			else:
				new1=dp[i]+0.1*delta_dp[j]
				new0=dp[0]
				new=np.array([new0,new1])
				
			if np.all(new>=0):
				dp_set.append(new)
	return dp_set

#------------------------------------------------------------------------------------------------------------
'''Input a initial dp for the first cell division; Return a dynmaic strategy (includes n dps). dp0->dp1->dp2->dp3->dp4..'''

def Dynamic_strategy(dp0,n):
	
	dpn=[dp0]
	initial_dp=dp0
	
	for i in range(n-1):		
		neighbour_dps=np.array(Neighbour_dps(initial_dp))
		random_num=random.randint(0,len(neighbour_dps)-1)
		random_neighbour=neighbour_dps[random_num]
		dpn.append(random_neighbour)
		initial_dp=random_neighbour
		
	return dpn
	
#------------------------------------------------------------------------------------------------------------
'''Return SAM dynamic strategies including n dps started from dp0 (for the first cell division). 13**(n=5)=371293'''

def Dynamic_stra_set(dp0,SAM,n):
	dy_sta_set=[]
	
	for i in range(SAM):
		one_sta=Dynamic_strategy(dp0,n)
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

#===========Senond part================Local dynamic optimal============================

'''Find the optimal dynamic strategy among the neighbourhood dps.
dp_squeue: is a list of dp-squeues. [np.array([dp1,dp2,....],np.array([dp1,dp2,....],... )]'''
											  
def Dp_lambda_set(neighbour_dp_squeue,grate,para,n):
		
	growth_table=[np.array([0,grate])]
	
	for i in range(1,len(neighbour_dp_squeue)):              # calculate growth rate of its neighbours dp shape: n * 2 *3
		dp=neighbour_dp_squeue[i]
		lambda0=Growth_num(n,dp,para)
		growth_table.append(np.array([i,lambda0]))           # insert lambda at the end
	
	growth_table=np.array(growth_table)
	# find the optimal dp_squeue	
	growth_max=np.nanmax(growth_table[:,1])			   	   # maximum lambda
	index_growth_max=np.where(growth_table==growth_max)       # index for maximum
	dp_max=neighbour_dp_squeue[index_growth_max[0][0]]        #  dp for maximum  :np.ndarray
		
	return [dp_max,growth_max]

#------------------------------------------------------------------------------------------------------------
'''Find the optimal one on one left or right, i.e. [dp_new, dp2,dp3,dp4,dp5,dp_new], where the dp_new is one; '''

# right neighbour 
def Optimal_one_step(ddp0,grate,para,n):            # the initial dynamic dp
	neighbour_dps=np.array(Neighbour_dps(ddp0[-1]))             # all neighbours of the rightmost dp
	
	neighbour_dy_set=[ddp0]                                     # all dynamic strategies include one more step to the right
	for i in range(len(neighbour_dps)):
		new=neighbour_dps[i]
		new_dyna=	np.append(ddp0[1:],[new],axis=0)
		neighbour_dy_set.append(new_dyna)
	
	optimal=Dp_lambda_set(neighbour_dy_set,grate,para,n)    # get the optimal dp_queue among neighbour
	return optimal                                             # [dp_max, growth_rate] 

# left neighbour 
def Optimal_one_step_l(ddp0,grate,para,n):        # the initial dynamic ddp0
			
	neighbour_dps=np.array(Neighbour_dps(ddp0[0]))             # the neighbour of the leftmost dp
	
	neighbour_dy_set=[ddp0]                                   # all dynamic strategies
	for i in range(len(neighbour_dps)):
		new=neighbour_dps[i]
		new_dyna=	np.append([new],ddp0[:-1],axis=0)
		neighbour_dy_set.append(new_dyna)
	
	optimal=Dp_lambda_set(neighbour_dy_set,grate,para,n)	                      # get the optimal dp_queue among neighbour
	return optimal                                                 # [dp_max, growth_rate] 

#------------------------------------------------------------------------------------------------------------
'''The optimal dynamic strategy at equilibrium; Return: [  [ddps] ,  [grate_list]  ]'''

def Local_Optimal(ddp0,grate,para,n):
	
	# ddy and grate
	dyna_stra_list=[ddp0]
	dyna_rate_list=[grate]
	
	# add into list for optimal ddp and growth rate after one step
	one_step=Optimal_one_step(ddp0,grate,para,n)                      # one step for initial dyna_dp in right direction
	one_step_l=Optimal_one_step_l(ddp0,grate,para,n)                  # one step for initial dyna_dp in left direction
	
	if one_step[1]>one_step_l[1]:                                     # go to the right direction
		dyna_stra_list.append(one_step[0])
		dyna_rate_list.append(one_step[1])
		
		while abs(dyna_rate_list[-1]- dyna_rate_list[-2])>10**(-8):
			new_optimal_dp=dyna_stra_list[-1]
			grate=dyna_rate_list[-1]
			new_optimal=Optimal_one_step(new_optimal_dp,grate,para,n)
			
			dyna_stra_list.append(new_optimal[0])
			dyna_rate_list.append(new_optimal[1])
	else:		                                                    # go to the left direction
		dyna_stra_list.append(one_step_l[0])
		dyna_rate_list.append(one_step_l[1])
	
		while abs(dyna_rate_list[-1]- dyna_rate_list[-2])>10**(-8):    # dyna_rate_list[-1]!=dyna_rate_list[-2]:
			new_optimal_dp=dyna_stra_list[-1]
			grate=dyna_rate_list[-1]
			new_optimal=Optimal_one_step_l(new_optimal_dp,grate,para,n)
			
			dyna_stra_list.append(new_optimal[0])
			dyna_rate_list.append(new_optimal[1])
	if len(dyna_rate_list)==2:
		return [[ddp0],[grate]]	
	else:		
		return [dyna_stra_list,dyna_rate_list]	

#===========Third part================Global optimal ddp============================
'''All dps' DP=np.array([[s_ss,s_gs,s_gg];[g_ss,g_gs,g_gg]]); totally 11*12/2=66*66=4356'''

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
						initial_para.append(np.array([[q[i],1-q[i]-q[j],q[j]],[q[k],1-q[k]-q[l],q[l]]]))
	return initial_para                              # return a list

n_dp=4356									
	
#------------------------------------------------------------------------------------------------------------
'''Choose M initial DP=np.array([[s_ss,s_gs,s_gg];[g_ss,g_gs,g_gg]]) from the totally 
11*12/2=66*66=4356 dps; Return an array which is the optial initial dynamic strategy 
and its growth rate np.array([ ddp_array, grate ])'''

def Golbal_initial_ds(M,para,n,SAM):
	M_index=np.random.randint(n_dp-1,size=M)               # choose M intial numbers
	
	initial_para=Dp_space()	
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
	
	if 	np.shape(op_index[0])[0]==1 or n==1:		
		return neighbour_optimal[op_index[0][0]]            # np.ndarray    np.aray([ddp0,grate])
	
	else:
		multiple_ddy0=[]
		for index in op_index[0]:
			ddp_flatten=np.array([neighbour_optimal[index][0][i].flatten()	 for i in range(n)])
			multiple_ddy0.append(ddp_flatten)
		multiple_ddy0=np.array(multiple_ddy0)
		
		if np.all(multiple_ddy0[:,:,-1].flatten()==1):
			return neighbour_optimal[op_index[0][0]]	      # np.ndarray    np.aray([ddp0,grate])
		else:
			result_list=[]
			for index in op_index[0]:
				result_list.append(neighbour_optimal[index][0])
				
			result_list.append(max_grate)
			return np.array(result_list)       # two or more ddp0 and with grate  np.aray([ddp0,ddp0,grate])

#------------------------------------------------------------------------------------------------------------
'''Find the local optimal for the global intial optimal dynamic strategy'''

def Golbal_ddp(M,para,SAM,n):
	
	local_ddp_grate=Golbal_initial_ds(M,para,n,SAM)
	grate=local_ddp_grate[-1]
	num=np.shape(local_ddp_grate)[0]
	if num==2:
		ddp0=local_ddp_grate[0]
		history_dynamic=Local_Optimal(ddp0,grate,para,n)
	
		return [history_dynamic[0][-1],history_dynamic[1][-1]]   # return a list [ddp,grate]
	else:
		ddp_list=[]
		grate_list=[]
		for i in range(num-1):
			ddp0=local_ddp_grate[i]
			history_dynamic=Local_Optimal(ddp0,grate,para,n)
			ddp_list.append(history_dynamic[0][-1])
			grate_list.append(history_dynamic[1][-1])
		grate_list=np.array(grate_list)		
		max_grate=np.nanmax(grate_list)
		op_index=np.where(max_grate==max_grate)
		
		if np.shape(op_index[0])[0]==1:                  # if there is only one global ddp		
			return [ddp_list[op_index],max_grate]        # return a list [ddp,grate]
		else:
			ddp_global_list=[]
			for item in op_index[0]:
				ddp_global_list.append(ddp_list[item])
			ddp_global_list.append(max_grate)	
			return ddp_global_list                      # return a list [ddp,ddp,ddp, grate]    

##------------------------------------------------------------------------------------------------------------
'''output result'''

def Output_global_ddp(M,para,SAM,n):
	data_ddp_grate=Golbal_ddp(M,para,SAM,n) 
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

with open('data/%s_%s_%s_%s.txt'%(num_duplicate,n_cluster,alpha_cluster,beta_cluster), 'wb') as f:
	np.savetxt(f, result, fmt='%1.8f')                     # demical number
	
