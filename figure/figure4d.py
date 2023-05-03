#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 16:41:30 2023

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

def data(x,y):
	for alpha_cluster in range(x,x+1):                               # how many figures or cell divisions
		for beta_cluster in range(y,y+1):                           # how many figures or cell divisions
	
			max_dup=[[] for i in range(num_dup)]
			for dup in range(0,num_dup):                               # how many figures or cell divisions
				max_dup[dup]=np.loadtxt('../data/data_fig4acd/%s_%s_%s_%s.txt'%(dup,n,alpha_cluster,beta_cluster))
				all_result[dup][alpha_cluster][beta_cluster]=max_dup[dup]
	
			"find the maximum over duplicates"
			grate_list=[]
			for dup in range(0,num_dup):
				grate_list.append(max_dup[dup][n,0])
	
			"save maximum"
			result[alpha_cluster][beta_cluster]=max_dup[np.argmax(grate_list)]
	return result[alpha_cluster][beta_cluster]

"choose the bc dots interested"
x=[3,20,38,30,38]
y=[1,1,1,28,28]

dp_list=[]
for i in range(len(x)):
	a=x[i]
	b=y[i]
	dp_list.append(data(a,b))

"get the transition probabilities of the random five strategies"
def Fration(n,dp):                              # dp is a sp series with number n
	p12=[np.array([0,0])]                       #p_g->s & p_s->g     total n+1 items
	fgs=[np.array([1,0])]                       #total n+1 items
	for i in range(n):
		p1=dp[i][3]+0.5*dp[i][4]
		p2=dp[i][2]+0.5*dp[i][1]
		p12.append(np.array([p1,p2]))
	for i in range(n+1):
		matrix=np.array([[1-p12[i][0],p12[i][1]],[p12[i][0],1-p12[i][1]]])
		fgs.append(matrix.dot(fgs[i]))
	del fgs[0]	
	return [p12,fgs]  

import matplotlib.pyplot as plt
import matplotlib.lines as mlines

label_list=['ND','RD','IGD','ISD','IGSD']
color_list=['#238b45','dimgrey','#2171b5',"k",'#8B4513'] ##984ea3','#ff7f00','#cb181d'
size=[i for i in range(0,6)]

fig, axs= plt.subplots(1,5, gridspec_kw={'width_ratios': [1,1,1,1,1]},figsize=(16, 3))
				
filled_o = mlines.Line2D([], [], color='#8B4513', marker='o', linestyle='-',
                          markersize=10, label=r'$g_{g \to s}$')
circle_o = mlines.Line2D([], [], color='#8B4513', marker='o', linestyle='--',markerfacecolor='w', 
								 markersize=10, label=r'$s_{s \to g}$')

"figs-lines"
frac_list=[]
fg_list=[]
i=0
for ax in axs.reshape(-1): 
	test1=dp_list[i]
	frac=np.array(Fration(n,test1)[0])
	fg=np.array(Fration(n,test1)[1])

	ax.plot(size,frac[:,0] , linestyle = '-',linewidth=1,color=color_list[i])	
	ax.scatter(size,frac[:,0],s=70,marker="o",color=color_list[i])
	ax.plot(size,frac[:,1] , linestyle = '--',linewidth=1,color=color_list[i])	
	ax.scatter(size,frac[:,1],s=70,marker="o",color=color_list[i],facecolors='none', edgecolors=color_list[i])

	ax.plot(size,fg[:,0] , linestyle = ':',linewidth=1,color=color_list[i])	
	ax.scatter(size,fg[:,0],s=60,marker="H",color=color_list[i])	

	ax.set_xlabel(r'Cell division times, $i$',fontsize=14, labelpad=0)	

	if i==0:
		ax.text(size[1]-0.05,fg[:,0][1]+0.05,r"$f_{g}$",fontsize=18,color=color_list[i])
	elif i==1:
		ax.text(size[4]-0.05,fg[:,0][4]+0.05,r"$f_{g}$",fontsize=18,color=color_list[i])
	elif i==2:
		ax.text(size[4]-0.05,fg[:,0][4]+0.08,r"$f_{g}$",fontsize=18,color=color_list[i])
	else:
		ax.text(size[2]-0.05,fg[:,0][2]+0.05,r"$f_{g}$",fontsize=18,color=color_list[i])

	if i==4:
		ax.legend(handles=[filled_o, circle_o],loc="center right",fontsize=14)

	i=i+1
plt.show()

#fig.savefig('./fig4d.pdf', bbox_inches = 'tight')   # save figures

