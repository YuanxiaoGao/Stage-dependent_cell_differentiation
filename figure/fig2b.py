#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 13:29:38 2022

@author: gao
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import random

np.random.seed(6)
#-------------------------------------------------------------------------------
# set display width
np.get_printoptions()['linewidth']
np.set_printoptions(linewidth=2000)
para=np.array([[1,1],[1,1]])                    # para=n	np.array([[b,alpha,...],[c,beta,...]])

#------------------------------------------------------------------------------------------------------------
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
						dpi=np.round([[q[i],1-q[i]-q[j],q[j]],[q[k],1-q[k]-q[l],q[l]]],decimals=1)						
						initial_para.append(np.array(dpi))
	return initial_para                              # return a list

n_dp=4356	

dp_pool=[[] for i in range(4)]			# RD; IGD; ISD; IGSD
for stra in Dp_space():
	if stra[0][0]==stra[1][2]==1:
		dp_pool[3].append(stra)       # IGSD
	elif stra[0][0]<1 and stra[1][2]==1:
		dp_pool[1].append(stra)        # IGD
	elif stra[0][0]==1 and stra[1][2]<1:
		dp_pool[2].append(stra)       # ISD
	else:
		dp_pool[0].append(stra)

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
			new2=np.round(new,decimals=1)
			
			if np.all(new2>=0):
				dp_set.append(new2)
				
	return dp_set

#------------------------------------------------------------------------------------------------------------
'''Input a initial dp for the first cell division; Return a dynmaic strategy (includes n dps). dp0->dp1->dp2->dp3->dp4..'''
def Dynamic_strategy(dp0,n,stra_cluster):
		
	dy_sta_set=[]
	while len(dy_sta_set)<1:
		
		dpn=[dp0]
		initial_dp=dp0
		for i in range(n-1):		
			neighbour_dps=np.array(Neighbour_dps(initial_dp))
			random_num=random.randint(0,len(neighbour_dps)-1)
			random_neighbour=neighbour_dps[random_num]
			dpn.append(random_neighbour)
			initial_dp=random_neighbour
		dpn=np.array(dpn)	
		
		if stra_cluster>0 and all(item==1 for item in dpn[:,1,2]): 
			pass
		else:
			dy_sta_set.append(dpn)
		
	return dpn[::-1]
"get the random five strategies"
n=10
run_num=20

num=[]
for i in range(4):
	num.append(len(dp_pool[i]))
	
nd=np.array([[[0,0,0],[0,0,1]],
			 [[0,0,0],[0,0,1]],
			 [[0,0,0],[0,0,1]],
			 [[0,0,0],[0,0,1]],
			 [[0,0,0],[0,0,1]]])

def stra():	
	# random dp0
	dp0_list=[]
	for i in range(4):
		dp0_list.append(random.randint(0,num[i]-1)	)	
	# dynamic dp	 
	dps=[nd]
	
	for i in range(4):
		test_dp0=dp_pool[i][dp0_list[i]]
		dps.append(np.array(Dynamic_strategy(test_dp0,n,i)))
	return dps

"get the frequency and time of the random five strategies"
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
	t_list=[0]
	t=0
	for i in range(1,n+1):
		cost=1+para[1][0]*(fgs[i-1][0]*p12[i][0]+para[1][1]*fgs[i-1][1]*p12[i][1])
		bene=1+para[0][0]*(fgs[i-1][1])**para[0][1]
		t_i=cost/bene
		t_list.append(t_i)
		t=t+t_i
	growth_rate=np.log(2**n*fgs[n][0])/t
	return [growth_rate,t_list]

"draw figures"
fg_data=[[] for stra in range(5)]
fg_mean=[[] for stra in range(5)]
fg_std=[[] for stra in range(5)]

r_data=[[] for stra in range(5)]
r_mean=[[] for stra in range(5)]
r_std=[[] for stra in range(5)]

grate_data=[[] for stra in range(5)]
grate_mean=[[] for stra in range(5)]
grate_std=[[] for stra in range(5)]

for run in range(run_num):
	dps=stra()
	grate_data[0].append(np.log(2))
	for i in range(1,5):
		dp=dps[i]
		fre=np.array(Fration(n,dp)[1])[:,0]
		time=np.array(Growth_num(n,dp,para)[1])
		grate=np.array(Growth_num(n,dp,para)[0])		
		fg_data[i].append(fre)
		r_data[i].append(1/time)
		r_data[i][run][0]=np.nan
		grate_data[i].append(grate)

"mean and std of frequency, glambda and division rate"		
fg_data=[np.array(i) for i in fg_data]	
for i in range(1,5):
	fg_mean[i].append(np.mean(fg_data[i],axis=0))
	fg_std[i].append(np.std(fg_data[i],axis=0))

r_data=[np.array(i) for i in r_data]	
for i in range(1,5):
	r_mean[i].append(np.nanmean(r_data[i],axis=0))
	r_std[i].append(np.nanstd(r_data[i],axis=0))


grate_data=[np.array(i) for i in grate_data]	
grate_data=np.array(grate_data)
grate_data[np.isinf(grate_data)]=np.nan
for i in range(1,5):
	grate_mean[i].append(np.nanmean(grate_data[i],axis=0))
	grate_std[i].append(np.nanstd(grate_data[i],axis=0))

#-------------DRAW FIGURE A B C D----------------------------------------------------------------
"figures"
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter

label_list=[r'$ND$',r'$RD$',r'$IGD$',r'$ISD$',r'$IGSD$']
color_list=['#238b45','#2171b5','#984ea3','#ff7f00','#cb181d']		
	
fig, (ax1,ax2,ax3) = plt.subplots(1, 3, gridspec_kw={'width_ratios': [2,2,1]},figsize=(7.9, 4.7))
plt.rcParams["figure.autolayout"] = True
plt.rcParams['axes.linewidth'] =1

from matplotlib.patches import Patch
# ARTIFICIAL legend
# add a big axes, hide frame
fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
legend_elements =[]

p_color_list=['w','#238b45','w','#2171b5','w','#984ea3','#8B4513','#ff7f00','w','#cb181d']	
p_label_list=['',r'$ND$', '', r'$RD$','',r'$IGD$', r'ID', r'$ISD$','',r'$IGSD$']

for i in range(10):
	legend_elements.append(Patch(facecolor=p_color_list[i], edgecolor='w',label=p_label_list[i]))
plt.legend(handles=legend_elements,ncol=5, loc='upper center', bbox_to_anchor=(0.5, 1.18),fontsize=11,frameon=False)

ax1.patch.set_edgecolor('r')  
ax1.patch.set_linewidth('0.01')  

# frame color
color_f="lightgrey"
ax1.spines['bottom'].set_color(color_f)
ax1.spines['top'].set_color(color_f) 
ax1.spines['right'].set_color(color_f)
ax1.spines['left'].set_color(color_f)

ax2.spines['bottom'].set_color(color_f)
ax2.spines['top'].set_color(color_f) 
ax2.spines['right'].set_color(color_f)
ax2.spines['left'].set_color(color_f)

ax3.spines['bottom'].set_color(color_f)
ax3.spines['top'].set_color(color_f) 
ax3.spines['right'].set_color(color_f)
ax3.spines['left'].set_color(color_f)

ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax3.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

ax1.tick_params(direction='in', length=2, width=1, colors='k',grid_color='r', grid_alpha=1)     # inout
ax2.tick_params(direction='in', length=2, width=1, colors='k',grid_color='r', grid_alpha=1)
ax3.tick_params(direction='in', length=2, width=1, colors='k',grid_color='r', grid_alpha=1)

ax1.locator_params(axis='y', nbins=5)
ax2.locator_params(axis='y', nbins=4)
#	ax1.tick_params(axis='y', which='major', pad=5)
#	ax2.tick_params(axis='y', which='major', pad=5)#	ax1.axes.get_xaxis().set_visible(False)
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
ax2.xaxis.set_major_locator(MaxNLocator(integer=True))

ax1.set_xlabel(r'Cell division times, $i$',fontsize=14, labelpad=0)
ax2.set_xlabel(r'Cell division times, $i$',fontsize=14, labelpad=0)
ax3.set_xlabel('Dynamic strategy',fontsize=14, labelpad=12)
ax3.set_xticks([] )

ax1.set_ylabel(r'Frequency of germ-like cells, $f_{g}^{(i)}$',fontsize=14, labelpad=0)
ax2.set_ylabel(r'Cell division rate, $r^{(i)}$',fontsize=14, labelpad=0)
ax3.set_ylabel(r'Population growth rate, $\lambda$',fontsize=14, labelpad=0)

alp_small=0.6
'fig1'
x = np.linspace(0,n,n+1)
# ND
nd_fra=np.array([1 for i in range(n+1)])
ax1.plot(x,nd_fra , linestyle = '-',linewidth=0.3,color=color_list[0],alpha=alp_small)	
ax1.scatter(x,nd_fra ,s=50,marker="o",label=label_list[0],color=color_list[0],alpha=1)
ax1.plot(x,nd_fra , linestyle = '-',linewidth=3,color=color_list[0],alpha=1)	
# lines
for i in range(1,5):
	for run in range(run_num):
		ax1.plot(x,fg_data[i][run] , linestyle = '-',linewidth=0.3,color=color_list[i],alpha=alp_small)	
		ax1.scatter(x,fg_data[i][run] ,s=1,marker="o",label=label_list[i],color=color_list[i],alpha=alp_small)
# shades
for i in range(1,5):
	ax1.scatter(x,fg_mean[i][0] ,s=50,marker="o",label=label_list[i],color=color_list[i],alpha=1)
	ax1.plot(x,fg_mean[i][0] , linestyle = '-',linewidth=3,color=color_list[i],alpha=1)		
	ax1.fill_between(x,fg_mean[i][0]-fg_std[i][0],fg_mean[i][0]+fg_std[i][0],
				 label=label_list[i],color=color_list[i],alpha=0.1)
#------x and y labels and limits--------
ax1.set_xlim(-0.1,10.5)
ax1.set_ylim(-0.1,1.05)

'fig2'
# ND
nd_fra=np.array([1 for i in range(n+1)])
ax2.plot(x[1:,],nd_fra[1:,] , linestyle = '-',linewidth=0.3,color=color_list[0],alpha=alp_small)	
ax2.scatter(x[1:,],nd_fra[1:,] ,s=1,marker="o",label=label_list[0],color=color_list[0],alpha=1)
ax2.plot(x[1:,],nd_fra[1:,] , linestyle = '-',linewidth=3,color=color_list[0],alpha=1)	
# lines
for i in range(1,5):
	for run in range(run_num):
		ax2.plot(x,r_data[i][run] , linestyle = '-',linewidth=0.3,color=color_list[i],alpha=alp_small)	
		ax2.scatter(x,r_data[i][run] ,s=1,marker="o",label=label_list[i],color=color_list[i],alpha=alp_small)
# shades
for i in range(1,5):
	ax2.scatter(x,r_mean[i][0] ,s=50,marker="o",label=label_list[i],color=color_list[i],alpha=1)
	ax2.plot(x,r_mean[i][0] , linestyle = '-',linewidth=3,color=color_list[i],alpha=1)	
	ax2.fill_between(x,r_mean[i][0]-r_std[i][0],r_mean[i][0]+r_std[i][0],
				 label=label_list[i],color=color_list[i],alpha=0.1)

'fig3'
x1 = np.linspace(0,5,5)
grate_mean[0].append(np.log(2))
grate_std[0].append(0)

# dots
for i in range(1,5):
	for run in range(run_num):
		ax3.scatter(x1[i],grate_data[i][run] ,s=2,marker="o",label=label_list[i],color=color_list[i],alpha=1)
# shades
for i in range(5):
	ax3.errorbar(x1[i],grate_mean[i][0],grate_std[i][0],linestyle='-',linewidth=1,marker="h",
			 markersize=12,label=label_list[i],
		   color=color_list[i],alpha=0.8 )
ax3.set_xlim(-0.5,5.5)
ax3.set_ylim(-3,0.99)

plt.show()

#fig.savefig('./fig2b.pdf', bbox_inches = 'tight')   # save figures

