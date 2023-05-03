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
dp_all=Dp_space()
n=10
#------------------------------------------------------------------------------------------------------------
"get the random five strategies"
n=10
run_num=20
nd0=np.array([[0. , 1. , 0. ],[0. , 0, 1]])	
nd=np.array([nd0 for i in range(n)])

def stra():	
	# random dp0
	dp_list=[]
	dp_list.append(nd)

	num_rd=np.random.randint(0,n_dp)
	dp0=dp_all[num_rd]
	while dp0[0][0]==1 or dp0[1][2]>=1:                     #rd
		num_rd1=np.random.randint(0,n_dp)
		dp0=dp_all[num_rd1]
	dp_list.append(np.array([dp0 for i in range(n)]))
	
	num_rd1=np.random.randint(0,n_dp)
	dp1=dp_all[num_rd1]
	while dp1[0][0]!=1 or dp1[1][2]==1:                     #isd
		num_rd1=np.random.randint(0,n_dp)
		dp1=dp_all[num_rd1]		
	dp_list.append(np.array([dp1 for i in range(n)]))

	return dp_list

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
fg_data=[[] for stra in range(3)]
fg_mean=[[] for stra in range(3)]
fg_std=[[] for stra in range(3)]

r_data=[[] for stra in range(3)]
r_mean=[[] for stra in range(3)]
r_std=[[] for stra in range(3)]

grate_data=[[] for stra in range(3)]
grate_mean=[[] for stra in range(3)]
grate_std=[[] for stra in range(3)]

for run in range(run_num):
	dps=stra()
	grate_data[0].append(np.log(2))
	for i in range(3):
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
for i in range(3):
	fg_mean[i].append(np.mean(fg_data[i],axis=0))
	fg_std[i].append(np.std(fg_data[i],axis=0))

r_data=[np.array(i) for i in r_data]	
for i in range(3):
	r_mean[i].append(np.nanmean(r_data[i],axis=0))
	r_std[i].append(np.nanstd(r_data[i],axis=0))


grate_data=[np.array(i) for i in grate_data]	
grate_data=np.array(grate_data)
for i in range(3):
	grate_mean[i].append(np.nanmean(grate_data[i],axis=0))
	grate_std[i].append(np.nanstd(grate_data[i],axis=0))

#-------------DRAW FIGURE A B C D----------------------------------------------------------------
"figures"
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter

color_list=['#238b45','#2171b5','#8B4513','#984ea3','#ff7f00','#cb181d']		
	
fig, (ax1,ax2,ax3) = plt.subplots(1, 3, gridspec_kw={'width_ratios': [2,2,1]},figsize=(7.9, 4.5))
plt.rcParams["figure.autolayout"] = True
plt.rcParams['axes.linewidth'] =1

from matplotlib.patches import Patch
# ARTIFICIAL legend
# add a big axes, hide frame
fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
label_list=[r'$ND$',r'static $RD$',r'static $ID$ ($ISD$)']
legend_elements =[]
for i in range(3):
	legend_elements.append(Patch(facecolor=color_list[i], edgecolor='w',label=label_list[i]))
plt.legend(handles=legend_elements,ncol=3, loc='upper center', bbox_to_anchor=(0.5, 1.12),fontsize=11,frameon=False)

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
ax3.set_xlabel('Static strategy',fontsize=14, labelpad=14)
ax3.set_xticks([] )

ax1.set_ylabel(r'Frequency of germ-like cells, $f_{g}^{(i)}$',fontsize=14, labelpad=0)
ax2.set_ylabel(r'Cell division rate, $r^{(i)}$',fontsize=14, labelpad=0)
ax3.set_ylabel(r'Population growth rate, $\lambda^{S}$',fontsize=14, labelpad=0)

alp_small=0.6
'fig1'
x = np.linspace(0,n,n+1)
# lines
for i in range(3):
	for run in range(run_num):
		ax1.plot(x,fg_data[i][run] , linestyle = '-',linewidth=0.3,color=color_list[i],alpha=alp_small)	
		ax1.scatter(x,fg_data[i][run] ,s=1,marker="o",color=color_list[i],alpha=alp_small)
# shades
for i in range(3):
	ax1.scatter(x,fg_mean[i][0] ,s=50,marker="o",color=color_list[i],alpha=1)
	ax1.plot(x,fg_mean[i][0] , linestyle = '-',linewidth=3,color=color_list[i],alpha=1)		
	ax1.fill_between(x,fg_mean[i][0]-fg_std[i][0],fg_mean[i][0]+fg_std[i][0],
				 label=label_list[i],color=color_list[i],alpha=0.1)
#------x and y labels and limits--------
ax1.set_xlim(-0.1,10.5)
ax1.set_ylim(-0.1,1.05)

'fig2'
# lines
for i in range(3):
	for run in range(run_num):
		ax2.plot(x,r_data[i][run] , linestyle = '-',linewidth=0.3,color=color_list[i],alpha=alp_small)	
		ax2.scatter(x,r_data[i][run] ,s=1,marker="o",color=color_list[i],alpha=alp_small)
# shades
for i in range(3):
	ax2.scatter(x,r_mean[i][0] ,s=50,marker="o",color=color_list[i],alpha=1)
	ax2.plot(x,r_mean[i][0] , linestyle = '-',linewidth=3,color=color_list[i],alpha=1)	
	ax2.fill_between(x,r_mean[i][0]-r_std[i][0],r_mean[i][0]+r_std[i][0],
				 label=label_list[i],color=color_list[i],alpha=0.1)

'fig3'
x1 = np.linspace(0,3,3)
grate_mean[0].append(np.log(2))
grate_std[0].append(0)
# dots
for i in range(3):
	for run in range(run_num):
		ax3.scatter(x1[i],grate_data[i][run] ,s=2,marker="o",label=label_list[i],color=color_list[i],alpha=1)
# shades
for i in range(3):
	ax3.errorbar(x1[i],grate_mean[i][0],grate_std[i][0],linestyle='-',linewidth=1,marker="h",
			 markersize=12,label=label_list[i],
		   color=color_list[i],alpha=0.8 )
ax3.set_xlim(-1,4)
ax3.set_ylim(-3,0.99)

plt.show()
fig.savefig('./freuency_rate_grate_random_static-01.pdf', bbox_inches = 'tight')   # save figures



