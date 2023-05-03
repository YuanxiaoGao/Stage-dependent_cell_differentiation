#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 14:40:40 2018

@author: gao
"""
import os
import sys
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
#-------------------------------------------------------------------------------
# set display width
np.get_printoptions()['linewidth']
np.set_printoptions(linewidth=2000)

#-------------------------------------------------------------------------------
'''Parameters'''
num_dup=20

b_points=[0.1,0.1,60,60]
c_points=[0.1,60,0.1,60]
#-------------------------------------------------------------------------------
'''read data '''
'''read data '''
grid_num=len(b_points)
result=[[[] for stra in range(4)] for b in range(len(b_points))]   # save the four categories

for b in range(0,grid_num):  
	for stra_cluster in range(4): 
		for dup in range(num_dup): 
			for n in range(5,6):  
		
				if os.path.exists('./data/%s_%s_%s_%s.txt'%(dup,n,b,stra_cluster)):
					data=np.loadtxt('./data/%s_%s_%s_%s.txt'%(dup,n,b,stra_cluster))						
					result[b][stra_cluster].append(data)
			
				else:
					print('%s_%s_%s_%s'%(dup,n,b,stra_cluster))


"get the frequency and time of the random five strategies"
  
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
#		print(t_i)
	growth_rate=np.log(2**n*fgs[n][0])/t
	return [growth_rate,t_list]	

#===================for extreme=======================
"""return array of switching probabilities p1 and p2; return the fraction of g and s after each i;
[    [np.array([p1,p2]),...],[[np.array([fg,fs]),....]]      ]"""
def Fration1(n,dp):                              # dp is a sp series with number n
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

#dp=[np.array([[1,0,0],[1,0,0]]),np.array([[1,0,0],[1,0,0]]),np.array([[1,0,0],[1,0,0]]),np.array([[1,0,0],[1,0,0]]),np.array([[0,0,1],[0,0,1]])]
def Growth_num1(n,dp,para):
	p12,fgs=Fration1(n,dp)
	t=0
	for i in range(1,n+1):
		cost=1+para[1][0]*(fgs[i-1][0]*p12[i][0]+para[1][1]*fgs[i-1][1]*p12[i][1])
		bene=1+para[0][0]*(fgs[i-1][1])**para[0][1]
		t_i=cost/bene
#		print(t_i)
		t=t+t_i
	growth_rate=np.log(2**n*fgs[n][0])/t
	return growth_rate


#-------------DRAW FIGURE A B C D----------------------------------------------------------------
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter

label_list=['ND','RD','IGD','ISD','IGSD','ED']
color_list=['#238b45','#2171b5','#984ea3','#ff7f00','#cb181d','k']

size=[i for i in range(0,6)]
		
for b_in in range(len(b_points)):  

	b=b_points[b_in]
	c=c_points[b_in]
	
	fig, (ax1,ax2,ax3) = plt.subplots(1, 3, gridspec_kw={'width_ratios': [2,2,1]},figsize=(7.5, 4))
	fig.suptitle(r"$0 \leq \delta \leq 1$, $b=%s$, $c=%s$"%(b,c),fontsize=14,y=1)

	plt.rcParams["figure.autolayout"] = True
	plt.rcParams['axes.linewidth'] =1
	
	ax1.patch.set_edgecolor('r')  
	ax1.patch.set_linewidth('0.01')  
	
	from matplotlib.patches import Patch
	# ARTIFICIAL legend
	# add a big axes, hide frame
	fig.add_subplot(111, frameon=False)
	# hide tick and tick label of the big axes
	plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
	plt.grid(False)
	legend_elements =[]
	for i in range(6):
		legend_elements.append(Patch(facecolor=color_list[i], edgecolor='w',label=label_list[i]))
	plt.legend(handles=legend_elements,ncol=6, loc='upper center', bbox_to_anchor=(0.5, 1.115),fontsize=10,frameon=False)
	
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
	ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
	ax3.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
	
	ax1.tick_params(direction='in', length=2, width=1, colors='k',
               grid_color='r', grid_alpha=1)     # inout
	ax2.tick_params(direction='in', length=2, width=1, colors='k',
               grid_color='r', grid_alpha=1)
	ax3.tick_params(direction='in', length=2, width=1, colors='k',
               grid_color='r', grid_alpha=1)
	
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
	marker_s=12
    # ND
	nd_fra=np.array([1 for i in range(6)])
	ax1.plot(size,nd_fra , linestyle = '-',linewidth=0.3,color=color_list[0],alpha=alp_small)	
	ax1.scatter(size,nd_fra ,s=50,marker="o",label=label_list[0],color=color_list[0],alpha=1)
	ax1.plot(size,nd_fra , linestyle = '-',linewidth=3,color=color_list[0],alpha=1)	
	
	nd_r=np.array([1,1,1,1,1])
	ax2.plot(size[1:],nd_r , linestyle = '-',linewidth=0.3,color=color_list[0],alpha=alp_small)	
	ax2.scatter(size[1:],nd_r ,s=50,marker="o",label=label_list[0],color=color_list[0],alpha=alp_small)
	ax2.plot(size[1:],nd_r , linestyle = '-',linewidth=3,color=color_list[0],alpha=1)	
	
	ax3.errorbar(0,np.log(2),0,linestyle='-',linewidth=1,marker="h",markersize=12,label=label_list[stra_cluster+1],
	   color=color_list[0],alpha=0.8 )
	
	# star for the extreme strategy
	ax1.plot(size,[1,0,0,0,0,1],linestyle = '-',linewidth=0.3,color="k")
	ax1.scatter(size,[1,0,0,0,0,1] ,s=3,marker="o",color="k")
	ax2.plot(size[1:],[1/(1+c),(1+b),(1+b),(1+b),(1+b)/(1+c)],linestyle = '-',linewidth=0.3,color="k")
	ax2.scatter(size[1:],[1/(1+c),(1+b),(1+b),(1+b),(1+b)/(1+c)],s=3,marker="o",color="k")
	dp_ex=[np.array([[1,0,0],[1,0,0]]),np.array([[1,0,0],[1,0,0]]),np.array([[1,0,0],[1,0,0]]),np.array([[1,0,0],[1,0,0]]),np.array([[0,0,1],[0,0,1]])]
	para=np.array([[b_points[b_in],1],[c_points[b_in],1]])                # para=n	np.array([[b,alpha,...],[c,beta,...]])
	extrem=Growth_num1(n,dp_ex,para)
	ax3.plot(2,extrem,marker="h",markersize=12,color="k",alpha=0.8 )

	# for others		
	for stra_cluster in range(4): 	
		all_stra=result[b_in][stra_cluster]

		fg_list=[]
		t_list=[]
		grate_list=[]
		for run in range(num_dup):
			one=all_stra[run]
			dp=one[0:5,:]
			grate=one[5][0]
			fg=np.array(Fration(n,dp)[1])[:,0]
			t=np.array(Growth_num(n,dp,np.array([[b,1],[c,1]]))[1])
			fg_list.append(fg)
			t_list.append(t)
			grate_list.append(grate)
#		"frquency of germ-like cells"	
		fg_np=np.array(fg_list)
		fg_mean=np.mean(fg_np,axis=0)
		fg_std=np.std(fg_np,axis=0)
#		"division rate r"		
		r_np=1/np.array(t_list)
		r_np[0]=np.nan
		r_mean=np.nanmean(r_np,axis=0)
		r_std=np.nanstd(r_np,axis=0)
#		"growthe rate lambda"
		grate_list=np.array(grate_list)
		grate_mean=np.nanmean(grate_list,axis=0)
		grate_std=np.nanstd(grate_list,axis=0)
		
		"figs-lines"
		for run in range(num_dup):
			ax1.plot(size,fg_np[run] , linestyle = '-',linewidth=0.3,color=color_list[stra_cluster+1],alpha=alp_small)	
			ax1.scatter(size,fg_np[run] ,s=1,marker="o",label=label_list[stra_cluster+1],color=color_list[stra_cluster+1],alpha=alp_small)
		
			ax2.plot(size,r_np[run] , linestyle = '-',linewidth=0.3,color=color_list[stra_cluster+1],alpha=alp_small)	
			ax2.scatter(size,r_np[run] ,s=1,marker="o",label=label_list[stra_cluster+1],color=color_list[stra_cluster+1],alpha=alp_small)
		
			ax3.scatter(size[stra_cluster+1],grate_list[run] ,s=1,marker="o",label=label_list[stra_cluster+1],color=color_list[stra_cluster+1],alpha=1)
		"shades"
		ax1.scatter(size,fg_mean ,s=50,marker="o",label=label_list[stra_cluster+1],color=color_list[stra_cluster+1],alpha=1)
		ax1.plot(size,fg_mean , linestyle = '-',linewidth=3,color=color_list[stra_cluster+1],alpha=1)		
		ax1.fill_between(size,fg_mean-fg_std,fg_mean+fg_std,label=label_list[stra_cluster+1],color=color_list[stra_cluster+1],alpha=0.1)
	
		ax2.scatter(size,r_mean ,s=50,marker="o",label=label_list[stra_cluster+1],color=color_list[stra_cluster+1],alpha=1)
		ax2.plot(size,r_mean , linestyle = '-',linewidth=3,color=color_list[stra_cluster+1],alpha=1)		
		ax2.fill_between(size,r_mean-r_std,r_mean+r_std,label=label_list[stra_cluster+1],color=color_list[stra_cluster+1],alpha=0.1)

		ax3.errorbar(size[stra_cluster+1],grate_mean,grate_std,linestyle='-',linewidth=1,marker="h",markersize=12,label=label_list[stra_cluster+1],
			   color=color_list[stra_cluster+1],alpha=0.8 )
#		
		ax3.set_xlim(-0.6,4.6)
#		ax3.set_ylim(0.65,0.7)

	plt.show()
	fig.savefig('./fg_r_grate_b%s_c%s.pdf'%(b,c), bbox_inches = 'tight')   # save figures


