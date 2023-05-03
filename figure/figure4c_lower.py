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
			max_dup[dup]=np.loadtxt('../data/data_fig4acd/%s_%s_%s_%s.txt'%(dup,n,alpha_cluster,beta_cluster))
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
		
		num_nd_rd=class_percent_nd[alpha_cluster][beta_cluster]+class_percent_rd[alpha_cluster][beta_cluster]
		id_num=dup_num-num_nd_rd
		class_percent_isd[alpha_cluster][beta_cluster]=np.count_nonzero(all_classify_narry[:,alpha_cluster,beta_cluster]==2)/id_num #isd
		class_percent_igd[alpha_cluster][beta_cluster]=np.count_nonzero(all_classify_narry[:,alpha_cluster,beta_cluster]==1)/id_num #isd
		class_percent_id[alpha_cluster][beta_cluster]=np.count_nonzero(all_classify_narry[:,alpha_cluster,beta_cluster]==3)/id_num #isd

class_percent_isd_narry=np.array([np.array(i) for i in class_percent_isd])
class_percent_igd_narry=np.array([np.array(i) for i in class_percent_igd])
class_percent_id_narry=np.array([np.array(i) for i in class_percent_id])

#-------------------------------------------------------------------------------
"draw figures"
#from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
# ----set fig frame--------
fig, ax = plt.subplots(1, 3, gridspec_kw={'width_ratios': [1,1,1]},figsize=(9, 3.2))
#-----colormap--------------
color2=['w','#efedf5','#dadaeb','#bcbddc','#9e9ac8','#807dba','#6a51a3','#54278f','#3f007d']		
color3=['w','#ffffe5','#fff7bc','#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506']
color4=['w','#fff5f0','#fee0d2','#fcbba1','#fc9272','#fb6a4a','#ef3b2c','#cb181d','#a50f15','#67000d']

cmap2 = LinearSegmentedColormap.from_list('mycmap', color2)
cmap3 = LinearSegmentedColormap.from_list('mycmap', color3)
cmap4 = LinearSegmentedColormap.from_list('mycmap', color4)

norm1 = plt.Normalize(0, 1)
norm2 = plt.Normalize(0, 0.3)

#----figure------------------
im2 = ax[0].imshow(class_percent_igd_narry, interpolation=None, origin='lower',cmap=cmap2, norm=norm1)
im3 = ax[1].imshow(class_percent_isd_narry, interpolation=None, origin='lower',cmap=cmap3, norm=norm1)
im4 = ax[2].imshow(class_percent_id_narry, interpolation=None, origin='lower',cmap=cmap4, norm=norm1)

# title
ax[0].set_title(r'$IGD$',fontsize=12)
ax[1].set_title(r'$ISD$',fontsize=12)
ax[2].set_title(r'$IGSD$',fontsize=12)

# x and y label
ax[0].set_xlabel(r'Cell differentiation benefit $b$',fontsize=12)
ax[0].set_ylabel(r'Cell differentiation cost $c$',fontsize=12)
ax[1].set_xlabel(r'Cell differentiation benefit $b$',fontsize=12)
ax[2].set_xlabel(r'Cell differentiation benefit $b$',fontsize=12)

# xy ticks
ax[0].tick_params(direction='inout', length=3, width=1, colors='k')
ax[1].tick_params(direction='inout', length=3, width=1, colors='k')
ax[2].tick_params(direction='inout', length=3, width=1, colors='k')

# artifical x and y ticks
test0=np.linspace(0,grid_num-1,4,endpoint=True)
ax[0].set_xticks( test0, minor=False)
test=np.linspace(alpha_expo_range[0],alpha_expo_range[1],4,endpoint=True)
x_label=[r'$10^{'+str(int(i))+'}$' for i in test]
ax[0].xaxis.set_major_formatter(mpl.ticker.FixedFormatter(x_label))

ax[0].set_yticks( test0, minor=False)
y_test=np.linspace(beta_expo_range[0],beta_expo_range[1],4,endpoint=True)
y_label=[r'$10^{'+str(int(i))+'}$' for i in y_test]
ax[0].yaxis.set_major_formatter(mpl.ticker.FixedFormatter(y_label))

# fig 2 x
ax[1].set_xticks( test0, minor=False)
ax[1].xaxis.set_major_formatter(mpl.ticker.FixedFormatter(x_label))
ax[2].set_xticks( test0, minor=False)
ax[2].xaxis.set_major_formatter(mpl.ticker.FixedFormatter(x_label))
ax[1].set_yticks([] )
ax[2].set_yticks([] )

# color bar
cbar_ax0 = fig.add_axes([0.358, 0.3, 0.0085, 0.5])
cbar0=fig.colorbar(im2, cax=cbar_ax0, orientation="vertical",norm=norm1, boundaries=None)
cbar0.ax.tick_params(labelsize=8,length=2,direction='in')
cbar0.outline.set_visible(False)

cbar_ax1 = fig.add_axes([0.668, 0.3, 0.0082, 0.5])
cbar1=fig.colorbar(im3, cax=cbar_ax1, orientation="vertical",norm=norm1, boundaries=None)
cbar1.ax.tick_params(labelsize=8,length=2,direction='in')
cbar1.outline.set_visible(False)

cbar_ax2 = fig.add_axes([0.98, 0.3, 0.008, 0.5])
cbar2=fig.colorbar(im4, cax=cbar_ax2, orientation="vertical",norm=norm1, boundaries=None)
cbar2.ax.set_ylabel(r'Fraction of strategy', rotation=90,fontsize=10)
cbar2.ax.tick_params(labelsize=8,length=2,direction='in')
cbar2.outline.set_visible(False)

plt.show()
#fig.savefig('./fig4c_lower_3minor_dup%s_n%s.pdf'%(num_dup,n), bbox_inches='tight')   # save figures

