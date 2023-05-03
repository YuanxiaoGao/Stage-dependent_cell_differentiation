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
'''read data '''
static=np.loadtxt('./static_m.txt')   # 0=ND, 1=RD, 2=ISD
dy_nd=np.loadtxt('./dynamic_nd.txt')
dy_rd=np.loadtxt('./dynamic_rd.txt')
dy_igd=np.loadtxt('./dynamic_igd.txt')
dy_isd=np.loadtxt('./dynamic_isd.txt')
dy_igsd=np.loadtxt('./dynamic_igsd.txt')

grid_num=41
result=[static,dy_nd,dy_rd,dy_igd,dy_isd,dy_igsd]
       #  0,    1   ,  2      3      4      5
np.shape(result)


# ND-> RD and ID        0 RD invade ND; 1 ID invade ND   
ND_result=[[np.nan for i in range(grid_num)] for i in range(grid_num)]
np.shape(ND_result)

for alpha_cluster in range(0,grid_num):
    for beta_cluster in range(0,grid_num):
        s=result[0][alpha_cluster][beta_cluster]           # static strategy
        dnd=result[1][alpha_cluster][beta_cluster]
        drd=result[2][alpha_cluster][beta_cluster]
        did=result[3][alpha_cluster][beta_cluster]+result[4][alpha_cluster][beta_cluster]+result[5][alpha_cluster][beta_cluster]
        
        # invade ND
        if s==0 and did>=dnd and did>=drd:
            ND_result[alpha_cluster][beta_cluster]=0    
        elif s==0 and drd>=dnd and drd>=did:
            ND_result[alpha_cluster][beta_cluster]=1   
        elif s==0 and dnd>=did and dnd>=drd:
            ND_result[alpha_cluster][beta_cluster]=-3           # no invasion and stable ND AREA   
        # INVADE rd
        if s==1 and did>=dnd and did>=drd:
            ND_result[alpha_cluster][beta_cluster]=2    
        elif s==1 and dnd>=drd and dnd>=did:
            ND_result[alpha_cluster][beta_cluster]=3
        elif s==1 and drd>=dnd and drd>=did:
            ND_result[alpha_cluster][beta_cluster]=-2         # no invasion and stable RD AREA
        # invade ID
        if s==2 and dnd>=drd and dnd>=did:
            ND_result[alpha_cluster][beta_cluster]=4    
        elif s==2 and drd>=dnd and drd>=did:
            ND_result[alpha_cluster][beta_cluster]=5
        elif s==2 and did>=dnd and did>=drd:
            ND_result[alpha_cluster][beta_cluster]=-1
ND_result_arr=np.array(ND_result) 
# -3  -2    -1     0     1X         2         3X         4X       5      X NON-EXIST
# ND  RD   ISD  ND->id  ND->rd    RD->ID   RD->ND   ID->ND    ID->RD
#-------------------------------------------------------------------------------
'''read data '''
'''Parameters'''
alpha=[-1,2]
beta=[-1,2]

grid_num=41
alpha_expo_range=np.array([alpha[0],alpha[1]])                   # exponential distribution of alpha
grid_alpha=np.linspace(alpha_expo_range[0],alpha_expo_range[1],num=grid_num,endpoint=True) # split log(a) into grid_num points
alpha_list=10**grid_alpha

beta_expo_range=np.array([beta[0],beta[1]])
grid_beta=np.linspace(beta_expo_range[0],beta_expo_range[1],num=grid_num,endpoint=True) # split log(a) into grid_num points
beta_list=10**grid_beta

"draw figures"
#from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch

# ----set fig frame--------
fig, ax = plt.subplots(1, 1,figsize=(5, 5))
# legend
legend_elements = [Patch(facecolor='#238b45', edgecolor='#238b45',label=r'$\lambda_{ID}>\lambda_{ND}^{S}$'),
                   Patch(facecolor='#2171b5', edgecolor='#2171b5',label=r'$\lambda_{ID}>\lambda_{RD}^{S}$'),
                   Patch(facecolor='#8B4513', edgecolor='#8B4513',label=r'$\lambda_{RD}>\lambda_{ID}^{S}$')]

#-----colormap--------------
color0=['#DCDCDC','#A9A9A9','#000000', '#238b45','grey',"#2171b5","grey","grey","#8B4513"]
cmap0 = LinearSegmentedColormap.from_list('mycmap', color0)
norm1 = plt.Normalize(-3, 5)
#----figure------------------
im0 = ax.imshow(ND_result_arr, interpolation=None, origin='lower',cmap=cmap0, norm=norm1,vmin=-3)
cmap0.set_bad('grey') 

# x and y label
ax.set_xlabel(r'Cell differentiation benefit $b$',fontsize=14)
ax.set_ylabel(r'Cell differentiation cost $c$',fontsize=14)
#
# artifical x and y ticks
ax.tick_params(direction='in', length=5, width=1, colors='k')
ax.locator_params(axis='x', nbins=4)
ax.locator_params(axis='y', nbins=4)

# fig 1 xy axes
test0=np.linspace(0,grid_num-1,4,endpoint=True)
ax.set_xticks( test0, minor=False)
test=np.linspace(alpha_expo_range[0],alpha_expo_range[1],4,endpoint=True)
x_label=[r'$10^{'+str(int(i))+'}$' for i in test]
ax.xaxis.set_major_formatter(mpl.ticker.FixedFormatter(x_label))

ax.set_yticks( test0, minor=False)
y_test=np.linspace(beta_expo_range[0],beta_expo_range[1],4,endpoint=True)
y_label=[r'$10^{'+str(int(i))+'}$' for i in y_test]
ax.yaxis.set_major_formatter(mpl.ticker.FixedFormatter(y_label))

# add text strategy
ax.text(12,25,r'$ND$',fontsize=15,color='k')
ax.text(22,4,r'$RD$',fontsize=15,color='k')
ax.text(26,22,r'$ID$',fontsize=15,color='lightgrey')

ax.legend(handles=legend_elements, loc='upper left',fontsize=11.5)

plt.show()
#fig.savefig('./fig4a_invasion00.pdf', bbox_inches='tight')   # save figures

