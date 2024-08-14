# -*- coding: utf-8 -*-
"""
Created on Wed May 22 13:52:26 2024

@author: hli490
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Setting working directory
path = r'C:\Users\hli490\OneDrive - Vrije Universiteit Amsterdam\Desktop\Paper2'

df_diff_ead = pd.read_csv(os.path.join(path,'Impact/Coastline/diff_ead_estimates_new.csv'))
df_diff_rp200 = pd.read_csv(os.path.join(path,'Impact/Coastline/diff_rp200_estimates.csv'))
df_co = pd.read_csv(os.path.join(path,'Coastline/countries-by-coastline-2024.csv'),encoding='windows-1252')
df_diff = pd.merge(df_diff_ead,df_diff_rp200,on='FID')

df = pd.merge(df_diff, df_co, on="FID")
df['Coastline'] = df['Coastline'].astype(int)
df.to_csv(os.path.join(path,'Impact/Coastline/effect_coastline.csv'),index=False)
df = df.dropna()
## plotting
fig=plt.figure(figsize=(13,5))
plt.rcParams.update({'font.size':10})
gs=gridspec.GridSpec(1,2,wspace=0.15,hspace=0.1)
ax=[[] for _ in range(2)] 
c1 = '#0571b0' # blue
c2 = '#ca0020' # red

# EAD
ax[0] = plt.subplot(gs[0,0])
x = df['Coastline'].astype(int)
y = df['ead_nopro_mean']
ax[0].scatter(x,y,marker='o',edgecolor=c1,facecolor='none',label='Without protection')

x = df['Coastline'].astype(int)
y = df['ead_pro_mean']
ax[0].scatter(x,y,marker='x',color=c2,label='With protection')

x = df['Coastline'].astype(int)
y = df['ead_nopro_mean']
z = np.polyfit(np.log10(x), y, 1)
t = np.log10(x)*z[0]+z[1]

ax[0].plot(x,t,color=c1,linewidth=1.2,label='Trendline without protection')

x = df['Coastline'].astype(int)
y = df['ead_pro_mean']
z = np.polyfit(np.log10(x), y, 1)
t = np.log10(x)*z[0]+z[1]
ax[0].plot(x,t,color=c2,linewidth=1.2,label='Trendline with protection')
ax[0].grid(color='gray',linestyle = 'solid', alpha = 0.2) 

ax[0].set_ylim([-80,80])
ax[0].set_xscale('log')
ax[0].set_title('a) Expected annual damage (EAD)',fontsize=12)
ax[0].legend(loc='lower left')


### RP200
ax[1] = plt.subplot(gs[0,1])
x = df['Coastline'].astype(int)
y = df['rp200_nopro_mean']
ax[1].scatter(x,y,marker='o',edgecolor=c1,facecolor='none')
z = np.polyfit(np.log10(x), y, 1)
t = np.log10(x)*z[0]+z[1]
ax[1].plot(x,t,linewidth=1.2,color=c1)

x = df['Coastline'].astype(int)
y = df['rp200_pro_mean']
x = x[y<=100]
y = y[y<=100]
ax[1].scatter(x,y,marker='x',color=c2)
z = np.polyfit(np.log10(x), y, 1)
t = np.log10(x)*z[0]+z[1]
ax[1].plot(x,t,linewidth=1.2,color=c2)
ax[1].grid(color='gray',linestyle = 'solid', alpha = 0.2) 

ax[1].set_ylim([-80,80])
ax[1].set_xscale('log')
ax[1].set_title('b) 1 in 200-year damage',fontsize=12)

fig.text(0.5, 0.02, 'Coastline length (km)',fontsize=12, ha='center')
fig.text(0.08, 0.5, 'Relative difference (%)',fontsize=12, va='center', rotation='vertical')
outputdir = os.path.join(path,'Figures')
plt.savefig(os.path.join(outputdir,'coastline_diff_new.png'),format='png',dpi=400,bbox_inches='tight')


# df.to_csv(r'C:\Users\hli490\Desktop\Paper2\Coastline\diff_coastline.csv',index=False)
