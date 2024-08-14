# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 12:58:46 2024

@author: hli490
"""

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cmplt
import cartopy.crs as ccrs
import re
import mapclassify
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.axes_divider import make_axes_area_auto_adjustable

# set the colormap and centre the colorbar
class MidpointNormalize(mpl.colors.Normalize):
    """Normalise the colorbar."""
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        mpl.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))
    
# Path and folder location
path = r'C:\Users\hli490\OneDrive - Vrije Universiteit Amsterdam\Desktop\Paper2'

# df_rp200['FID'] = df_rp200['iso']

rp200_nopro = pd.read_csv(os.path.join(path,'Impact/new_rp200_diff_nopro.csv'))
rp200_nopro['diff'] = (rp200_nopro['Event_based_mean']-rp200_nopro['RP_based'])/rp200_nopro['RP_based']*100
rp200_nopro['count'] = np.nan # showing whether the difference is statistically insignificant
for i in range(len(rp200_nopro)):
    if np.isnan(rp200_nopro['RP_based'][i]): continue
    if rp200_nopro['RP_based'][i] >= rp200_nopro['Event_based_low'][i] and rp200_nopro['RP_based'][i] <= rp200_nopro['Event_based_high'][i]:
        rp200_nopro.iloc[i,7]=1
    else:
        rp200_nopro.iloc[i,7]=0

rp200_pro = pd.read_csv(os.path.join(path,'Impact/new_rp200_diff_pro.csv'))
rp200_pro['diff'] = (rp200_pro['Event_based_mean']-rp200_pro['RP_based'])/rp200_pro['RP_based']*100
rp200_pro['count'] = np.nan# showing whether the difference is statistically insignificant
for i in range(len(rp200_pro)):
    if np.isnan(rp200_pro['diff'][i]): continue
    if rp200_pro['RP_based'][i] >= rp200_pro['Event_based_low'][i] and rp200_pro['RP_based'][i] <= rp200_pro['Event_based_high'][i]:
        rp200_pro.iloc[i,7]=1
    else:
        rp200_pro.iloc[i,7]=0


df_diff = pd.DataFrame(data=dict(FID=rp200_nopro['FID'],
                                 ISO=rp200_nopro['ISO'],
                                 diff_nopro=rp200_nopro['diff'],
                                 count_nopro=rp200_nopro['count'],
                                 diff_pro=rp200_pro['diff'],
                                 count_pro=rp200_pro['count'],
                                 change=rp200_pro['diff']-rp200_nopro['diff']))

#### Plotting
# read the shapefile using geopandas
shpfilename = os.path.join(path,'Geogunit/geogunit_0.shp')
shp = gpd.read_file(shpfilename,engine="pyogrio")
shp = (shp.sort_values(by='ISO')).reset_index(drop=True)
shp['FID'] = shp.index
df = shp.merge(df_diff,how='left',on='FID')

# world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))

#### EAD

### Difference
# set a variable that will call whatever column we want to visualise on the map
# variable = 'EAD_diff_re'
crg = ccrs.PlateCarree()
fig=plt.figure(figsize=(10,8)) # width, height
gs=gridspec.GridSpec(2,1,wspace=0.15,hspace=0.15)
# gs=gridspec.GridSpec(2,1,height_ratios=[1.5,1.5,1,1])
# gs=gridspec.GridSpec(4,2,wspace=0.3,hspace=0.25)
ax=[[] for _ in range(0,6)] 
plt.rcParams.update({'font.size':10})
plt.rcParams['hatch.linewidth'] = 0.5
plt.rcParams['hatch.color'] = 'green'


# for idx,gs1,gs2,label,var in zip(range(0,2),[0,1],[0,0],['a)','b)'],['EAD_diff_nopro','EAD_diff_pro']):
for idx,label,var,var_hatch,title in zip(range(0,2),['a)','b)'],['diff_nopro','diff_pro'],['count_nopro','count_pro'],
                                                 ['RP200 damage difference without flood protection','RP200 damage difference with flood protection']):
    # print(1)
    
    if idx == 0:
        ax[idx] = plt.subplot(gs[0,:], projection=crg)
    else:
        ax[idx] = plt.subplot(gs[1,:], projection=crg)
    ax[idx].text(-0.05,1.03,label,fontsize=11,transform=ax[idx].transAxes) 
    
    cmap = plt.get_cmap('PiYG_r')
    #         # cmap1 = mpl.colormaps.get_cmap("RdBu_r")
    bins = np.arange(-80,45,5)
        
    norm = MidpointNormalize(min(bins), max(bins), 0.)
    # df.plot(ax=ax[idx],column='RP200_diff_pro', cmap=cmap1, norm=norm, legend='True', scheme='user_defined',classification_kwds={'bins':bins},
    #             missing_kwds={
    #         "color": "gainsboro",
    #         "edgecolor": "k",
    #         "label": "No Flood Damages",
    #     },
    #     linewidth=0.3, edgecolor='k')
    df.plot(ax=ax[idx],column=var, cmap=cmap, norm=norm, linewidth=0.3, edgecolor='k',missing_kwds={
            "color": "darkgray",
            "edgecolor": "k",
            "label": "No Flood Damages",
        })
    
    df.loc[df[var_hatch]==1].plot(ax=ax[idx],column=var, cmap=cmap, norm=norm, linewidth=0.3, edgecolor='k',hatch='/////')
    ax[idx].set_extent([-180, 180, -60, 90], crs=crg)
    ax[idx].set_title(title,fontsize=11)
    # cax1= fig.add_axes([0.12,0.3,0.2,0.02]) #[left, bottom, width, height]
    
    ax[idx].spines['top'].set_visible(False)
    ax[idx].spines['right'].set_visible(False)
    ax[idx].spines['bottom'].set_visible(False)
    ax[idx].spines['left'].set_visible(False)
    # ax[idx].axis("off")

    
# cax1= fig.add_axes([0.93,0.517,0.015,0.28]) #[left, bottom, width, height]
cax1= fig.add_axes([0.38,0.065,0.28,0.015]) #[left, bottom, width, height]
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
cbar= plt.colorbar(sm,cax=cax1,cmap=cmap,spacing='uniform',
                    orientation='horizontal',extend='both',boundaries=bins,
                    ticks=[-80,-60,-40,-20,0,20,40]) #relative difference
cbar.set_label('Relative difference (%)',size=11)
cbar.ax.set_xticklabels(['-80',  '-60', '-40','-20',   '0',  '20',  '40']) 

filename='rp200_diff_map_new_try.png'
outputdir = os.path.join(path,'Figures\EAD')
plt.savefig(os.path.join(outputdir,filename),format='png',dpi=300,bbox_inches='tight') 

df_new = df_diff.drop(df_diff[df_diff['diff_nopro']>200].index)
df_new = df_diff.drop(df_diff[df_diff['diff_pro']>200].index)
abs(df_new['diff_nopro']).mean()
abs(df_new['diff_nopro']).std()

abs(df_new['diff_pro']).mean()
abs(df_new['diff_pro']).std()

sum(df_diff['diff_nopro']<=0)


nopro = abs(df_diff['diff_nopro'])
pro = abs(df_diff['diff_pro'])
t = pro-nopro
sum(t>0)
# set a variable that will call whatever column we want to visualise on the map
# variable = 'RP200_diff_re'

# # read the shapefile using geopandas
# shpfilename = r'C:\Users\hli490\Desktop\Paper2\Geogunit\geogunit_0.shp'
# shp = gpd.read_file(shpfilename,engine="pyogrio")
# shp = (shp.sort_values(by='ISO')).reset_index(drop=True)
# shp['FID'] = shp.index
# # shp = shp.loc[ind]
# df = shp.merge(df_rp200,how='left',on='FID')

# # create figure and axes for Matplotlib

# crg = ccrs.PlateCarree()
# fig=plt.figure(figsize=(12,6))
# ax=plt.subplot(projection=crg)
# cmap = mpl.colormaps.get_cmap("Greens")
# # bins = np.arange(0.0,1.1,0.1)
# bins = np.arange(0.0,2.2,0.2)
# norm = mpl.colors.BoundaryNorm(bins, cmap.N)
# df.plot(ax=ax,column=variable, cmap=cmap,scheme='user_defined',classification_kwds={'bins':bins},
#             missing_kwds={
#         "color": "gainsboro",
#         "edgecolor": "black",
#         "label": "No Data",
#     },
#     linewidth=0.3, edgecolor='k')

# # df.plot(ax=ax,column=variable, cmap=cmap,scheme='user_defined',classification_kwds={'bins':bins},
# #             missing_kwds={
# #         "color": "gainsboro",
# #         "edgecolor": "k",
# #         "hatch": "/",
# #         "label": "No Data",
# #     },
# #     linewidth=1, edgecolor='k')

# ax.set_extent([-180, 180, -60, 90], crs=crg)
# # cax1= fig.add_axes([0.12,0.3,0.2,0.02]) #[left, bottom, width, height]
# cax1= fig.add_axes([0.46,0.2,0.2,0.02]) #[left, bottom, width, height]
# sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
# cbar= plt.colorbar(sm,ax=ax,cax=cax1,cmap=cmap,spacing='uniform',
#                     orientation='horizontal',extend='max',boundaries=bins,
#                     ticks=[0,0.4,0.8,1.2,1.6,2.0]) #relative difference
# cbar.set_label('Relative difference (-)',size=10)
# cbar.ax.tick_params(labelsize=10)
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# ax.spines['bottom'].set_visible(False)
# ax.spines['left'].set_visible(False)
# ax.axis("off")

# filename='diff_RP200.png'
# outputdir = r'C:\Users\hli490\Desktop\Paper2\Figures\EAD'
# plt.savefig(os.path.join(outputdir,filename),format='png',dpi=600,bbox_inches='tight') 

# plt.hist(df_rp200['RP200_diff_re'])

# ead_iso = np.array(df_rp200.loc[df_rp200['EAD-GTSM']>0,'FID'])
# rp200_iso = np.array(df_rp200.loc[df_rp200['RP200-GTSM']>0,'FID'])
# np.setdiff1d(ead_iso,rp200_iso)
# from matplotlib.ticker import PercentFormatter
# fig=plt.figure(figsize=(10,6))
# bins=np.arange(0.0,2.6,0.2)
# ax=plt.subplot()
# ax.hist([df_rp200['EAD_diff_re'],df_rp200['RP200_diff_re']], bins=bins,label=['EAD','RP200'])
