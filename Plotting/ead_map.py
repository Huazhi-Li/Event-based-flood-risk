# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 16:46:27 2023

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
from scipy import stats, interpolate
import mapclassify
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
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


# EAD
# df_ead = pd.read_csv(r'C:\Users\hli490\Desktop\Paper2\Impact\EAD_diff.csv')
# df_ead_pro = pd.read_csv(r'C:\Users\hli490\Desktop\Paper2\Impact\EAD_diff_pro.csv')
# df_ead['diff'] = (df_ead['Complete_dependence']-df_ead['Annual_max'])/df_ead['Annual_max']*100
# df_ead_pro['diff'] = (df_ead_pro['Complete_dependence']-df_ead_pro['Annual_max'])/df_ead_pro['Annual_max']*100

# # RP200 damage
# df_rp200 = pd.read_csv(r'C:\Users\hli490\Desktop\Paper2\Impact\rp200_diff.csv')
# df_rp200['diff'] = (df_rp200['Complete_dependence']-df_rp200['Annual_max'])/df_rp200['Annual_max']*100
# df_rp200_pro = pd.read_csv(r'C:\Users\hli490\Desktop\Paper2\Impact\rp200_diff_pro.csv')
# df_rp200_pro['diff'] = (df_rp200_pro['Complete_dependence']-df_rp200_pro['Annual_max'])/df_rp200_pro['Annual_max']*100

# difference dataframe
# df_diff = pd.DataFrame(data=dict(FID=df_ead['FID'],
#                                   ISO=df_ead['ISO'],
#                                   EAD_diff_nopro=df_ead['diff'],
#                                   EAD_diff_pro=df_ead_pro['diff'],
#                                   RP200_diff_nopro=df_rp200['diff'],
#                                   RP200_diff_pro=df_rp200_pro['diff']))
# df_diff.to_csv(r'C:\Users\hli490\Desktop\Paper2\Impact\diff_in_estimates.csv',index=False)

ead_nopro = pd.read_csv(os.path.join(path,'Impact/new_EAD_diff_nopro.csv'))
ead_nopro['diff'] = (ead_nopro['Event_based_mean']-ead_nopro['RP_based'])/ead_nopro['RP_based']*100
ead_nopro['count'] = np.nan # showing whether the difference is statistically insignificant
for i in range(len(ead_nopro)):
    if np.isnan(ead_nopro['RP_based'][i]): continue
    if ead_nopro['RP_based'][i] >= ead_nopro['Event_based_low'][i] and ead_nopro['RP_based'][i] <= ead_nopro['Event_based_high'][i]:
        ead_nopro.iloc[i,7]=1
    else:
        ead_nopro.iloc[i,7]=0

ead_pro = pd.read_csv(os.path.join(path,'Impact/new_EAD_diff_pro.csv'))
ead_pro['diff'] = (ead_pro['Event_based_mean']-ead_pro['RP_based'])/ead_pro['RP_based']*100
ead_pro['count'] = np.nan# showing whether the difference is statistically insignificant
for i in range(len(ead_pro)):
    if np.isnan(ead_pro['diff'][i]): continue
    if ead_pro['RP_based'][i] >= ead_pro['Event_based_low'][i] and ead_pro['RP_based'][i] <= ead_pro['Event_based_high'][i]:
        ead_pro.iloc[i,7]=1
    else:
        ead_pro.iloc[i,7]=0

df_diff = pd.DataFrame(data=dict(FID=ead_nopro['FID'],
                                 ISO=ead_nopro['ISO'],
                                 diff_nopro=ead_nopro['diff'],
                                 count_nopro=ead_nopro['count'],
                                 diff_pro=ead_pro['diff'],
                                 count_pro=ead_pro['count'],
                                 change=ead_pro['diff']-ead_nopro['diff']))


# df_diff = pd.read_csv(os.path.join(path,'Impact/new_EAD_diff_nopro.csv'))
#### Plotting
# read the shapefile using geopandas
shpfilename = os.path.join(path,'Geogunit/geogunit_0.shp')
shp = gpd.read_file(shpfilename,engine="pyogrio")
shp = (shp.sort_values(by='ISO')).reset_index(drop=True)
shp['FID'] = shp.index
df = shp.merge(df_diff,how='left',on='FID')



### Plotting
# set a variable that will call whatever column we want to visualise on the map
# variable = 'EAD_diff_re'

### horizontal
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
                                                 ['EAD difference without flood protection','EAD difference with flood protection']):
    # print(1)
    
    if idx == 0:
        ax[idx] = plt.subplot(gs[0,:], projection=crg)
    else:
        ax[idx] = plt.subplot(gs[1,:], projection=crg)
    ax[idx].text(-0.05,1.03,label,fontsize=11,transform=ax[idx].transAxes) 
    
    cmap = mpl.colormaps.get_cmap("RdBu_r")
    bins = np.arange(-40,45,5)
        
    norm = mpl.colors.BoundaryNorm(bins, cmap.N)
    # df.plot(ax=ax[idx],column='RP200_diff_pro', cmap=cmap1, norm=norm, legend='True', scheme='user_defined',classification_kwds={'bins':bins},
    #             missing_kwds={
    #         "color": "gainsboro",
    #         "edgecolor": "k",
    #         "label": "No Flood Damages",
    #     },
    #     linewidth=0.3, edgecolor='k')
    df.plot(ax=ax[idx],column=var, cmap=cmap, norm=norm, linewidth=0.3, edgecolor='k',missing_kwds={
            "color": "darkgrey",
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
                    ticks=[-40,-20,0,20,40]) #relative difference
cbar.set_label('Relative difference (%)',size=11)
cbar.ax.set_xticklabels(['-40',  '-20',   '0',  '20',  '40']) 

filename='EAD_diff_map_new_try.png'
outputdir = os.path.join(path,'Figures\EAD')
plt.savefig(os.path.join(outputdir,filename),format='png',dpi=300,bbox_inches='tight') 

abs(df_diff['diff_nopro']).mean()
abs(df_diff['diff_nopro']).std()

abs(df_diff['diff_pro']).mean()
abs(df_diff['diff_pro']).std()
# cax1= fig.add_axes([0.32,0.53,0.35,0.02]) #[left, bottom, width, height]
# sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
# cbar= plt.colorbar(sm,ax=ax[idx],cax=cax1,cmap=cmap,spacing='uniform',
#                     orientation='horizontal',extend='both',boundaries=bins,
#                     ticks=[-40,-20,0,20,40]) #relative difference
# cbar.set_label('EAD difference (%)')
# cbar.ax.set_xticklabels(['-40', '-20', '0',  '20', '40']) 


## add risk curves
# def get_return_periods(x, yr=10000):
#         # weibull plotting formula a=0, b=1
#         a = 0.0
#         b = 1.0
#         ranks = (len(x) + 1) - stats.rankdata(x, method="ordinal")
#         lam = len(x)/yr
#         freq = ((ranks - a) / (len(x) + b)) * lam
#         rps = 1 / freq
#         # return ranks, rps
#         return rps
    

# Getting iso information
# geog0 = pd.read_csv(os.path.join(path,'Geogunit/geogunit_0_list.csv'),encoding = 'unicode_escape')
# df_iso = pd.read_csv(os.path.join(path,'Geogunit/Station/station_info.csv'))
# df_events = pd.read_csv(os.path.join(path,'Eventset/num_event_99th.csv'))

# # color
# c1 = '#ca0020' # red RP-based no protection
# c2 = '#0571b0' # blue Event-based no protection
# c3 = 'purple' # purple RP-based  protection
# c4 = '#35afd7' # light blue Event-based protection
# # c4 = '#3fa671' # green Event-based protection
# isos= [236,46,87,115]
# countries=['United States','Republic of the Congo','Guinea-Bissau','Japan']
# scale = 1e6

# # for idx,gs1,gs2,label,iso,country in zip(range(0,8),[0,0,1,1,2,2,3,3],[0,1,0,1,0,1,0,1],['a)','b)','c)','d)','e)','f)','g)','h)'],isos,countries):
# for idx,gs1,gs2,label,iso,country in zip(range(2,6),[2,2,3,3],[0,1,0,1],['c)','d)','e)','f)'],isos,countries):
    
#     iso_in = df_iso[df_iso['iso']==iso].reset_index(drop=True)
#     basins = np.unique(iso_in['basin'])
#     dmg_all = pd.read_csv(os.path.join(path, 'Impact/Glofris/glofris_Built-up_Area_geogunit0.csv'))
    
#     file = os.path.join(path, 'Impact/complete_dependence_Built-up_Area_geogunit0.csv')
#     if not os.path.isfile(file) or dmg_all['1000'][iso]==0: continue
    
#     # No protection
#     dmg_all_nopro = pd.read_csv(os.path.join(path, 'Impact/complete_dependence_Built-up_Area_geogunit0.csv'))
#     annual_max_all_nopro = pd.read_csv(os.path.join(path,'Impact/Annual_damage/Geogunit0/nopro/annual_max_dmg_unit{:03d}.csv'.format(iso)))
#     annual_max_nopro = annual_max_all_nopro[annual_max_all_nopro.columns[0:3]]
#     dmg_unit_nopro = np.concatenate([[0,0],np.array(dmg_all_nopro.loc[dmg_all_nopro['FID']==iso])[0][2:]])
#     rp = np.array([1,2,5,10,25,50,100,250,500,1000]).astype(int)
#     dmg_cd_nopro = pd.DataFrame(data=dict(rp=rp,dmg=dmg_unit_nopro))
    
#     # With protection
#     dmg_all_pro = pd.read_csv(os.path.join(path, 'Impact/complete_dependence_Built-up_Area_geogunit0_pro.csv'))
#     annual_max_all_pro = pd.read_csv(os.path.join(path,'Impact/Annual_damage/Geogunit0/pro/annual_max_dmg_unit{:03d}.csv'.format(iso)))
#     annual_max_pro = annual_max_all_pro[annual_max_all_pro.columns[0:3]]
#     dmg_cd_pro = pd.DataFrame(data=dict(rp= np.arange(1,1001),dmg=np.array(dmg_all_pro.loc[dmg_all_pro['FID']==iso])[0][1:]))
    
#     rp_lowest_dmg = np.where(dmg_cd_pro['dmg']>0)[0][0].min()
#     if not rp_lowest_dmg in rp.astype(int):
#         rp_nearest_idx = np.argmin(abs(rp_lowest_dmg-rp))
#         if rp[rp_nearest_idx]<rp_lowest_dmg: 
#             rp_insert = np.arange(rp_lowest_dmg, rp[rp_nearest_idx+1])
#         else:
#             rp_insert = np.arange(rp_lowest_dmg, rp[rp_nearest_idx])
#         rp_new = rp
#         rp_new = np.sort(np.concatenate([rp,rp_insert]))
#         dmg_cd_pro = dmg_cd_pro.loc[dmg_cd_pro['rp'].isin(rp_new),:]
      
#     else:
#         dmg_cd_pro = dmg_cd_pro.loc[dmg_cd_pro['rp'].isin(rp),:]
    
#     # Bootstrapping
#     rp_bootstrap = np.arange(1,1001)
#     yr = 1000
#     n = 100
#     # event_bootstrap = pd.DataFrame(index=range(n),columns=rp_bootstrap)
#     max_pro_bootstrap = pd.DataFrame(index=range(n),columns=rp_bootstrap)
#     max_nopro_bootstrap = pd.DataFrame(index=range(n),columns=rp_bootstrap)
    
#     for i in range(n):
    
#         # No protection
#         max_dmg_nopro = np.array(annual_max_nopro['dmg'].sample(n=int(len(annual_max_nopro)/(10000/yr)),replace=True))
#         max_rp_nopro = get_return_periods(max_dmg_nopro, yr=1000)
#         max_annual_nopro = pd.DataFrame(data=dict(dmg=max_dmg_nopro,rp=max_rp_nopro))
#         max_annual_nopro = max_annual_nopro.sort_values(by='rp')
#         max_nopro_bootstrap.iloc[i,:] = np.interp(rp_bootstrap,max_annual_nopro['rp'],max_annual_nopro['dmg'])
        
#         # With protection
#         max_dmg_pro = np.array(annual_max_pro['dmg'].sample(n=int(len(annual_max_pro)/(10000/yr)),replace=True))
#         max_rp_pro = get_return_periods(max_dmg_pro, yr=1000)
#         max_annual_pro = pd.DataFrame(data=dict(dmg=max_dmg_pro,rp=max_rp_pro))
#         max_annual_pro = max_annual_pro.sort_values(by='rp')
#         max_pro_bootstrap.iloc[i,:] = np.interp(rp_bootstrap,max_annual_pro['rp'],max_annual_pro['dmg'])
    
    
#     # Medians and confidence intervals
#     max_nopro_mean = max_nopro_bootstrap.mean().astype(float)
#     max_nopro_ci_lower = max_nopro_bootstrap.quantile(0.05, numeric_only=False).astype(float)
#     max_nopro_ci_upper = max_nopro_bootstrap.quantile(0.95, numeric_only=False).astype(float)
    
#     max_pro_mean = max_pro_bootstrap.mean().astype(float)
#     max_pro_ci_lower = max_pro_bootstrap.quantile(0.05, numeric_only=False).astype(float)
#     max_pro_ci_upper = max_pro_bootstrap.quantile(0.95, numeric_only=False).astype(float)
    
#     ### plotting

#     ax[idx] = plt.subplot(gs[gs1,gs2])
#     ax[idx].text(-0.05,1.03,label,fontsize=11,transform=ax[idx].transAxes) 
    
#     ax[idx].plot(dmg_cd_nopro['rp'],dmg_cd_nopro['dmg'], color = c1, linewidth=2, label='RP-based no protection')
#     ax[idx].plot(rp_bootstrap,max_nopro_mean,color=c2,linewidth=2,label='Event-based no protection')
#     ax[idx].fill_between(rp_bootstrap,max_nopro_ci_lower,max_nopro_ci_upper,facecolor=c2,alpha=0.3)
    
#     ax[idx].plot(dmg_cd_pro['rp'],dmg_cd_pro['dmg'], color = c3, alpha=0.75, linewidth=2, linestyle = "--",label='RP-based with protection')
#     ax[idx].plot(rp_bootstrap,max_pro_mean,color=c4,linewidth=2,linestyle = "--",label='Event-based with protection')
#     ax[idx].fill_between(rp_bootstrap,max_pro_ci_lower,max_pro_ci_upper,facecolor=c4,alpha=0.3)
    
#     # ax[idx].plot(dmg_cd_nopro['rp'],dmg_cd_nopro['dmg']/scale, color = c1, linewidth=2, label='RP-based without protection',path_effects=[pe.Stroke(linewidth=5, foreground='w'), pe.Normal()])
#     # ax[idx].plot(rp_bootstrap,max_nopro_mean/scale,color=c2,linewidth=2,label='Event-based without protection')
#     # ax[idx].fill_between(rp_bootstrap,max_nopro_ci_lower/scale,max_nopro_ci_upper/scale,facecolor=c2,alpha=0.3)
    
#     # ax[idx].plot(dmg_cd_pro['rp'],dmg_cd_pro['dmg']/scale, color = c3, alpha=0.75, linewidth=2, linestyle = "--",label='RP-based with protection',path_effects=[pe.Stroke(linewidth=5, foreground='w'), pe.Normal()])
#     # ax[idx].plot(rp_bootstrap,max_pro_mean/scale,color=c4,linewidth=2,linestyle = "--",label='Event-based with protection')
#     # ax[idx].fill_between(rp_bootstrap,max_pro_ci_lower/scale,max_pro_ci_upper/scale,facecolor=c4,alpha=0.3)

    
#     ax[idx].set_xlim([1, 1000])
#     ax[idx].set_ylim(ymin=0)
#     ax[idx].set_xscale('log')
#     ax[idx].grid(color='gray',linestyle = 'solid', alpha = 0.2)
#     ax[idx].set_xticks([1, 2, 5, 10, 25, 50, 100, 250, 500, 1000])
    
#     if idx==2:
#         ax[idx].legend(loc ='upper left',handlelength=2.2,prop={'size': 8})
    
#     # if idx==2:
#     #     ax[idx].set_ylabel('Economic damage (USD2005)')
#     # if gs2==0:
#     #     ax[idx].set_ylabel('Economic damage (USD2005)')
    
#     if gs1==2:
        
#         ax[idx].set_xticklabels([])
#         # ax[idx].set_xlabel('Return Period (year)')
    
#     # ax[idx].set_title(label+' '+country,fontsize=14)
#     ax[idx].set_title(country,fontsize=11)

# fig.text(0.5, 0.075, 'Return period (years)',fontsize=11, ha='center')
# fig.text(0.06, 0.26, 'Flood economic damage (USD2005)',fontsize=11,  va='center', rotation='vertical')



sum((df_diff['diff_nopro'] >= -10) & (df_diff['diff_nopro'] <= 10))


### Horizontal
# crg = ccrs.PlateCarree()
# fig=plt.figure(figsize=(10,8)) # width, height
# # gs=gridspec.GridSpec(1,2,wspace=0.15,hspace=0.1)
# # ax=[[] for _ in range(0,2)] 
# plt.rcParams.update({'font.size':10})
# ax=plt.subplot(projection=crg)
# # for idx,gs1,gs2,label,var in zip(range(0,2),[0,0],[0,1],['A','B'],['EAD_diff_nopro','EAD_diff_pro']):
#     # print(1)

# cmap = mpl.colormaps.get_cmap("RdBu_r")
# bins = np.arange(-40,45,5)
    
# norm = mpl.colors.BoundaryNorm(bins, cmap.N)
# # df.plot(ax=ax[idx],column='RP200_diff_pro', cmap=cmap1, norm=norm, legend='True', scheme='user_defined',classification_kwds={'bins':bins},
# #             missing_kwds={
# #         "color": "gainsboro",
# #         "edgecolor": "k",
# #         "label": "No Flood Damages",
# #     },
# #     linewidth=0.3, edgecolor='k')
# var = 'EAD_diff_pro'
# df.plot(ax=ax,column=var, cmap=cmap, norm=norm, linewidth=0.3, edgecolor='k',missing_kwds={
#         "color": "gainsboro",
#         "edgecolor": "k",
#         "label": "No Flood Damages",
#     })
# ax.set_extent([-180, 180, -90, 90], crs=crg)
# cax1= fig.add_axes([0.12,0.3,0.2,0.02]) #[left, bottom, width, height]

# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# ax.spines['bottom'].set_visible(False)
# ax.spines['left'].set_visible(False)
# ax.axis("off")

    

# cax1= fig.add_axes([0.32,0.03,0.35,0.02]) #[left, bottom, width, height]
# sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
# cbar= plt.colorbar(sm,ax=ax[idx],cax=cax1,cmap=cmap,spacing='uniform',
#                     orientation='horizontal',extend='both',boundaries=bins,
#                     ticks=[-40,-20,0,20,40]) #relative difference
# cbar.set_label('EAD difference (%)',size=12)
# cbar.ax.set_xticklabels(['-40', '-20', '0',  '20', '40']) 


# filename='EAD_diff_pro_new.png'
# outputdir = r'C:\Users\hli490\Desktop\Paper2\Figures\EAD'
# plt.savefig(os.path.join(outputdir,filename),format='png',dpi=300,bbox_inches='tight') 



### EAD and RP200 damage together
# for idx,gs1,gs2,label,var in zip(range(0,4),[0,0,1,1],[0,1,0,1],['A','B','C','D'],['EAD_diff_nopro','RP200_diff_nopro','EAD_diff_pro','RP200_diff_pro']):
#     # print(1)
    
#     ax[idx] = plt.subplot(gs[gs1,gs2], projection=crg)
#     ax[idx].text(-0.05,0.95,label,fontweight='bold',transform=ax[idx].transAxes) 
    
#     if idx==0 or idx==2: 
#         cmap = mpl.colormaps.get_cmap("RdBu_r")
#         bins = np.arange(-40,45,5)
        
#     else: 
#         # cmap = mpl.colormaps.get_cmap("RdPu")
#         # cmap = plt.get_cmap('PuOr')
#         cmap = plt.get_cmap('PiYG_r')
#         # cmap1 = mpl.colormaps.get_cmap("RdBu_r")
#         bins = np.arange(-80,210,10)
        

#     # norm = mpl.colors.BoundaryNorm(bins, cmap.N)
#     norm = MidpointNormalize(min(bins), max(bins), 0.)
#     # df.plot(ax=ax[idx],column='RP200_diff_pro', cmap=cmap1, norm=norm, legend='True', scheme='user_defined',classification_kwds={'bins':bins},
#     #             missing_kwds={
#     #         "color": "gainsboro",
#     #         "edgecolor": "k",
#     #         "label": "No Flood Damages",
#     #     },
#     #     linewidth=0.3, edgecolor='k')
#     df.plot(ax=ax[idx],column=var, cmap=cmap, norm=norm, linewidth=0.3, edgecolor='k',missing_kwds={
#             "color": "gainsboro",
#             "edgecolor": "k",
#             "label": "No Flood Damages",
#         })
#     ax[idx].set_extent([-180, 180, -60, 90], crs=crg)
#     # cax1= fig.add_axes([0.12,0.3,0.2,0.02]) #[left, bottom, width, height]
    
#     ax[idx].spines['top'].set_visible(False)
#     ax[idx].spines['right'].set_visible(False)
#     ax[idx].spines['bottom'].set_visible(False)
#     ax[idx].spines['left'].set_visible(False)
#     ax[idx].axis("off")
    
#     # add colorbar
#     if idx==2:
#         # cax1= fig.add_axes([0.46,0.2,0.2,0.02]) #[left, bottom, width, height]
#         cax1= fig.add_axes([0.2,0.02,0.2,0.02]) #[left, bottom, width, height]
#         sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
#         cbar= plt.colorbar(sm,ax=ax[idx],cax=cax1,cmap=cmap,spacing='uniform',
#                             orientation='horizontal',extend='both',boundaries=bins,
#                             ticks=[-40,-20,0,20,40]) #relative difference
#         cbar.set_label('EAD difference (%)',size=12)
#         cbar.ax.set_xticklabels(['-40', '-20', '0',  '20', '40']) 
        
#     if idx==3:
#         # cax1= fig.add_axes([0.46,0.2,0.2,0.02]) #[left, bottom, width, height]
#         cax1= fig.add_axes([0.62,0.02,0.2,0.02]) #[left, bottom, width, height]
#         sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
#         cbar= plt.colorbar(sm,ax=ax[idx],cax=cax1,cmap=cmap,spacing='uniform',
#                             orientation='horizontal',extend='max',boundaries=bins,
#                             ticks=[-80,-40, 0, 40, 80, 120, 160, 200]) #relative difference
#         cbar.set_label('200-year damage difference (%)',size=12)
#         cbar.ax.set_xticklabels([-80,-40, 0, 40, 80, 120, 160, 200]) 

# filename='diff_map.png'
# outputdir = r'C:\Users\hli490\Desktop\Paper2\Figures\EAD'
# plt.savefig(os.path.join(outputdir,filename),format='png',dpi=600,bbox_inches='tight') 


# Path and folder location
# path = r'C:\Users\hli490\Desktop\Paper2'

# # df_ead['FID'] = df_ead['iso']

# df_ead = pd.read_csv(r'C:\Users\hli490\Desktop\Paper2\Impact\EAD_diff_107.csv')
# df_ead['diff'] = -(df_ead['Annual_max']-df_ead['Complete_dependence'])/df_ead['Annual_max']*100

# #### Plotting
# # read the shapefile using geopandas
# shpfilename = r'C:\Users\hli490\Desktop\Paper2\Geogunit\geogunit_107.shp'
# shp = gpd.read_file(shpfilename,engine="pyogrio")
# shp = (shp.sort_values(by='ISO')).reset_index(drop=True)
# shp['FID'] = shp.index
# # shp = shp.loc[ind]
# df = shp.merge(df_ead,how='left',on='FID')

# # world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))

# #### EAD

# ### Difference
# # set a variable that will call whatever column we want to visualise on the map
# # variable = 'EAD_diff_re'
# variable = 'diff'
# # create figure and axes for Matplotlib

# crg = ccrs.PlateCarree()
# fig=plt.figure(figsize=(12,6))
# ax=plt.subplot(projection=crg)
# cmap = mpl.colormaps.get_cmap("RdBu_r")
# # bins = np.arange(0.0,1.1,0.1)
# # bins = np.arange(0.0,1.1,0.1)
# bins = np.arange(-20,21,2)
# bins = [-20,-15,-10,-5,0,5,10,15,20]
# norm = mpl.colors.BoundaryNorm(bins, cmap.N)
# df.plot(ax=ax,column=variable, cmap=cmap,scheme='user_defined',classification_kwds={'bins':bins},
#             missing_kwds={
#         "color": "gainsboro",
#         "edgecolor": "k",
#         "label": "No Flood Damages",
#     },
#     linewidth=0.3, edgecolor='k')

# df.plot(ax=ax,column=variable, cmap=cmap,scheme='user_defined',classification_kwds={'bins':bins},
#             missing_kwds={
#         "color": "gainsboro",
#         "edgecolor": "k",
#         "hatch": "///",
#         "label": "No Data",
#     },
#     linewidth=1, edgecolor='k')

# ax.set_extent([-180, 180, -60, 90], crs=crg)
# # cax1= fig.add_axes([0.12,0.3,0.2,0.02]) #[left, bottom, width, height]
# cax1= fig.add_axes([0.46,0.2,0.2,0.02]) #[left, bottom, width, height]
# sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
# cbar= plt.colorbar(sm,ax=ax,cax=cax1,cmap=cmap,spacing='uniform',
#                     orientation='horizontal',extend='both',boundaries=bins,
#                     ticks=[-20,-10,0,10,20]) #relative difference
# cbar.set_label('EAD difference (%)',size=10)
# cbar.ax.set_xticklabels(['-20', '-10', '0', '10', '20']) 
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# ax.spines['bottom'].set_visible(False)
# ax.spines['left'].set_visible(False)
# ax.axis("off")

# filename='diff_EAD.png'
# outputdir = r'C:\Users\hli490\Desktop\Paper2\Figures\EAD'
# plt.savefig(os.path.join(outputdir,filename),format='png',dpi=600,bbox_inches='tight') 

# ############## rp200 damage
# df_rp200 = pd.read_csv(r'C:\Users\hli490\Desktop\Paper2\Impact\rp200_diff.csv')
# df_rp200['rp200_diff'] = -(df_rp200['Annual_max']-df_rp200['Complete_dependence'])/df_rp200['Annual_max']*100
# df = shp.merge(df_rp200,how='left',on='FID')
# # set a variable that will call whatever column we want to visualise on the map
# # variable = 'EAD_diff_re'
# variable = 'rp200_diff'
# # create figure and axes for Matplotlib

# crg = ccrs.PlateCarree()
# fig=plt.figure(figsize=(12,6))
# ax=plt.subplot(projection=crg)
# cmap = mpl.colormaps.get_cmap("Greens")
# # bins = np.arange(0.0,1.1,0.1)
# # bins = np.arange(0.0,1.1,0.1)
# # bins = np.arange(-20,21,2)
# bins = [0,10,20,30,40,50,60,70,80,90,100]
# norm = mpl.colors.BoundaryNorm(bins, cmap.N)
# df.plot(ax=ax,column=variable, cmap=cmap,scheme='user_defined',classification_kwds={'bins':bins},
#             missing_kwds={
#         "color": "gainsboro",
#         "edgecolor": "k",
#         "label": "No Flood Damages",
#     },
#     linewidth=0.3, edgecolor='k')

# df.plot(ax=ax,column=variable, cmap=cmap,scheme='user_defined',classification_kwds={'bins':bins},
#             missing_kwds={
#         "color": "gainsboro",
#         "edgecolor": "k",
#         "hatch": "///",
#         "label": "No Data",
#     },
#     linewidth=1, edgecolor='k')

# ax.set_extent([-180, 180, -60, 90], crs=crg)
# # cax1= fig.add_axes([0.12,0.3,0.2,0.02]) #[left, bottom, width, height]
# cax1= fig.add_axes([0.46,0.2,0.2,0.02]) #[left, bottom, width, height]
# sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
# cbar= plt.colorbar(sm,ax=ax,cax=cax1,cmap=cmap,spacing='uniform',
#                     orientation='horizontal',extend='max',boundaries=bins,
#                     ticks=[0,20,40,60,80,100]) #relative difference
# cbar.set_label('Difference (%)',size=10)
# cbar.ax.set_xticklabels(['0', '20', '40', '60', '80', '100']) 
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# ax.spines['bottom'].set_visible(False)
# ax.spines['left'].set_visible(False)
# ax.axis("off")

# filename='diff_rp200.png'
# outputdir = r'C:\Users\hli490\Desktop\Paper2\Figures\EAD'
# plt.savefig(os.path.join(outputdir,filename),format='png',dpi=600,bbox_inches='tight') 
