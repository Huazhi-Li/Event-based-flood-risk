# -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:28:38 2024

@author: hli490
"""

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cmplt
import cartopy
import cartopy.crs as ccrs
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# Setting working directory
path = r'C:\Users\hli490\OneDrive - Vrije Universiteit Amsterdam\Desktop\Paper2'

continent = 'South America'
dmg = pd.read_csv(os.path.join(path,'Impact/Worst_year/{:s}.csv'.format(continent)))
dmg[dmg==0] = np.nan
geo = pd.read_csv(os.path.join(path,'Geogunit/geogunit_107_0_list.csv'))
dmg = dmg.merge(geo,how='left',on='FID')

[dmg['dmg_nopro'].quantile(q=0.95)/1e6,
dmg['dmg_nopro'].quantile(q=0.75)/1e6,
dmg['dmg_nopro'].median()/1e6,
dmg['dmg_nopro'].quantile(q=0.25)/1e6,
dmg['dmg_nopro'].quantile(q=0.05)/1e6]

# read the shapefile using geopandas
shpfilename = os.path.join(path,'Geogunit/geogunit_107.shp')
shp = gpd.read_file(shpfilename,engine="pyogrio")


country_filename = os.path.join(path,'Geogunit/geogunit_0.shp')
country = gpd.read_file(country_filename,engine="pyogrio")
# country = gpd.clip(country,[0, 50, 25, 100])
df = shp.merge(dmg,how='left',on='FID')
# df = gpd.clip(df,[0, 50, 25, 100])
    
# x = np.array(dmg['dmg_nopro'])
# np.quantile(x[x>0],0.75)

## plotting
# crg = ccrs.PlateCarree()
# crg = ccrs.Mercator()
# fig=plt.figure(figsize=(10,10)) # width, height
# plt.rcParams.update({'font.size':10})
# ax=plt.subplot(projection=crg)
# # ax[idx].text(0.01,0.95,label,fontweight='bold',transform=ax[idx].transAxes,bbox=dict(facecolor='none', edgecolor='k')) 
# # ax.text(0.01,0.95,label,fontweight='bold',transform=ax[idx].transAxes)
# cmap = mpl.colormaps.get_cmap('RdYlBu_r')
# bins = np.array([0,1,10,100,1000,10000,100000,200000])*1e6 
# norm = mpl.colors.BoundaryNorm(bins, cmap.N)

# df.plot(ax=ax,column='dmg', cmap=cmap, norm=norm, linewidth=0.2,  missing_kwds={
#         "color": "darkgray",
#         "edgecolor": "none"
#     })

# country.plot(ax=ax, facecolor='none', linewidth=0.2, edgecolor='k')
# ax.add_feature(cartopy.feature.LAND, edgecolor='black')

# ax[idx].set_extent([-9.5, 16.5, 42, 62], crs=crg)
# ax.set_extent([-12, 45, 30, 72], crs=ccrs.PlateCarree()) # europe
# ax.set_extent([25, 160, -20, 50], crs=crg) # asia
# ax.set_extent([110,180,-50,0 ], crs=crg) # ocenia
# ax.set_extent([-30,70,-40,40 ], crs=crg) # africa
# ax.set_extent([-130,-50,0,70], crs=crg) # north america
# ax.set_extent([-100,-30,-60,20], crs=crg) # south america

# gl = ax.gridlines(crs=ccrs.PlateCarree(),linewidth=0.5, color='none', alpha=0.5, 
#                         draw_labels=True,x_inline=False, y_inline=False)
# gl.xlocator = mticker.FixedLocator([-10, 0 ,10,20,30,40])
# gl.ylocator = mticker.FixedLocator([35,45,55,65])
# gl.bottom_labels = True
# gl.left_labels   = True
# gl.top_labels    = False
# gl.right_labels  = False
# gl.xformatter = LONGITUDE_FORMATTER
# gl.yformatter = LATITUDE_FORMATTER

# cax1= fig.add_axes([0.2,0.08,0.65,0.02]) #[left, bottom, width, height]
# sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
# cbar= plt.colorbar(sm,ax=ax,cax=cax1,cmap=cmap,spacing='uniform',
#                     orientation='horizontal',extend='max',boundaries=bins,
#                     ticks=np.array([0,1,10,100,1000,10000,100000,200000])*1e6)
# cbar.set_label('Flood damage (mil. USD2005)',size=12,labelpad=5)
# cbar.ax.set_xticklabels(['0', '1', '10',  '100', '1000','10000','100000','>200000'])

# filename='Worst_year_{:s}.png'.format(continent)
# outputdir = r'C:\Users\hli490\Desktop\Paper2\Figures\Worst_year'
# plt.savefig(os.path.join(outputdir,filename),format='png',dpi=300,bbox_inches='tight')


crg = ccrs.PlateCarree()
fig=plt.figure(figsize=(10,6.3)) # width, height
gs=gridspec.GridSpec(1,2,wspace=0.1,hspace=0.1)
ax=[[] for _ in range(0,2)] 
plt.rcParams.update({'font.size':10})
mpl.rcParams['text.color']='k'

for idx,gs1,gs2,label,var in zip(range(0,2),[0,0],[0,1],['a)','b)'],['dmg_nopro','dmg_pro']):
    
    ax[idx]=plt.subplot(gs[gs1,gs2],projection=crg)
    # ax[idx].text(0.01,0.95,label,fontweight='bold',transform=ax[idx].transAxes,bbox=dict(facecolor='none', edgecolor='k')) 
    if continent=='Asia':
        ax[idx].text(0.95,0.01,label,transform=ax[idx].transAxes)
    elif continent=='Oceania' or continent=='South America':
        ax[idx].text(0.01,0.01,label,transform=ax[idx].transAxes)
    else:
        ax[idx].text(0.01,0.95,label,transform=ax[idx].transAxes)
    cmap = mpl.colormaps.get_cmap('RdYlBu_r')
    # bins = np.array([0,1,10,100,1000,10000,100000,200000])*1e6 
    bins = np.array([0,1,5,10,25,50,100,250,500,1000])*1e6 
    norm = mpl.colors.BoundaryNorm(bins, cmap.N)

    df.plot(ax=ax[idx],column=var, cmap=cmap, norm=norm, linewidth=0.1, missing_kwds={
            "color": "darkgray",
            "edgecolor": "none"
        })
    
    country.plot(ax=ax[idx], facecolor='none', linewidth=0.2, edgecolor='k')
    
    
    # ax[idx].set_extent([-9.5, 16.5, 42, 62], crs=crg)
    # ax[idx].set_extent([-12, 45, 30, 72], crs=crg) #eu
    # ax[idx].set_extent([-25, 65, -40, 40], crs=crg) #africa
    ax[idx].set_extent([-95, -30, -60, 15], crs=crg) #africa
    gl = ax[idx].gridlines(crs=ccrs.PlateCarree(),linewidth=0.5, color='none', alpha=0.5, 
                            draw_labels=True,x_inline=False, y_inline=False)
    gl.xlocator = mticker.FixedLocator([-90,-75,-60,-45,-30])
    gl.ylocator = mticker.FixedLocator([-60,-45,-30,-15,0,10])
    gl.bottom_labels = True
    gl.left_labels   = True
    gl.top_labels    = False
    gl.right_labels  = False
    
    if idx==1:
        gl.left_labels   = False
    

    
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER


# colorbar
cax1= fig.add_axes([0.315,0.08,0.4,0.02]) #[left, bottom, width, height]
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
cbar= plt.colorbar(sm,ax=ax[0],cax=cax1,cmap=cmap,spacing='uniform',
                    orientation='horizontal',extend='max',boundaries=bins,
                    ticks=bins)

cbar.set_label('Flood damage (mil. USD2005)',size=11,labelpad=5)
# cbar.ax.set_xticklabels(['0','1','10','50','100','500','1000','5000','10000'])
cbar.ax.set_xticklabels(['0','1','5','10','25','50','100','250','500','1000'])

filename='Worst_year_eu_continent_{:s}.png'.format(continent)
outputdir = os.path.join(path,'Figures/Worst_year')
plt.savefig(os.path.join(outputdir,filename),format='png',dpi=300,bbox_inches='tight')


# fig = plt.figure()
# crg = ccrs.Mercator(central_longitude=10)
# # crg = ccrs.Robinson()
# ax = plt.axes(projection=crg)
# ax.add_feature(cartopy.feature.LAND, edgecolor='black')
# ax.add_feature(cartopy.feature.OCEAN)
# ax.set_extent([-25, 45, 30, 72], crs=ccrs.PlateCarree())
# gl = ax.gridlines(crs=ccrs.PlateCarree(),linewidth=0.3, color='k', alpha=0.5,zorder=3)
# gl.xlocator = mticker.FixedLocator([-15, 0,15,30,45])
# gl.ylocator = mticker.FixedLocator([30,40,50,60,70])


# stacked bar chart
# out = dmg.groupby('FID_geog0').sum()
# FID_in = np.array(out.index).astype(int)
# country=[]
# for i in range(len(FID_in)):
#     idx=np.where(dmg['FID_geog0']==FID_in[i])[0][0]
#     country.append(dmg['ISO'][idx])
    
# df_bar = pd.DataFrame(index=range(2), columns=['Scenario']+country)
# df_bar['Scenario'] = ['No protection','Protection']
# df_bar.iloc[0,1:] = np.array(out['dmg_nopro'])/out['dmg_nopro'].sum()*100
# df_bar.iloc[1,1:] = np.array(out['dmg_pro'])/out['dmg_pro'].sum()*100

# ax[2]=plt.subplot(gs[0,2])
# df_bar.plot(ax=ax[2],x='Scenario',kind='bar',stacked=True)



