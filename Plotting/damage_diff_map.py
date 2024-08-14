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
import mapclassify

# Path and folder location
basin = 'NEA'
path = r'C:\Users\hli490\Desktop\Paper2'
iso_cluster_loc = os.path.join(path,'Geogunit\Station\iso_cluster_{:s}.pkl'.format(basin))
cluster_link_loc = os.path.join(path,'Hazard/Cluster_link/{:s}.csv'.format(basin))

    
## Event set damage
df_event = pd.read_csv(r'C:\Users\hli490\Desktop\Paper2\Risk\Nopro\geogunit_0\country\dmg_country_NEA.csv')
df_event.index=df_event.iso

## Glofris damage
df_glo = pd.read_csv(r'C:\Users\hli490\Desktop\Paper2\Risk\Nopro\glofris\glofris_Built-up_Area_geogunit0.csv')
df_glo.index=df_glo.iso

# make diff map
rp_diff = 250
ind = np.where(np.in1d(df_glo.index,df_event.index))
df_diff = pd.concat([df_event[str(rp_diff)], df_glo.loc[ind,str(rp_diff)]], axis=1)
df_diff.columns=['Event','Glofris']
df_diff['diff']=(df_diff['Glofris']-df_diff['Event'])/df_diff['Event']*100
df_diff['FID'] = df_event.index
df_diff.index=range(len(df_diff))
df_diff.loc[df_diff['FID']==194,'diff']=np.nan

#### Plotting
# # get natural earth data (http://www.naturalearthdata.com/)

# # get country borders
# # resolution = '10m'
# # category = 'cultural'
# # name = 'admin_0_countries'

# # shpfilename = shapereader.natural_earth(resolution, category, name)

# read the shapefile using geopandas
shpfilename = r'C:\Users\hli490\Desktop\Paper2\Geogunit\geogunit_0.shp'
shp = gpd.read_file(shpfilename,engine="pyogrio")
shp = (shp.sort_values(by='ISO')).reset_index(drop=True)
shp['FID'] = shp.index
# shp = shp.loc[ind]
df = shp.merge(df_diff,how='left',on='FID')
df['Event']=df['Event']/1e6


# world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))

### Difference
# set a variable that will call whatever column we want to visualise on the map
variable = 'diff'

# create figure and axes for Matplotlib

crg = ccrs.PlateCarree()
fig=plt.figure(figsize=(4*3.13,4*3.13))
ax=plt.subplot(projection=crg)
# cmap = mpl.colormaps.get_cmap('OrRd')
cmap = mpl.colormaps.get_cmap('Greens')
bins = [0,2,5,10,20,30,40,50,60,70,80,90,100]
# bins = [-5,0,5,10,25,50,75,80,90,100]
norm = mpl.colors.BoundaryNorm(bins, cmap.N)
df.plot(ax=ax,column=variable, cmap=cmap,scheme='user_defined',classification_kwds={'bins':bins},
            missing_kwds={
        "color": "gainsboro",
        "edgecolor": "black",
        "label": "No Data",
    },
    linewidth=1, edgecolor='k')

# df.plot(ax=ax,column=variable, cmap=cmap,scheme='user_defined',classification_kwds={'bins':bins},
#             missing_kwds={
#         "color": "gainsboro",
#         "edgecolor": "k",
#         "hatch": "/",
#         "label": "No Data",
#     },
#     linewidth=1, edgecolor='k')

# ax=df.plot(column=variable, cmap=cmap,scheme='quantiles', 
#         missing_kwds={
#         "color": "grey",
#         "edgecolor": "black",
#         "hatch": "///",
#         "label": "No Data",
#     },
#     edgecolor='k',legend=True, legend_kwds=dict(loc='upper left'),ax=ax)

ax.set_extent([-25, 45, 15, 75], crs=crg)
# ax.set_title('Flood damage RP50',fontsize=18,fontweight='bold')
cax1= fig.add_axes([0.152,0.12,0.7,0.02]) #[left, bottom, width, height]
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
cbar= plt.colorbar(sm,ax=ax,cax=cax1,cmap=cmap,spacing='uniform',
                    orientation='horizontal',extend='max',boundaries=bins,
                    ticks=[0,2,5,10,20,30,40,50,60,70,80,90,100]) #relative difference
cbar.set_label('Damage difference (%)',size=14)
cbar.ax.tick_params(labelsize=12)
filename='diff_rp100.png'
outputdir = r'C:\Users\hli490\Desktop\Paper2\Figures\Diff_maps'
plt.savefig(os.path.join(outputdir,filename),format='png',dpi=600,bbox_inches='tight') 


#### Damage
# set a variable that will call whatever column we want to visualise on the map
variable = 'Event'

# create figure and axes for Matplotlib

crg = ccrs.PlateCarree()
fig=plt.figure(figsize=(4*3.13,4*3.13))
ax=plt.subplot(projection=crg)
# cmap = mpl.colormaps.get_cmap('OrRd')
cmap = mpl.colormaps.get_cmap('Reds')
cmap.set_under(color='gainsboro')  
bins = [0,10,100,1000,10000,100000,250000,500000]
norm = mpl.colors.BoundaryNorm(bins, cmap.N)
df.plot(ax=ax,column=variable, cmap=cmap,scheme='user_defined',classification_kwds={'bins':bins},
            missing_kwds={
        "color": "gainsboro",
        "edgecolor": "black",
        "label": "No Data",
    },vmin=0.0000001,
    linewidth=1, edgecolor='k')

# df.plot(ax=ax,column=variable, cmap=cmap,scheme='user_defined',classification_kwds={'bins':bins},
#             missing_kwds={
#         "color": "gray",
#         "edgecolor": "gainsboro",
#         "hatch": "/",
#         "label": "No Data",
#     },
#     linewidth=1, edgecolor='gainsboro')

# ax=df.plot(column=variable, cmap=cmap,scheme='quantiles', 
#         missing_kwds={
#         "color": "grey",
#         "edgecolor": "black",
#         "hatch": "///",
#         "label": "No Data",
#     },
#     edgecolor='k',legend=True, legend_kwds=dict(loc='upper left'),ax=ax)

ax.set_extent([-25, 45, 15, 75], crs=crg)
# ax.set_title('Flood damage RP50',fontsize=18,fontweight='bold')
cax1= fig.add_axes([0.152,0.12,0.7,0.02]) #[left, bottom, width, height]
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
cbar= plt.colorbar(sm,ax=ax,cax=cax1,cmap=cmap,spacing='uniform',
                    orientation='horizontal',extend='max',boundaries=bins,
                    ticks=[0,10,100,1000,10000,100000,250000,500000]) #relative difference
cbar.set_label('Economic damage (mil. USD)',size=14)
cbar.ax.tick_params(labelsize=12)
filename='dmg_rp100.png'
outputdir = r'C:\Users\hli490\Desktop\Paper2\Figures\Diff_maps'
plt.savefig(os.path.join(outputdir,filename),format='png',dpi=600,bbox_inches='tight')

    
fig=plt.figure(figsize=(4*3.13,4*3.13))
ax=plt.subplot(projection=crg)
# cmap = mpl.colormaps.get_cmap('OrRd')
ax.hist()

ax.set_extent([-25, 45, 15, 75], crs=crg)
# ax.set_title('Flood damage RP50',fontsize=18,fontweight='bold')
cax1= fig.add_axes([0.152,0.12,0.7,0.02]) #[left, bottom, width, height]
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
cbar= plt.colorbar(sm,ax=ax,cax=cax1,cmap=cmap,spacing='uniform',
                    orientation='horizontal',extend='max',boundaries=bins,
                    ticks=[0,10,100,1000,10000,100000,250000,500000]) #relative difference
cbar.set_label('Economic damage (mil. USD)',size=14)
cbar.ax.tick_params(labelsize=12)
filename='dmg_rp100.png'
outputdir = r'C:\Users\hli490\Desktop\Paper2\Figures\Diff_maps'
plt.savefig(os.path.join(outputdir,filename),format='png',dpi=600,bbox_inches='tight')




####
# world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
# ax2 = world.plot(figsize=(15,15), edgecolor=u'white', color='gray')


# countries = ['Germany', 'Norway', 'Russia', 'China', 'Japan','Netherlands']
# for name in (countries):
#     world.loc[world['name'].eq(name)].plot(edgecolor=u'white', color=np.random.rand(3,), ax=ax2)

# ax2.axis('scaled')
# plt.show()