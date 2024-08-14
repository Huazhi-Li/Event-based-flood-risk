# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 11:33:53 2023
Aggregate damage per cell to geogunit 107

@author: hli490
"""

import sys
import os
import os.path
import pandas as pd
import numpy as np
import xarray as xr

#### make the geogunit map/impact map upside down
#### crop the impact map of glofris to the same cells

# path = r'c:\Users\hli490\Desktop\Spatial_dependence\Paper2\\'
# basin = 'NEA'
# cluster = 1

basin = str(sys.argv[1])
cluster = int(sys.argv[2])
path = '/projects/0/FWC2/Spatial_dependence/Paper2/'

# Read the events
event_folder = os.path.join(path,'Eventset/'+basin+'/rp_event_set_cluster{:d}.pkl'.format(cluster))
event_rp =  pd.read_pickle(event_folder)
num = len(event_rp)

# geogunit to aggregate
pro_status = 'Nopro' # with protection standard or not
geogunit = 0
geo_loc = os.path.join(path,'Geogunit/geogunit_{:d}_all.nc'.format(geogunit))
geo = xr.load_dataset(geo_loc)

# clip to extent
# extent boundary
extent = pd.read_csv(os.path.join(path, 'Hazard/Extent/'+basin+'.csv'))
cluster_link = np.array(pd.read_csv(os.path.join(path, 'Hazard/Cluster_link/'+basin+'.csv')))
cluster_all = cluster_link[cluster-1,~np.isnan(cluster_link[cluster-1,:])] # all linked clusters

# get the extent coordinates
for i in range(len(cluster_all)):
    cluster0 = int(cluster_all[i])-1
    minx0, maxx0, miny0, maxy0 = extent.iloc[cluster0,1], extent.iloc[cluster0,2], extent.iloc[cluster0,3], extent.iloc[cluster0,4]
    
    if i==0: 
        minx, maxx, miny, maxy = minx0, maxx0, miny0, maxy0
    else:
        minx, maxx, miny, maxy = min(minx0,minx), max(maxx0,maxx), min(miny0,miny), max(maxy0,maxy)

lats0 = geo.lat.values
lons0 = geo.lon.values

# latitude lower and upper index
latli = np.argmin(np.abs(lats0 - miny))
latui = np.argmin(np.abs(lats0 - maxy))

# longitude lower and upper index
lonli = np.argmin(np.abs(lons0 - minx))
lonui = np.argmin(np.abs(lons0 - maxx))

# variable 
geo_clip = geo.Geogunits.values.copy()

if (len(geo_clip.shape) == 3):
    geo_clip = geo_clip[0, latli:latui+1, lonli:lonui+1]
else:
    geo_clip = geo_clip[latli:latui+1, lonli:lonui+1]

np_et = geo_clip.copy()
ids = np.unique(np_et)  # subnational area ids

# aggragate impact
for event in range(num):
    # Impact folder
    dmg_fo = os.path.join(path,'Impact_nopro/'+basin+'/cluster{:d}'.format(cluster)+'/Inun_event{:06d}_Built-up_Area.nc'.format(event))
    # if not os.path.isfile(dmg_fo): sys.exit()
    if not os.path.isfile(dmg_fo): continue
    dmg = xr.open_dataset(dmg_fo)
    np_dmg = np.array(dmg.Impact_Results.values)
       
    # Mask layer extent (needs to be updated using the coordinates of stations later), WGS84 coordinates
    df = pd.DataFrame(0,index=np.arange(len(ids)),columns=['FID', 'Dmg'])
    df['FID']=ids
    
    for i in range(len(ids)):
        df.iloc[i,1]=np_dmg[np.where(np_et==ids[i])].sum()
       # df.iloc[i,2]=np_pop[np.where(np_et==ids[i])].sum()
    
    # impact results for all events
    if not 'df_dmg' in locals():
        df_dmg=pd.DataFrame(0,index=range(num),columns=df['FID'])
        df_dmg.index.name='eventid'
        
    df_dmg.iloc[event,:] = df.iloc[:,1]        
    
# save
out_fo = os.path.join(path,'Risk/'+pro_status+'/geogunit_{:d}/'.format(geogunit)+basin)
if not os.path.isdir(out_fo): os.makedirs(out_fo)
df_dmg.to_csv(out_fo+'/damage_cluster{:d}.csv'.format(cluster))





