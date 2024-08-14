# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 11:33:53 2023
Aggregate damage per cell to geogunit 107

@author: hli490
"""

import os
import os.path
import sys
import pandas as pd
import numpy as np
import xarray as xr

#### make the geogunit map/impact map upside down
#### crop the impact map of glofris to the same cells

# path = r'c:\Users\hli490\Desktop\Spatial_dependence\Paper2\\'

# basin = str(sys.argv[1])
# cluster = int(sys.argv[2])
# start = int(sys.argv[3])
path = '/projects/0/FWC2/Spatial_dependence/Paper2/'
# event = int(sys.argv[3])
# event_num_loc = os.path.join(path, 'Eventset/Event_num/num_events_'+basin+'.csv')
# event_num = pd.read_csv(event_num_loc)
# num = event_num.iloc[cluster-1,1]
#num = 11952
num = 119514
start=0

# geogunit to aggregate
pro_status = 'Nopro' # with protection standard or not
geogunit = 0
geo_loc = os.path.join(path,'Geogunit/geogunit_0_all.nc')
geo = xr.load_dataset(geo_loc)

# crop to the extent
extent = pd.read_csv(os.path.join(path, 'Hazard/Extent/NEA.csv'))
minx = min(extent.iloc[0, 1],extent.iloc[1, 1],extent.iloc[25, 1])
maxx = max(extent.iloc[0, 2],extent.iloc[1, 2],extent.iloc[25, 2]) 
miny = min(extent.iloc[0, 3],extent.iloc[1, 3],extent.iloc[25, 3])
maxy = max(extent.iloc[0, 4],extent.iloc[1, 4],extent.iloc[25, 4]) 

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

for event in range(start,num):
    # Impact folder
    dmg_fo = os.path.join(path,'Impact_nopro/splittest/Inun_event{:06d}_Built-up_Area.nc'.format(event))
    # if not os.path.isfile(dmg_fo): sys.exit()
    if not os.path.isfile(dmg_fo): continue
    dmg = xr.open_dataset(dmg_fo)
    np_dmg = np.array(dmg.Impact_Results.values)
    #pop_fo = os.path.join(path,'Impact_nopro/splittest/Inun_event{:06d}_POPexp.nc'.format(event))
    #pop = xr.open_dataset(pop_fo)
    #np_pop = np.array(pop.Impact_Results.values)
       
    # Mask layer extent (needs to be updated using the coordinates of stations later), WGS84 coordinates
    np_et = geo_clip.copy()
    
    # np_et[np.isnan(np_et)]=-9999
    ids = np.unique(np_et)  # subnational area ids
    #df = pd.DataFrame(0,index=np.arange(len(ids)),columns=['FID', 'Dmg','POPexp'])
    df = pd.DataFrame(0,index=np.arange(len(ids)),columns=['FID', 'Dmg'])
    df['FID']=ids
    
    for i in range(len(ids)):
        df.iloc[i,1]=np_dmg[np.where(np_et==ids[i])].sum()
       # df.iloc[i,2]=np_pop[np.where(np_et==ids[i])].sum()
    
    # impact results for all events
    if not 'df_dmg' in locals():
        df_dmg=pd.DataFrame(0,index=range(num),columns=df['FID'])
        df_dmg.index.name='eventid'
        
    #if not 'df_pop' in locals():
        #df_pop=pd.DataFrame(0,index=range(num),columns=df['FID'])
        #df_pop.index.name='eventid'
        
    df_dmg.iloc[event,:] = df.iloc[:,1]        
    #df_pop.iloc[event,:] = df.iloc[:,2]    
    
# save
out_fo = os.path.join(path,'Risk/'+pro_status+'/geogunit_{:d}/splittest'.format(geogunit))
if not os.path.isdir(out_fo): os.makedirs(out_fo)
df_dmg.to_csv(out_fo+'/damage_splittest.csv')
#df_pop.to_csv(out_fo+'/POPexp_splittest01.csv')




