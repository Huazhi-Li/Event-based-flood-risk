# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 13:22:03 2024

@author: hli490
"""

import os
import sys
import numpy as np
import pandas as pd
import xarray as xr

# path =  '/projects/0/FWC2/Spatial_dependence/Paper2'
path = 'c:/Users/hli490/Desktop/Paper2'
geounit = 107
rps = [2,5,10,25,50,100,250,500,1000]
j = 0
rp = rps[j]

# Inuput file locations
src_station_loc = os.path.join(path,'Hazard/Source_station/src_station_rp1000_globe.nc')
dmg_file_loc = os.path.join(path, 'Impact/Glofris')
geo_loc = os.path.join(path,'Geogunit/geogunit_107_all.nc')

# Read file
src_rp1000_file = xr.load_dataset(src_station_loc)
src_rp1000 = src_rp1000_file.src.values
geo_globe = xr.load_dataset(geo_loc)
geo = geo_globe.Geogunits.values
ds_geo_station = pd.read_pickle(os.path.join(path,'Station_geogunit107/geogunit_station_globe.pkl'))
ds_station_geo = pd.read_pickle(os.path.join(path,'Station_geogunit107/station_geoguint_globe.pkl'))

# read damages
dmg_globe = xr.load_dataset(os.path.join(dmg_file_loc,'inuncoast_historical_nosub_hist_rp{:04d}_0_Built-up_Area.nc'.format(rp)))
dmg = dmg_globe.Impact_Results.values # clip to the basin

# Damage percentage file
sts = ds_station_geo.iloc[:,0]
df_dmg = pd.DataFrame(index=range(len(sts)),columns = ['Station','FID','Dmg_perc'])
df_dmg['Station'] = sts
src = src_rp1000.copy()
# Damage percentage file
perc = []
FIDs = []

for i in range(len(sts)):
    
    st = sts[i]
    r_src, c_src = np.where(src==st)
    
    # get all units affected by this station
    geo_st = geo[r_src,c_src]
    units = np.unique(geo_st)
    FIDs.append(list(units))
    
    # calculate the total damage of this station
    dmg_total = dmg[r_src, c_src].sum()
    
    if dmg_total == 0:
        perc.append([0.0])
        continue
    
    # set up damage percentages
    dmg_perc = []
    for unit in units:
        # Get the points where both files show positive values
        idx = np.where(geo_st==unit)[0]
        r_unit, c_unit = r_src[idx], c_src[idx]
        dmg_unit = dmg[ r_unit, c_unit].sum()
        dmg_perc.append(dmg_unit/dmg_total)
        
    perc.append(dmg_perc)
    print(i)

df_dmg['FID'] = FIDs
df_dmg['Dmg_perc'] = perc
out_dir = os.path.join(path,'Station_geogunit107/dmg_perc_rp{:04d}_globe.pkl'.format(rp))
df_dmg.to_pickle(out_dir)

df_dmg_perc = pd.read_pickle(r'C:\Users\hli490\Desktop\Paper2\Station_geogunit107\dmg_perc_rp1000_globe_old.pkl')
n=0
for i in range(len(sts)):
    a1 = np.array(df_dmg['Dmg_perc'][i]).astype(float)
    # a2 = np.array(df_dmg_perc['Dmg_perc'][i]).astype(float)
    a3 = np.array(df_dmg['FID'][i]).astype(float)
    # a4 = np.array(df_dmg_perc['FID'][i]).astype(float)
    # if not np.array_equal(a1, a2): 
    #     print(i)
    #     n=n+1
    if np.array_equal(a1,[0.0]): 
        print(i)
        n=n+1
    

dmg_perc = pd.read_pickle(r'C:\Users\hli490\Desktop\Paper2\Impact\Input\dmg_perc_globe.pkl')
dmg_perc['FID'] = df_dmg['FID']
    
out_dir = os.path.join(path,'Impact/Input/dmg_perc_globe.pkl')
dmg_perc.to_pickle(out_dir)
    
    