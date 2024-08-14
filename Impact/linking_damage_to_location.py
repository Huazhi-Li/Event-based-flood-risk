# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 14:09:21 2024

@author: hli490
"""


import os
import numpy as np
import pandas as pd
import xarray as xr
import sys

def nc_clip(minx, maxx, miny, maxy, data, variable):
    
    # clipping Mask layer extent, WGS84 coordinates

    lats0 = data['lat'].values
    lons0 = data['lon'].values

    # latitude lower and upper index
    latli = np.argmin(np.abs(lats0 - miny))
    latui = np.argmin(np.abs(lats0 - maxy))

    # longitude lower and upper index
    lonli = np.argmin(np.abs(lons0 - minx))
    lonui = np.argmin(np.abs(lons0 - maxx))

    # get the inundation
    data_clip = data[variable].values

    if (len(data_clip.shape) == 3):
        data_clip = data_clip[0, latli:latui+1, lonli:lonui+1]
    else:
        data_clip = data_clip[latli:latui+1, lonli:lonui+1]

    # ocean cells, pixel value nodata   
    nodata =-9999
    data_clip=np.nan_to_num(data_clip, nan=nodata)

    # non-inundated land cells, pixel value -1
    if variable=='inun':
        data_clip[data_clip==0] = -1
    
    return data_clip

path = 'c:/Users/hli490/Desktop/Paper2'
# path =  '/projects/0/FWC2/Spatial_dependence/Paper2'
# basin = 'NEA'
# Set up
rps = [2,5,10,25,50,100,250,500,1000]
# basins = ['IO','NA','NEA','NWA','NWP','Oceania','SA','SEA','SWA','SWP']
# rp_ind = int(sys.argv[1])
# Inuput file locations
dmg_file_loc = os.path.join(path, 'Impact_nopro/GLOFRIS')
src_station_loc = os.path.join(path,'Hazard/Source_station')

# Source station file for inundation layers of other rps
# rp_ind=8
# rp = rps[rp_ind]

rp=1000
# Read data
src_rp_file = xr.load_dataset(os.path.join(src_station_loc,'src_station_rp{:04d}_globe.nc'.format(rp)))
src_rp = src_rp_file.src.values
# minx, maxx, miny, maxy = src_rp_file.lon.values.min(), src_rp_file.lon.values.max(), src_rp_file.lat.values.min(), src_rp_file.lat.values.max()

# Global damage file (built-up area)
dmg_globe = xr.load_dataset(os.path.join(dmg_file_loc,'inuncoast_historical_nosub_hist_rp{:04d}_0_Built-up_Area.nc').format(rp))
# dmg = nc_clip(minx, maxx, miny, maxy, data=dmg_globe, variable='Impact_Results') # clip to the basin
dmg = dmg_globe.Impact_Results.values

# Aggragate damages to station-wise cells
src = src_rp.copy()
# src = np.where(src<0, np.nan, src)
out = pd.DataFrame(np.hstack((src.reshape(-1, 1), dmg.reshape(-1, 1))),columns=['Station','Damage'])
dmg_rp = out.groupby('Station')['Damage'].sum()

# Saving to basin damage files
# if rp==2: 
#     df = pd.DataFrame(index=dmg_rp.index,columns = ['Station','2','5','10','25','50','100','250','500','1000'])
#     df['Station'] = dmg_rp.index.astype(int)

# df[str(rp)] = dmg_rp

out_dir = r'c:\Users\hli490\Desktop\Paper2\Station_geogunit107\\'
if not os.path.isdir(out_dir): os.mkdir(out_dir)
dmg_rp.to_csv(os.path.join(out_dir,'damage_per_station_rp{:04d}_globe.csv'.format(rp)))


    






