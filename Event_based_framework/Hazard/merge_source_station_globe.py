# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 09:59:32 2024

@author: hli490
"""

import os
import numpy as np
import pandas as pd
import xarray as xr
import geopy.distance as gd
# import matplotlib.pyplot as plt

path = 'c:/Users/hli490/Desktop/Paper2'

basins = ['IO','NEA','NA','NWA','NWP','Oceania','SA','SEA','SWA','SWP']
df_station = pd.read_csv(os.path.join(path,'Geogunit\Station\station_info.csv'))
inun_file_loc = os.path.join(path, 'Hazard/deltares_2021/HIST/inuncoast_historical_nosub_hist_rp1000_0.nc')
inun_globe = xr.open_dataset(inun_file_loc)
src_value = inun_globe.inun.values[0,:,:]
src_value[:,:] = np.nan

lats0 = inun_globe['lat'].values
lons0 = inun_globe['lon'].values

n=0
overlap = 0

for basin in basins:
    # read the basin source station file
    src_loc = os.path.join(path,'Hazard/Source_station/'+basin+'/src_station_rp1000_{:s}.nc'.format(basin))
    # src_loc = os.path.join(path,'Hazard/Source_station/'+basin+'/src_station_rp0005_{:s}.nc'.format(basin)) # for testing
    src_file = xr.open_dataset(src_loc)
    src_value_basin = src_file.src.values
    src_value_basin = np.where(src_value_basin<0, np.nan, src_value_basin)
    
    # get the extend of the basin file
    minx, maxx, miny, maxy = src_file.lon.values.min(), src_file.lon.values.max(), src_file.lat.values.min(), src_file.lat.values.max()
       
    # clip the global src file into the basin area
    # latitude lower and upper index
    latli = np.argmin(np.abs(lats0 - miny))
    latui = np.argmin(np.abs(lats0 - maxy))

    # longitude lower and upper index
    lonli = np.argmin(np.abs(lons0 - minx))
    lonui = np.argmin(np.abs(lons0 - maxx))
    
    src_globe_basin = src_value[latli:latui+1, lonli:lonui+1]
    
    # update the src value in the global src file
    src_globe_basin[np.isnan(src_globe_basin)] = src_value_basin[np.isnan(src_globe_basin)] # replace the nans with the basin src values
    
    # get the cells of positive src values
    r_src, c_src = np.where(src_value_basin>=0)
    idx = np.where(src_value_basin[r_src,c_src]!=src_globe_basin[r_src,c_src])[0]
    
    # check if there are overlaps between two src files
    if len(idx) > 0: 
    
        for i in idx:
            # Get the coordinates of the cell, and two stations from two src files
            point_cell = src_file.lat.values[r_src[i]], src_file.lon.values[c_src[i]]
            station1 = src_value_basin[r_src[i],c_src[i]]
            station2 = src_globe_basin[r_src[i],c_src[i]]
            point_s1 = df_station.loc[df_station['station']==station1,'lat'].iloc[0],df_station.loc[df_station['station']==station1,'lon'].iloc[0]
            point_s2 = df_station.loc[df_station['station']==station2,'lat'].iloc[0],df_station.loc[df_station['station']==station2,'lon'].iloc[0]
            # calculate the distances
            d1 = gd.geodesic(point_cell, point_s1)
            d2 = gd.geodesic(point_cell, point_s2)
            # select the closer station
            if d1<=d2: 
                src_globe_basin[r_src[i],c_src[i]] = src_value_basin[r_src[i],c_src[i]]
            else:
                src_globe_basin[r_src[i],c_src[i]] = src_globe_basin[r_src[i],c_src[i]]
    
    # update the src value in the global src file       
    src_value[latli:latui+1, lonli:lonui+1] = src_globe_basin

src_globe = inun_globe.copy()   
src_globe = src_globe.drop_dims('time')
src_globe['src']=(['lat','lon'], src_value)
src_globe.attrs['Title'] = 'Source station globe'
src_globe.attrs['Author'] = 'Huazhi Li'
src_globe.attrs['Institution'] = 'VU Amsterdam'
src_globe.attrs['Source'] = 'GLOFRIS global coastal inundation layer rp1000'
del src_globe.attrs['config_file']
del src_globe.attrs['project']
del src_globe.attrs['references']
del src_globe.attrs['history']
del src_globe.attrs['title']
del src_globe.attrs['institution']
del src_globe.attrs['source']

# Save to nc 
outdir = os.path.join(path, 'Hazard/Source_station/src_station_rp1000_globe.nc')
src_globe.to_netcdf(path=outdir,mode='w')
src_globe.close()



