# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 15:13:02 2023

@author: hli490
"""

# -*- coding: utf-8 -*-
"""
Creating the cluster inundation layer from the global inundation layer

@author: Huazhi Li (huazhi.li@vu.nl)
"""

import sys
import os.path
import numpy as np
import pandas as pd
import xarray as xr
    
# path = '/projects/0/FWC2/Spatial_dependence/Paper2/'
path = r'c:\Users\hli490\Desktop\Paper2\\'

# cluster info. this will be changed for parallel simulation on Snellius
# basin = str(sys.argv[1])
# cluster = int(sys.argv[2])
basin = 'NEA'
cluster = 27
# rp = int(sys.argv[3]) 
# basin = 'NEA'
# cluster = 1
rps = [2,5,10,25,50,100,250,500,1000]

# get the boundary coordinates of the potentially affected clusters
extent = pd.read_csv(os.path.join(path, 'Hazard/Extent/'+basin+'.csv'))
cluster_link = np.array(pd.read_csv(os.path.join(path, 'Hazard/Cluster_link/'+basin+'.csv')))
cluster_all = cluster_link[cluster-1,~np.isnan(cluster_link[cluster-1,:])] # all linked clusters

# get the extent coordinates
for i in range(len(cluster_all)):
    cluster0 = int(cluster_all[i])-1
    if basin=='NEA' and cluster0==27: continue
    
    minx0, maxx0, miny0, maxy0 = extent.iloc[cluster0,1], extent.iloc[cluster0,2], extent.iloc[cluster0,3], extent.iloc[cluster0,4]
    
    if i==0: 
        minx, maxx, miny, maxy = minx0, maxx0, miny0, maxy0
    else:
        minx, maxx, miny, maxy = min(minx0,minx), max(maxx0,maxx), min(miny0,miny), max(maxy0,maxy)

for rp in rps:
    # global inundation file
    inun_globe_file = os.path.join(path, 'Hazard/deltares_2021/Interpolation/inuncoast_historical_nosub_hist_rp{:04d}_0.nc'.format(rp))
    
    # Open inundation raster
    inun_globe = xr.open_dataset(inun_globe_file)
    
    lats0 = inun_globe['lat'].values
    lons0 = inun_globe['lon'].values
    
    # latitude lower and upper index
    latli = np.argmin(np.abs(lats0 - miny))
    latui = np.argmin(np.abs(lats0 - maxy))
    
    # longitude lower and upper index
    lonli = np.argmin(np.abs(lons0 - minx))
    lonui = np.argmin(np.abs(lons0 - maxx))
    
    # variable 
    inun_clip = inun_globe.inun.values
    lons = lons0[lonli:lonui+1]
    lats = lats0[latli:latui+1]
    
    if (len(inun_clip.shape) == 3):
        inun_clip = inun_clip[0, latli:latui+1, lonli:lonui+1]
    else:
        inun_clip = inun_clip[latli:latui+1, lonli:lonui+1]
            
    # Make nc file
    ds = inun_globe.copy()
    ds = ds.drop_dims('time')
    ds['lat']=lats[:]
    ds['lon']=lons[:]
    ds['inun']=(['lat','lon'], inun_clip)
    ds.attrs['title'] = basin+' inundation file'
    ds.attrs['Author'] = 'Huazhi Li'
    ds.attrs['institution'] = 'VU Amsterdam'
    
    del ds.attrs['config_file']
    del ds.attrs['project']
    del ds.attrs['references']
    del ds.attrs['history']
    
    # Save to nc 
    inun_output = os.path.join(path,'Hazard/Spatial_inun/GLOFRIS/'+basin+'/cluster{:d}'.format(cluster))
    if not os.path.isdir(inun_output): os.makedirs(inun_output)
    outdir = os.path.join(inun_output, 'Inun_rp{:04d}.nc'.format(rp))
    ds.to_netcdf(path=outdir,mode='w')
    ds.close()

    

