"""
Assign inundation cells to 'nearest' stations. There must be a water pathway from the cell to its assigned station.
This script was built upon from the coastal inundation model for the Aqueduct Coastal Flooding project. 

@author: Huazhi Li (huazhi.li@vu.nl)
"""

# Import modules
import sys
import os.path
import numpy as np
import pandas as pd
import rasterio
import xarray as xr
from rasterio.mask import mask
from shapely.geometry import box
import geopandas as gpd
from fiona.crs import from_epsg
import pycrs
from rasterio.transform import Affine
from rasterio.fill import fillnodata
from numba import njit, prange
import matplotlib.pyplot as plt
from collections import deque
from affine import identity as IDENTITY
from rasterio.plot import show

def degree_metres_y(lat):
    """ "returns the verical length of a degree in metres at
    a given latitude."""
    m1 = 111132.92  # latitude calculation term 1
    m2 = -559.82  # latitude calculation term 2
    m3 = 1.175  # latitude calculation term 3
    m4 = -0.0023  # latitude calculation term 4
    # # Calculate the length of a degree of latitude and longitude in meters
    radlat = np.radians(lat)
    latlen = (
        m1
        + (m2 * np.cos(2.0 * radlat))
        + (m3 * np.cos(4.0 * radlat))
        + (m4 * np.cos(6.0 * radlat))
    )
    return latlen


def degree_metres_x(lat):
    """ "returns the horizontal length of a degree in metres at
    a given latitude."""
    p1 = 111412.84  # longitude calculation term 1
    p2 = -93.5  # longitude calculation term 2
    p3 = 0.118  # longitude calculation term 3
    # # Calculate the length of a degree of latitude and longitude in meters
    radlat = np.radians(lat)
    longlen = (
        (p1 * np.cos(radlat))
        + (p2 * np.cos(3.0 * radlat))
        + (p3 * np.cos(5.0 * radlat))
    )
    return longlen

def spread2d(id_s,obs, msk, frc, latlon=False, transform=IDENTITY):
    nrow, ncol = obs.shape
    xres, yres, north = transform[0], transform[4], transform[5]
    dx, dy = xres, yres
    
    out = obs.copy()
    src = np.full(obs.shape, -1, dtype=np.int32)   # linear index of source
    dst = np.full(obs.shape, -1, dtype=np.float32) # distance from source
    nxt = obs > 0    
    nxt1 = nxt.copy()
    
    for r, c in zip(*np.where(nxt)):
        src[r, c] = id_s[r,c]
        dst[r, c] = 0
        
    while np.any(nxt):
        for  r in prange(nrow):
            for c in prange(ncol):
                if nxt[r, c]:
                    nxt1[r, c] = False
                    i0 = src[r, c]
                    d0 = dst[r, c]
                    f0 = frc[r, c]
                    if latlon:
                        lat = north + (r + 0.5) * yres
                        dy = degree_metres_y(lat)*yres
                        dx = degree_metres_x(lat)*xres
                    for dr in range(-1,2):
                        for dc in range(-1,2):
                            if dr == 0 and dc == 0:
                                continue
                            r1, c1 = r+dr, c+dc
                            if r1 < 0 or r1 >= nrow or c1 < 0 or c1 >= ncol or msk[r1,c1] != 1:
                                continue
                            dd = np.hypot(dx*dr, dy*dc)
                            d = d0 + dd*f0
                            if src[r1, c1] == -1 or d < dst[r1, c1]:
                                src[r1, c1] = i0
                                dst[r1, c1] = d
                                out[r1, c1] = obs.flat[i0]
                                nxt1[r1, c1] = True
        nxt = nxt1.copy()
    return out, src, dst
                        
def spread2d_pnts(id_s, rs, cs, zs, msk, frc=None, nodata=0, latlon=False, transform=IDENTITY):
    
    z = np.full(msk.shape, nodata, dtype=zs.dtype)
    z[rs, cs] = zs
    
    if frc is None:
        frc = np.ones_like(msk, dtype=np.float32)
    
    return spread2d(id_s, z, msk, frc, latlon=latlon, transform=transform)

#####################################################################################################################################################
## MAIN SCRIPT ##

# Make the working diretory
path = '/projects/0/FWC2/Spatial_dependence/Paper2/Hazard/'
# path = 'c:/Users/hli490/Desktop/Spatial_dependence/Paper2/Hazard/'
basin = str(sys.argv[1])
cluster = int(sys.argv[2])
rp_id = int(sys.argv[3])
# basin = 'NEA'
# cluster = 1
rps = [2,5,10,25,50,100,250,500,1000]

# get the boundary coordinates of the potentially affected clusters
extent = pd.read_csv(os.path.join(path, 'Extent/'+basin+'.csv'))
cluster_link = np.array(pd.read_csv(os.path.join(path, 'Cluster_link/'+basin+'.csv')))
cluster_all = cluster_link[cluster-1,~np.isnan(cluster_link[cluster-1,:])] # all linked clusters

# get the extent coordinates
for i in range(len(cluster_all)):
    cluster0 = int(cluster_all[i])-1
    minx0, maxx0, miny0, maxy0 = extent.iloc[cluster0,1], extent.iloc[cluster0,2], extent.iloc[cluster0,3], extent.iloc[cluster0,4]
    
    if i==0: 
        minx, maxx, miny, maxy = minx0, maxx0, miny0, maxy0
    else:
        minx, maxx, miny, maxy = min(minx0,minx), max(maxx0,maxx), min(miny0,miny), max(maxy0,maxy)

# source station file
src_basin_file = os.path.join(path, 'Source_station/src_station_rp1000_{:s}.nc'.format(basin))
src_basin = xr.open_dataset(src_basin_file)

# Open station file
df_station_all = pd.read_csv(os.path.join(path, 'Stations/'+basin+'.csv'))
df_station = (df_station_all.loc[df_station_all['cluster'].isin(cluster_all)]).reset_index()
df_station_else = (df_station_all.loc[~(df_station_all['cluster'].isin(cluster_all))]).reset_index()



for rp in rps[rp_id:]:

    # global inundation file
    inun_globe_file = os.path.join(path, 'deltares_2021/HIST/inuncoast_historical_nosub_hist_rp{:04d}_0.nc'.format(rp))
    
    # Open inundation raster
    inun_globe = xr.open_dataset(inun_globe_file)
    
    # Raster clipping
    # Mask layer extent, WGS84 coordinates
    
    lats0 = inun_globe['lat'].values
    lons0 = inun_globe['lon'].values
    
    # latitude lower and upper index
    latli = np.argmin(np.abs(lats0 - miny))
    latui = np.argmin(np.abs(lats0 - maxy))
    
    # longitude lower and upper index
    lonli = np.argmin(np.abs(lons0 - minx))
    lonui = np.argmin(np.abs(lons0 - maxx))
    
    # get the inundation
    inun = inun_globe.inun.values
    
    if (len(inun.shape) == 3):
        inun = inun[0, latli:latui+1, lonli:lonui+1]
    else:
        inun = inun[latli:latui+1, lonli:lonui+1]
        
    nodata =-9999
    inun=np.nan_to_num(inun, nan=nodata)
    
    # Make mask layer
    msk = inun!=0 # boolean array with the inundation and water body pixels showing Ture and the non-inundation pixels showing False
    
    # cells inundated by stations from other clusters are aslo marked False
    lats_basin = src_basin['lat'].values
    lons_basin = src_basin['lon'].values
    # latitude lower and upper index
    latli_basin = np.argmin(np.abs(lats_basin - miny))
    latui_basin = np.argmin(np.abs(lats_basin - maxy))
    # longitude lower and upper index
    lonli_basin = np.argmin(np.abs(lons_basin- minx))
    lonui_basin = np.argmin(np.abs(lons_basin - maxx))
    
    np_src_clip = src_basin.src.values[latli_basin:latui_basin+1, lonli_basin:lonui_basin+1]
    station_all = np.unique(np_src_clip)
    station_known = np.append(np.array(df_station['id']),[-1,-9999])
    station_else = np.delete(station_all,np.isin(station_all,station_known))
    msk[np.isin(np_src_clip,station_else)]=False
    
    ## Make stations in the mask layer
    rs = np.zeros(len(df_station['id']),dtype=np.int32)
    cs = np.zeros(len(df_station['id']),dtype=np.int32)
    zs = np.zeros(len(df_station['id']),dtype=np.float32)
    
    lons = lons0[lonli:lonui+1]
    lats = lats0[latli:latui+1]
    
    for i in range(len(rs)):
        rs[i] = np.argmin(np.abs(lats - df_station['lat'][i]))
        cs[i] = np.argmin(np.abs(lons - df_station['lon'][i]))
        zs[i] = round(np.random.uniform(0, 5), 2)    
    
    
    id_s = np.full(msk.shape, nodata, dtype=np.int32)
    id_s[rs, cs] = df_station['id']
    
    out, src, dst = spread2d_pnts(id_s, rs, cs, zs, msk, nodata=nodata)
    
    # Make source station raster where every inundated pixel is assigned to its 'nearest' station; the pixel value is the station id
    src_station = src.copy()
    
    r_noninun, c_noninun = np.where(inun==0) # non-inundated cells, pixel value -1
    src_station[r_noninun, c_noninun] = -1
    r_ocean, c_ocean = np.where(inun<0)  # ocean cells, pixel value nodata
    src_station[r_ocean, c_ocean] = nodata
    
    # Make nc file
    ds = inun_globe.copy()
    ds = ds.drop_dims('time')
    ds['lat']=lats[:]
    ds['lon']=lons[:]
    ds['src']=(['lat','lon'], src_station)
    ds.attrs['title'] = 'Source station'
    ds.attrs['Author'] = 'Huazhi Li'
    ds.attrs['institution'] = 'VU Amsterdam'
    
    del ds.attrs['config_file']
    del ds.attrs['project']
    del ds.attrs['references']
    del ds.attrs['history']
    
    # Save to nc 
    src_output = os.path.join(path,'Source_station/'+basin+'/cluster{:d}'.format(cluster))
    if not os.path.isdir(src_output): os.makedirs(src_output)
    outdir = os.path.join(src_output, 'src_station_rp{:04d}.nc'.format(rp))
    ds.to_netcdf(path=outdir,mode='w')
    ds.close()

