# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 14:09:07 2024

@author: hli490
"""

import pandas as pd
import numpy as np
import xarray as xr
import os
import os.path
import sys

def dmg2rp(dmg, t=10000):
    """Function to translate water levels to return periods using Weibull's plotting formula"""
    dmg = dmg.sort_values(ascending=0)
    rank = dmg.rank(ascending=0)
    m = len(dmg) 
    p = (rank*m/t/(m+1))
    rp = 1/p
    return rp

path = '/projects/0/FWC2/Spatial_dependence/Paper2/'

# Station and geounit files
station_folder = os.path.join(path,'Geogunit/Station/station_info.csv')
station_geo = pd.read_pickle(os.path.join(path,'Impact/Input/station_geo107_globe_unique.pkl'))
df_sts = pd.read_csv(station_folder)
df_dmg_globe = pd.read_csv(os.path.join(path,'Impact/Input/damage_per_station_globe.csv'))


i = int(sys.argv[1])

# get the station name
st = station_geo.iloc[i,0] # unique station id, unique station id is divided from the gtsm station id
st_int = int(st) # gtsm station id

# if no damages are caused from this station, skip this loop
dmg_all = np.array(df_dmg_globe.loc[df_dmg_globe['station']==st_int,:])[0,:]
if sum(dmg_all[1:])==0: sys.exit()

# get the basin and cluster information of this station
basin = df_sts.loc[df_sts['station']==st_int,'basin'].iloc[0]
cluster = df_sts.loc[df_sts['station']==st_int,'cluster'].iloc[0]
if basin == 'NA_basin': basin = 'NA'

# get all the clusters that might have events affecting this station
cluster_link_loc = os.path.join(path,'Hazard/Cluster_link/{:s}.csv'.format(basin))
cluster_link = np.array(pd.read_csv(cluster_link_loc))
cluster_all = np.where(cluster_link==cluster)[0]+1


for cluster in cluster_all:
    df_dmg_nopro = pd.read_pickle(os.path.join(path,'Impact/station_unique/nopro/'+basin+'/damage_event_per_station_cluster{:d}.pkl'.format(cluster)))
    df_dmg_nopro.columns = [float(col) for col in df_dmg_nopro.columns]
    dmg_c = df_dmg_nopro[st]
    if not 'dmg' in locals():
        dmg = dmg_c
    else:
        dmg = pd.concat([dmg,dmg_c])

dmg = dmg.dropna()
dmg=dmg[dmg>0]
dmg = dmg.sort_values(ascending=0)
rp = dmg2rp(dmg)
t = range(len(rp))
data = {
        "t":{"dims": ("t"), "data": t},
        'dmg':{"dims": ("t"), "data": dmg},
        'rp':{"dims": ("t"), "data": rp}
        }
data = xr.Dataset.from_dict(data)
out_dir = os.path.join(path,'Impact/Event_damage/nopro')
data.to_netcdf(path=os.path.join(out_dir,'dmg_rp_gtsm{:07d}.nc'.format(int(st*100))),mode='w')

del data,dmg,rp



    
    
