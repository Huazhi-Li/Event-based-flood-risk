# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 14:56:22 2024

@author: hli490
"""
import os
import sys
import numpy as np
import pandas as pd
import glob 
import re

def agg_year(df1,df2,variable='dmg',method='max'):
    out = pd.concat([df1,df2],axis=0)
    if method == 'max':
        res = out.groupby('year')[variable].max().astype(float)
    if method == 'total':
        res = out.groupby('year')[variable].sum().astype(float)
    return np.array(res)
        
# Set up directory, basin, cluster
path = 'c:/Users/hli490/Desktop/Paper2'
# path =  '/projects/0/FWC2/Spatial_dependence/Paper2'

# basin = 'NEA'
# cluster = 1

# Station and geounit files
station_folder = os.path.join(path,'Geogunit/Station/station_info.csv')
station_geo = pd.read_csv(os.path.join(path,'Impact/Input/station_id_new_geo107_globe_unique.csv'))

df_sts = pd.read_csv(station_folder)
df_dmg_globe = pd.read_csv(os.path.join(path,'Impact/Input/damage_per_station_globe.csv'))

# i = int(sys.argv[1])

i=0 
id_dmg = np.array([])
for i in range(len(station_geo)):
    # get the station name
    st = station_geo.iloc[i,0] # unique station id, unique station id is divided from the gtsm station id
    st_int = int(st) # gtsm station id
    
    # if no damages are caused from this station, skip this loop
    dmg_all = np.array(df_dmg_globe.loc[df_dmg_globe['station']==st_int,:])[0,:]
    
    if sum(dmg_all[1:])==0: continue
    else: 
        id_new = np.where(station_geo['gtsm'] == st)[0][0]
        id_dmg = np.append(id_dmg, id_new)

file = '/Impact/Annual_damage/nopro/annual_max_damage_station16420.csv'
regex = re.compile(r'\d+')
id_exist = np.array([])
for file in glob.glob(os.path.join(path+'/Impact/Annual_damage/nopro','*.csv')):
    id_exist = np.append(id_exist, (np.array(regex.findall(file)).astype(int))[2])
        
non_station = np.setxor1d(id_dmg,id_exist)

df = pd.DataFrame({'station': non_station})
df.to_csv(r'c:\Users\hli490\Desktop\Paper2\Impact\sts_left.csv',index=False)

basins = []
clusters = np.array([])
for st in non_station:
    st_int = int(st)
    basin = df_sts.loc[df_sts['station']==st_int,'basin'].iloc[0]
    cluster = df_sts.loc[df_sts['station']==st_int,'cluster'].iloc[0]
    basins.append(basin)
    
    if basin == 'SA': 
        clusters = np.append(clusters, cluster)
        print(st)



# get the basin and cluster information of this station
basin = df_sts.loc[df_sts['station']==st_int,'basin'].iloc[0]
cluster = df_sts.loc[df_sts['station']==st_int,'cluster'].iloc[0]
if basin == 'NA_basin': basin = 'NA'

# get all the clusters that might have events affecting this station
cluster_link_loc = os.path.join(path,'Hazard/Cluster_link/{:s}.csv'.format(basin))
cluster_link = np.array(pd.read_csv(cluster_link_loc))
# cluster_all = cluster_link[cluster-1,~np.isnan(cluster_link[cluster-1,:])]
cluster_all = np.where(cluster_link==cluster)[0]+1

# annual damage
# annual_damage = pd.DataFrame(data=dict(year=np.arange(1,10001),dmg = np.zeros(10000)))
annual_damage = np.zeros(10000)
cluster_all = [1,2]

# loop over these clusters
for j in range(len(cluster_all)):
    
    cluster_id = cluster_all[j]
    
    # the year index for this particular cluster
    yr_index = pd.read_csv(os.path.join(path,'Eventset/'+basin+'/year_event_set_cluster{:d}.csv'.format(cluster_id)))
    df_dmg_nopro = pd.read_pickle(os.path.join(path,'Impact/station_unique/nopro/'+basin+'/damage_event_per_station_cluster{:d}.pkl'.format(cluster_id)))
    df_dmg_nopro.columns = [float(col) for col in df_dmg_nopro.columns]
    # df_dmg_pro = pd.read_pickle(os.path.join(path,'Impact/station_unique/pro/'+basin+'/damage_event_per_station_cluster{:d}.pkl'.format(cluster_id)))

    
    # create the dmg df
    df = pd.DataFrame(data=dict(year = yr_index['year'], dmg =  df_dmg_nopro[st]))
    out = df.groupby('year')['dmg'].max()
    dmg = np.zeros(10000)
    dmg[out.index.astype(int)-1] = np.array(out).astype(float)
    annual_damage = np.maximum(annual_damage, dmg)    


year = np.arange(1,10000+1)
data = {
        "year":{"dims": ("year"), "data": year},
        'dmg':{"dims": ("year"), "data": annual_damage}
        }
    
data = xr.Dataset.from_dict(data)
out_dir = os.path.join(path,'Impact/Event_damage/nopro')
data.to_netcdf(path=os.path.join(out_dir,'annual_max_dmg_gtsm{:07d}.nc'.format(int(st*100))),mode='w')

del data,dmg

a = [1,3,5,7,9,2]
b = [2,2,4,8,7,5]

np.maximum(a,b)
    
