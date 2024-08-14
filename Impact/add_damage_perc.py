# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 13:44:11 2024

@author: hli490
"""

import os
import sys
import numpy as np
import pandas as pd

def dmg_split(dmg,rp,perc,st_new,pro,pro_flag=0):
    
    rps=np.array([2,5,10,25,50,100,250,500,1000])
    df = pd.DataFrame(index=range(len(dmg)), columns=st_new)
    df.iloc[:,:] = 0.0
    
    idx = np.where(dmg>0)[0]
    dmg0 = np.array(dmg[idx].astype(float))
    rp0 = np.array(rp[idx].astype(float))
    rp0[rp0>1000] = 1000
    
    for i in range(len(idx)):
        idx0 = idx[i]
        rp_emp = rps[np.argmin(abs(rp0[i] -rps))]
        perc0 = np.array(perc[str(rp_emp)])
        df.iloc[idx0,:] = perc0*dmg0[i]
    
    # add protection standard
    if pro_flag == 1: 
        for i in range(len(df.columns)):
            st = df.columns[i].astype(float)
            pro_st = pro.loc[pro['gtsm']==st,'pro'].values[0]
            idx_pro = idx[rp0<=pro_st]
            if len(idx_pro)!=0: df.iloc[idx_pro,i] = 0.0
     
    return df


# Set up directory, basin, cluster
path = 'c:/Users/hli490/Desktop/Paper2'
# path =  '/projects/0/FWC2/Spatial_dependence/Paper2'
# basin = str(sys.argv[1])
# cluster = int(sys.argv[2])
# basin = 'NEA'
# cluster = 1

# Locations of damage calculation files 
# event_folder = os.path.join(path,'Eventset/'+basin+'/rp_event_set_cluster{:d}.pkl'.format(cluster))
# event_damage_nopro_folder = os.path.join(path,'Impact/station/nopro/'+basin+'/damage_event_per_station_cluster{:d}.pkl'.format(cluster))
# event_damage_pro_folder = os.path.join(path,'Impact/station/pro/'+basin+'/damage_event_per_station_cluster{:d}.pkl'.format(cluster))
station_geo_folder = os.path.join(path,'Impact/Input/station_id_new_geo107_globe_unique_pro.csv')

event_damage_pro_folder = r'C:\Users\hli490\Desktop\Paper2\test\damage_event_per_station_cluster1.pkl'
event_folder =  r'C:\Users\hli490\Desktop\Paper2\test\rp_event_set_cluster1.pkl'

# Read files
event_rp = pd.read_pickle(event_folder)
# df_dmg_nopro = pd.read_pickle(event_damage_nopro_folder)
df_dmg_pro = pd.read_pickle(event_damage_pro_folder)
dmg_perc = pd.read_pickle(os.path.join(path,'Impact/Input/dmg_perc_globe.pkl'))
station_geo = pd.read_csv(station_geo_folder)

# Get the stations that are across several geogunits
sts_perc = np.array(dmg_perc['Station'])
sts_dmg = np.array(df_dmg_pro.columns.astype(float))
sts_intersect = np.intersect1d(sts_dmg, sts_perc)

# if no station is across more than one units
if len(sts_intersect)==0: 
    # out_dir = os.path.join(path,'Impact/station_unique/nopro/'+basin)
    # if not os.path.isdir(out_dir): os.mkdir(out_dir)
    # df_dmg_nopro.to_pickle(os.path.join(out_dir,'damage_event_per_station_cluster{:d}.pkl'.format(cluster)))
    out_dir = os.path.join(path,'Impact/station_unique/pro/'+basin)
    if not os.path.isdir(out_dir): os.mkdir(out_dir)
    df_dmg_pro.to_pickle(os.path.join(out_dir,'damage_event_per_station_cluster{:d}.pkl'.format(cluster)))
    sys.exit()

# split into unique station
sts_new = np.array([]).astype(float)
for i in range(len(sts_dmg)):
    st = sts_dmg[i]
    if not st in sts_perc:
        sts_new = np.append(sts_new, st)
    else:
        n_units = len(dmg_perc.loc[dmg_perc['Station']==st,'FID'].iloc[0])
        sts_new = np.append(sts_new, st+np.arange(0,n_units,1)/100)
   
non_inter = sts_dmg[~np.isin(sts_dmg,sts_intersect)]
# dmg_nopro_new = pd.DataFrame(index=range(len(df_dmg_nopro)),columns=sts_new)
# dmg_nopro_new[non_inter] = df_dmg_nopro[non_inter.astype(int).astype(str)]

dmg_pro_new = pd.DataFrame(index=range(len(df_dmg_pro)),columns=sts_new)
dmg_pro_new[non_inter] = df_dmg_pro[non_inter.astype(int).astype(str)]


for st in sts_intersect:
    st = 3960
    idx = np.where(sts_dmg == st)[0][0]
    # dmg_nopro = df_dmg_nopro.iloc[:,idx]
    dmg_pro = df_dmg_pro.iloc[:,idx]
    print(sum(dmg_pro>0))
    rp = event_rp.iloc[:,idx]
    print(sum(rp>1000))
    perc = dmg_perc.loc[dmg_perc['Station']==st].iloc[0]
    n_units = len(dmg_perc.loc[dmg_perc['Station']==st,'FID'].iloc[0])
    st_new = st+np.arange(0,n_units,1)/100
    
    # df_nopro = dmg_split(dmg=dmg_nopro, rp=rp, perc=perc, st_new=st_new)
    # dmg_nopro_new[st_new] = df_nopro.iloc[:,:]
    
    df_pro = dmg_split(dmg=dmg_pro, rp=rp, perc=perc, st_new=st_new, pro=station_geo, pro_flag=1)
    dmg_pro_new[st_new] = df_pro.iloc[:,:]


# saving
# out_dir = os.path.join(path,'Impact/station_unique/nopro/'+basin)
# if not os.path.isdir(out_dir): os.mkdir(out_dir)
# dmg_nopro_new.to_pickle(os.path.join(out_dir,'damage_event_per_station_cluster{:d}.pkl'.format(cluster)))

out_dir = os.path.join(path,'Impact/station_unique/pro/'+basin)
if not os.path.isdir(out_dir): os.mkdir(out_dir)
dmg_pro_new.to_pickle(os.path.join(out_dir,'damage_event_per_station_cluster{:d}.pkl'.format(cluster)))   
