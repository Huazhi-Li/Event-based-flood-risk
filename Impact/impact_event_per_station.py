# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 13:28:03 2024

@author: hli490
"""

import os
import sys
import numpy as np
import pandas as pd

# Set up directory, basin, cluster
path = 'c:/Users/hli490/Desktop/Paper2'
# path =  '/projects/0/FWC2/Spatial_dependence/Paper2'
# basin = str(sys.argv[1])
# cluster = int(sys.argv[2])
basin = 'NEA'
cluster = 1

# Read damage calculation files
dmg_station = pd.read_csv(os.path.join(path,'Impact/Input/damage_per_station_globe.csv'))
flopros = pd.read_csv(os.path.join(path,'Flopros/FLOPROS_station.csv'))

# Event damage at the station scale
# Event folder
event_folder = os.path.join(path,'Eventset/'+basin+'/rp_event_set_cluster{:d}.pkl'.format(cluster))
event_rp =  pd.read_pickle(event_folder)

# Damage files with and without protection standard
df_dmg_nopro = pd.DataFrame(index=range(len(event_rp)),columns=event_rp.columns)
df_dmg_nopro.iloc[:,:] = 0.0
df_dmg_pro = df_dmg_nopro.copy()
rps = np.array([2,5,10,25,50,100,250,500,1000])
station_all = np.array(event_rp.columns.astype(int))

# for i in range(len(event_rp)):
for i in range(100):    
    # skip the event where no rp is above 2 years, i.e. the event does not cause any damages
    if event_rp.max(axis=1)[i]<=2: continue 

    # extract rps at each station
    event = np.array(event_rp.iloc[i,:])
    dmg_all = np.zeros(len(event))
    dmg_all_pro = np.zeros(len(event))
    
    # calculate damage for stations with rp>2
    for j in np.where(event>2)[0]:
        rp = event[j]
        station = station_all[j]
        idx = np.where(dmg_station['station']==station)[0][0]
        dmg_glofris = np.array(dmg_station.iloc[idx,1:])
        
        # without protection
        dmg_all[j] = np.interp(rp, rps, dmg_glofris)
        
        # with protection
        if rp<=flopros.iloc[idx,2]: continue
        dmg_all_pro[j] = np.interp(rp, rps, dmg_glofris)
    
    df_dmg_nopro.iloc[i,:] = dmg_all
    df_dmg_pro.iloc[i,:] = dmg_all_pro
    
# out_dir = os.path.join(path,'Impact/station/nopro/'+basin)
# if not os.path.isdir(out_dir): os.mkdir(out_dir)
# df_dmg_nopro.to_pickle(os.path.join(out_dir,'damage_event_per_station_cluster{:d}.csv'.format(cluster)))

# out_dir = os.path.join(path,'Impact/station/pro/'+basin)
# if not os.path.isdir(out_dir): os.mkdir(out_dir)
# df_dmg_pro.to_pickle(os.path.join(out_dir,'damage_event_per_station_cluster{:d}.csv'.format(cluster)))






















