# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 09:33:30 2024

This script calculates the annual damages of the stochastic events in two manner:
    1) Annual total damages
    2) Annual maxima damages (per station): this is for avoiding double-counting in temporally clustered events

@author: hli490
"""

import os
import sys
import numpy as np
import pandas as pd

def annual_damage(data,yr_index,method,yr_sim=10000):
    
    # add the year information to the damage data
    data0 = data.copy()
    data0['year'] = yr_index
    
    # damage variable storing annual damages
    res = pd.DataFrame(index=range(yr_sim),columns=['year']+list(data.columns))
    res['year'] = np.arange(1,yr_sim+1)
    res[res.columns[1:]] = 0.0
    
    if method == 'total':
        res0 = data0.groupby('year').sum()
        idx = np.array(res0.index.astype(int))-1
        res.iloc[idx,1:] = np.array(res0.iloc[:,:])
    elif method == 'max':
        res0 = data0.groupby('year').max()
        idx = np.array(res0.index.astype(int))-1
        res.iloc[idx,1:] = np.array(res0.iloc[:,:])
    
    return res


# Set up directory, basin, cluster
path =  '/projects/0/FWC2/Spatial_dependence/Paper2'
basin = str(sys.argv[1])
cluster = int(sys.argv[2])

# Locations of damage calculation files 
yr_index_folder = os.path.join(path,'Eventset/'+basin+'/year_event_set_cluster{:d}.csv'.format(cluster))
event_damage_nopro_folder = os.path.join(path,'Impact/station_unique/nopro/'+basin+'/damage_event_per_station_cluster{:d}.pkl'.format(cluster))
event_damage_pro_folder = os.path.join(path,'Impact/station_unique/pro/'+basin+'/damage_event_per_station_cluster{:d}.pkl'.format(cluster))

# Read damage calculation files
yr_index = pd.read_csv(yr_index_folder)
df_dmg_nopro = pd.read_pickle(event_damage_nopro_folder)
df_dmg_pro = pd.read_pickle(event_damage_pro_folder)

# Calculate annual damages
# annual total damages
dmg_total_nopro = annual_damage(data=df_dmg_nopro, yr_index=yr_index['year'],method='total')
dmg_total_pro = annual_damage(data=df_dmg_pro, yr_index=yr_index['year'],method='total')

# annual maximum damages
dmg_max_nopro = annual_damage(data=df_dmg_nopro, yr_index=yr_index['year'],method='max')
dmg_max_pro = annual_damage(data=df_dmg_pro, yr_index=yr_index['year'],method='max')

# Saving files
out_dir = os.path.join(path,'Impact/Annual_damage/nopro/{:s}'.format(basin))
if not os.path.isdir(out_dir): os.mkdir(out_dir)
dmg_total_nopro.to_pickle(os.path.join(out_dir,'annual_total_damage_cluster{:d}.pkl'.format(cluster)))
dmg_max_nopro.to_pickle(os.path.join(out_dir,'annual_max_damage_cluster{:d}.pkl'.format(cluster)))

out_dir = os.path.join(path,'Impact/Annual_damage/pro/{:s}'.format(basin))
if not os.path.isdir(out_dir): os.mkdir(out_dir)
dmg_total_pro.to_pickle(os.path.join(out_dir,'annual_total_damage_cluster{:d}'.format(cluster)))
dmg_max_pro.to_pickle(os.path.join(out_dir,'annual_max_damage_cluster{:d}'.format(cluster)))

# the number of events that caused damages
# df_nr_nopro = np.sign(df_dmg_nopro).astype(int)
# df_nr_pro = np.sign(df_dmg_pro).astype(int)
# nr_nopro = annual_damage(data=df_nr_nopro, yr_index=yr_index['year'],method='total')
# nr_pro = annual_damage(data=df_nr_pro, yr_index=yr_index['year'],method='total')



