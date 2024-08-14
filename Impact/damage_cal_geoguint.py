# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 11:41:41 2024

@author: hli490
"""

import os
import sys
import numpy as np
import pandas as pd

        
# Set up directory, basin, cluster
path = 'c:/Users/hli490/Desktop/Paper2'
# path =  '/projects/0/FWC2/Spatial_dependence/Paper2'

# Station and geounit files
df_sts = pd.read_csv(os.path.join(path,'Geogunit/Station/station_info.csv'))
station_geo = pd.read_pickle(os.path.join(path,'Impact/Input/station_geo107_globe_unique.pkl'))
df_dmg_globe = pd.read_csv(os.path.join(path,'Impact/Input/damage_per_station_globe.csv'))
dmg_perc = pd.read_pickle(os.path.join(path,'Impact/Input/dmg_perc_globe.pkl'))

fid = 1283

sts_in = np.array(station_geo.loc[station_geo['FID']==fid,'station'])

