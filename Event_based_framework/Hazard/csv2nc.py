# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 10:57:01 2023

@author: hli490
"""

import sys
import os.path
import numpy as np
import pandas as pd
import xarray as xr


def nearest_rp(rp):
    rps = np.array([2, 5, 10, 25, 50, 100, 250, 500, 1000])
    idx = np.argmin(abs(rp-rps))
    return rps[idx]

# path = '/projects/0/FWC2/Spatial_dependence/Paper2/'

path = r'c:\Users\hli490\Desktop\Spatial_dependence\Paper2\\'# Event folder

basin='NEA'
cluster = 1
event_folder = os.path.join(path,'Eventset/'+basin+'/rp_event_set_cluster{:d}.csv'.format(cluster))
event_rp =  pd.read_csv(event_folder)

event_rp.to_pickle(event_folder.replace('.csv','.pkl'))

df= pd.read_pickle(event_folder.replace('.csv','.pkl'))
data = xr.open_dataset(r'c:\Users\hli490\Desktop\Spatial_dependence\Paper1\RP\rp_gtsm06003.nc')
