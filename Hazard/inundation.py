# -*- coding: utf-8 -*-
"""
Creating the spatially-dependent inundation layer

@author: Huazhi Li (huazhi.li@vu.nl)
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

path = '/projects/0/FWC2/Spatial_dependence/Paper2/'

# path = r'c:\Users\hli490\Desktop\Paper2\\'
basin = str(sys.argv[1])
cluster = int(sys.argv[2])
i = int(sys.argv[3]) # event id

# Event folder
event_folder = os.path.join(path,'Eventset/'+basin+'/rp_event_set_cluster{:d}.pkl'.format(cluster))
thr_folder = os.path.join(path,'Eventset/'+basin+'/thr.csv')

# Read the threshold
thr = pd.read_csv(thr_folder)

# Read the events
event_rp =  pd.read_pickle(event_folder)
if event_rp.max(axis=1)[i]<=2: sys.exit() # skip the event when the largest rp < 2 yrs
event = (event_rp.iloc[i,:]).sort_values()  # extract rps for a single event
event = event[event>2]
if (len(event)==0): sys.exit()

# Other input files
hazard_layer = os.path.join(path,'Hazard/Spatial_inun/GLOFRIS/'+basin+'/cluster{:d}'.format(cluster))
source_station_folder = os.path.join(path,'Hazard/Source_station/'+basin+'/cluster{:d}'.format(cluster))

# Read the flood protection standard
#flopros = pd.read_csv(os.path.join(path,'Flopros/FLOPROS_{:s}.csv'.format(basin)))

# Initialize the inundation layer
max_rp = nearest_rp(max(event))
inun_init_file = xr.load_dataset(os.path.join(hazard_layer,'Inun_rp{:04d}.nc'.format(max_rp)))
inun_event0 = inun_init_file.inun.values

inun_event = inun_event0.copy()
nodata = -9999
inun_event[np.isnan(inun_event)] = nodata 
inun_event[np.where(inun_event>=0)] = 0     


prostand = 2 # protection standard of natural dyke; set to 2 years

# Assemble the spatiallt-dependent inundation layer
for ii in range(len(event)):
    id_s = int(event.index[ii]) # station id 
    s_cluster = thr.loc[thr['station']==id_s,'cluster'].values[0] # get the station cluster
    if basin=='NEA' and s_cluster==28: continue
    
    rp_station = event[ii] # get the station rp
    rp_station = nearest_rp(rp_station) # assign this rp to the nearest known rp
    
    if (rp_station <= prostand): continue # skip the station less than the local protection standard    
    src_station_nc = xr.load_dataset(os.path.join(source_station_folder,'src_station_rp{:04d}.nc'.format(rp_station)))
    src_station = src_station_nc.src.values
    inun_file = xr.load_dataset(os.path.join(hazard_layer,'Inun_rp{:04d}.nc'.format(rp_station)))
    inun = inun_file.inun.values
    r_station, c_station = np.where(src_station==int(id_s))
    inun_event[r_station, c_station] = inun[r_station, c_station]
    
    src_station_nc.close()
    inun_file.close()

# Make nc file if there are inundated cells
if inun_event.max()>0: 
    ds = inun_init_file.copy()
    # ds = ds.drop_dims('rp')
    ds['inun'] = (['lat','lon'], inun_event)
    ds.attrs['title'] = 'Spatially-dependent inundation file'
    ds.attrs['source'] = 'GLOFRIS'
    
    # Save to nc 
    fn = os.path.join(path, 'Hazard/Spatial_inun/'+basin+'/cluster{:d}'.format(cluster))
    if not os.path.isdir(fn): os.makedirs(fn)
    outdir = os.path.join(fn, 'Inun_event{:06d}.nc'.format(i))
    ds.to_netcdf(path=outdir,mode='w')
    ds.close()

    

    
    
    
    
    


