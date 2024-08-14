# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 20:23:17 2024

@author: hli490
"""

import os
import sys
import numpy as np
import pandas as pd

def agg_dmg(geo_unit, df, geogunit):
    df_agg =  pd.DataFrame(index=df.columns[1:],columns=['FID']+df['year'].to_list())
    df_agg['FID'] = geo_unit
    df_agg[df_agg.columns[1:]] = df.iloc[:,1:].T
    df_agg = df_agg.groupby('FID').sum()
    
    if geogunit==107:
        df_agg_0 = df_agg.transpose().drop(columns=-1)
    else:
        df_agg_0 = df_agg.transpose()
    df_agg_final = pd.DataFrame(index=range(len(df_agg_0)), columns = ['year']+df_agg_0.columns.to_list())
    df_agg_final['year'] = np.arange(1,10001)
    df_agg_final[df_agg_final.columns[1:]] = np.array(df_agg_0.iloc[:,:])
    
    return df_agg_final

# Set up directory, basin, cluster
path = 'c:/Users/hli490/Desktop/Paper2'
# path =  '/projects/0/FWC2/Spatial_dependence/Paper2'
iso_107 = 2020
# basin = 'NEA'
# cluster = 1

# Station and geounit files
station_folder = os.path.join(path,'Geogunit/Station/station_info_107.csv')
station_geo = pd.read_pickle(os.path.join(path,'Impact/Input/station_geo107_globe_unique.pkl'))
df_sts = pd.read_csv(station_folder)

# geounit
station_geo.loc[station_geo['FID']==iso_107,'station']

# Get the basins and clusters
basins = np.unique(df_sts.loc[df_sts['iso107']==iso_107,'basin'])

units = np.unique(station_geo['FID'])

annual_damage_nopro_folder = os.path.join(path,'Impact/Annual_damage/nopro/'+basin+'/annual_total_damage_cluster{:d}'.format(cluster))
annual_damage_pro_folder = os.path.join(path,'Impact/Annual_damage/pro/'+basin+'/annual_total_damage_cluster{:d}'.format(cluster))

# Read files
annual_dmg_nopro = pd.read_pickle(annual_damage_nopro_folder)
annual_dmg_pro = pd.read_pickle(annual_damage_pro_folder)

geog107_0 = pd.read_csv(os.path.join(path,'Geogunit/geogunit_107_0_list.csv'),encoding = "ISO-8859-1")

# geogunit 107
geo_unit = np.array([]).astype(float)
for i in range(1,len(annual_dmg_nopro.columns)):
    st = annual_dmg_nopro.columns[i]
    unit = np.array(station_geo.loc[station_geo['station']==st,'FID'])
    if len(unit)==0: geo_unit = np.append(geo_unit, -1)
    else: geo_unit = np.append(geo_unit, unit)

df_107_nopro = agg_dmg(geo_unit = geo_unit, df = annual_dmg_nopro, geogunit=107)
df_107_pro = agg_dmg(geo_unit = geo_unit, df = annual_dmg_pro, geogunit=107)

# save
out_dir = os.path.join(path,'Impact/Geogunit107/Total//nopro/'+basin)
if not os.path.isdir(out_dir): os.mkdir(out_dir)
df_107_nopro.to_csv(os.path.join(out_dir,'annual_dmg_cluster{:d}.csv'.format(cluster)),index=False)

out_dir = os.path.join(path,'Impact/Geogunit107/Total//pro/'+basin)
if not os.path.isdir(out_dir): os.mkdir(out_dir)
df_107_pro.to_csv(os.path.join(out_dir,'annual_dmg_cluster{:d}.csv'.format(cluster)),index=False)


# geogunit 0
geo_unit = np.array([]).astype(float)
for i in range(1,len(df_107_nopro.columns)):
    geog107 = df_107_nopro.columns[i]
    unit = np.array(geog107_0.loc[geog107_0['FID']==geog107,'FID_geog0'])
    geo_unit = np.append(geo_unit, unit)

df_0_nopro = agg_dmg(geo_unit = geo_unit, df = df_107_nopro, geogunit=0)
df_0_pro = agg_dmg(geo_unit = geo_unit, df = df_107_pro, geogunit=0)

# save
out_dir = os.path.join(path,'Impact/Geogunit0/Total//nopro/'+basin)
if not os.path.isdir(out_dir): os.mkdir(out_dir)
df_0_nopro.to_csv(os.path.join(out_dir,'annual_dmg_cluster{:d}.csv'.format(cluster)),index=False)

out_dir = os.path.join(path,'Impact/Geogunit0/Total//pro/'+basin)
if not os.path.isdir(out_dir): os.mkdir(out_dir)
df_0_pro.to_csv(os.path.join(out_dir,'annual_dmg_cluster{:d}.csv'.format(cluster)),index=False)

# # geogunit 0
# geog107 = pd.read_csv(r'C:\Users\hli490\Desktop\Paper2\Geogunit\geogunit_107_list.csv',encoding = "ISO-8859-1")
# geog0 = pd.read_csv(r'C:\Users\hli490\Desktop\Paper2\Geogunit\geogunit_0_list.csv',encoding = "ISO-8859-1")

# geog107_0 = geog107.copy()
# geog107_0['FID_geog0'] = np.nan
# for i in range(len(geog107_0)):
#     iso = geog107_0.iloc[i,1]
#     geog107_0.iloc[i,2] = geog0.loc[geog0['ISO']==iso,'FID'].iloc[0]
    
# geog107_0.to_csv(r'C:\Users\hli490\Desktop\Paper2\Geogunit\geogunit_107_0_list.csv',index=False)
