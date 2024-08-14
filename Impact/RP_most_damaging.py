# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 13:39:50 2024

@author: hli490
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Path and folder location
# path = r'C:\Users\hli490\Desktop\Paper2'
path =  '/projects/0/FWC2/Spatial_dependence/Paper2'
basin='NEA'
# Read files
geog0 = pd.read_csv(os.path.join(path,'Geogunit/geogunit_0_list.csv'),encoding = 'unicode_escape')
df_iso = pd.read_csv(os.path.join(path,'Geogunit/Station/station_info.csv'))
num_cluster = pd.read_csv(os.path.join(path,'Eventset/basin_cluster_num.csv'))
iso_cluster = df_iso[df_iso['basin']==basin]
iso_avail = np.unique(iso_cluster['iso'])
dmg_gloris = pd.read_csv(os.path.join(path,'Risk/Nopro/glofris/glofris_Built-up_Area_geogunit0.csv'))

df=pd.DataFrame(index=range(len(iso_avail)),columns=['iso','cluster','eventid','dmg_event','dmg_1000yr'])
df['iso']=iso_avail

for ttt in range(len(iso_avail)):
    # get the basins and clusters where this particular country crosses
    iso = iso_avail[ttt]
    iso_in = df_iso[df_iso['iso']==iso].reset_index(drop=True)
    basins = np.unique(iso_in['basin'])
    
    for basin in basins:
        iso_in_basin = iso_in[iso_in['basin']==basin].reset_index(drop=True)
        clusters = np.unique(iso_in_basin['cluster'])
        if basin == 'NA_basin': cluster_link = np.array(pd.read_csv(os.path.join(path,'Hazard/Cluster_link/NA.csv')))
        else: cluster_link = np.array(pd.read_csv(os.path.join(path,'Hazard/Cluster_link/{:s}.csv'.format(basin))))
        
        for i in range(len(clusters)):
            cluster = int(clusters[i])
            cluster0 = cluster_link[cluster-1,~np.isnan(cluster_link[cluster-1,:])] # all linked clusters
            if i==0: cluster_all = cluster0
            else: cluster_all = np.concatenate((cluster_all , cluster0))
            
        cluster_all = np.sort(np.unique(cluster_all))
        imp_max = 0
        
        for j in range(len(cluster_all)):
            cluster = cluster_all[j]
            impact = pd.read_csv(os.path.join(path,'Risk/Nopro/geogunit_0/'+basin+'/damage_cluster{:d}.csv'.format(int(cluster))))
            imp0 = impact.loc[:,str(float(iso))]
            
            if max(imp0)>imp_max:
                idx = imp0.idxmax()
                cluster_max = cluster
                imp_max = max(imp0)
        
        df.iloc[ttt,1] = cluster_max
        df.iloc[ttt,2] = idx
        df.iloc[ttt,3] = imp_max
        df.iloc[ttt,4] = dmg_gloris.iloc[iso,10]
        
        if basin == 'NA_basin': eventset = pd.read_pickle(os.path.join(path,'Eventset/NA/rp_event_set_cluster{:d}.pkl'.format(cluster_max)))
        else: eventset =  pd.read_pickle(os.path.join(path,'Eventset/'+basin+'/rp_event_set_cluster{:d}.pkl'.format(cluster_max)))

        stations = np.array(iso_in['station']).astype(int)
        ds_event = (eventset.columns.values).astype(int)
        event_cluster = eventset.loc[:,(ds_event[np.nonzero(np.in1d(ds_event,stations))[0]]).astype(str)]
        
        data = event_cluster .iloc[idx,:]
        bins = [0,2,5,10,25,50,100,250,500,1000,10000]
        
        hist, bin_edges = np.histogram(data,bins)
        hist = hist/len(stations)
        
        #plotting
        fig,ax = plt.subplots(figsize=(10, 7))

        # Plot the histogram heights against integers on the x axis
        ax.bar(range(len(hist)),hist,width=1,edgecolor='black', linestyle="-",linewidth=0.5) 
        
        # Set the ticks to the middle of the bars
        ax.set_xticks([-0.5+t1 for t1,t2 in enumerate(hist)])
        ax.set_xlim([-1,10])
        ax.grid(color='gray', linestyle='dashed', linewidth=0.5)
        ax.set_axisbelow(True)
        # ax.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.5)
        # Set the xticklabels to a string that tells us what the bin edges were
        ax.set_xticklabels(['0','2','5','10','250','50','100','250','500','1000'])
        ax.set_xlabel('Return period (yr)',fontsize=11)
        ax.set_ylabel('Portion of stations (-)',fontsize=11)
        ax.set_title('Distribution of station return periods of the most damaging event',fontsize=14)
        out_dir = os.path.join(path,'Risk/RP_most_damaging/{:d}.png'.format(iso))
        plt.savefig(out_dir,format='png',dpi=600,bbox_inches='tight') 
        
df.to_csv(os.path.join(path,'Risk/RP_most_damaging/rp_most_damaging.csv'),index=False)



         