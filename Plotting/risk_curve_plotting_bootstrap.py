# -*- coding: utf-8 -*-
"""
Created on Sun Apr  7 13:59:42 2024

@author: hli490
"""

# Import modules

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats, interpolate
import matplotlib.colors as mcolors
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.axes_divider import make_axes_area_auto_adjustable
import matplotlib.patheffects as pe

def get_return_periods(x, yr=10000):
        # weibull plotting formula a=0, b=1
        a = 0.0
        b = 1.0
        ranks = (len(x) + 1) - stats.rankdata(x, method="ordinal")
        lam = len(x)/yr
        freq = ((ranks - a) / (len(x) + b)) * lam
        rps = 1 / freq
        # return ranks, rps
        return rps
    
    

# Setting working directory
path = r'C:\Users\hli490\OneDrive - Vrije Universiteit Amsterdam\Desktop\Paper2'

# Getting iso information
geog0 = pd.read_csv(os.path.join(path,'Geogunit/geogunit_0_list.csv'),encoding = 'unicode_escape')
df_iso = pd.read_csv(os.path.join(path,'Geogunit/Station/station_info.csv'))
df_events = pd.read_csv(os.path.join(path,'Eventset/num_event_99th.csv'))

# original
c1 = '#0571b0' # blue
# c2='#999ACD' # purple
c2 = 'green'
# c2 ='BBFFFF'
# c2='#ECA47C' # orange
# c3='purple' # purple
# c2 = 'orangered'
c3 = '#e88035'
# c3 = 'orange'
c4 = '#ca0020' # red

# new
c1 = '#0571b0' # blue
c2 = c1
c3 = '#ca0020'
c4 = c3

# new
c1 = '#ca0020' # red RP-based no protection
c2 = '#0571b0' # blue Event-based no protection
c3 = 'purple' # purple RP-based  protection
c4 = '#35afd7' # light blue Event-based protection
# c4 = '#3fa671' # green Event-based protection

# fig=plt.figure(figsize=(4*4.13,2*4.13)) # width, height
# gs=gridspec.GridSpec(2,2,wspace=0.15,hspace=0.25)
fig=plt.figure(figsize=(4*4.13,4.5)) # width, height
gs=gridspec.GridSpec(1,2,wspace=0.15,hspace=0.25)
ax=[[] for _ in range(6)] 
plt.rcParams.update({'font.size':10})

# isos= [147,65,43,49,236,151,122,176]
# countries=['Myanmar','Ecuador','China','United States','Mozambique','Senegal','Philippines']
isos= [236,46,87,115]
countries=['United States','Republic of the Congo','Guinea-Bissau','Japan']
scale = 1e6
isos= [147,160,181]
countries=['Myanmar','New Caledonia','Puerto Rico']
# isos= [181]
# countries=['Puerto Rico']
# for idx,gs1,gs2,label,iso,country in zip(range(0,8),[0,0,1,1,2,2,3,3],[0,1,0,1,0,1,0,1],['a)','b)','c)','d)','e)','f)','g)','h)'],isos,countries):
# for idx,gs1,gs2,label,iso,country in zip(range(0,6),[0,0,1,1,2,2],[0,1,0,1,0,1],['a)','b)','c)','d)','e)','f)'],isos,countries):
for idx,gs1,gs2,label,iso,country in zip(range(0,2),[0,0],[0,1],['a)','b)'],isos,countries):
# for idx,gs1,gs2,label,iso,country in zip(range(0,1),[0,0],[0,1],['c)','b)'],isos,countries):
    iso_in = df_iso[df_iso['iso']==iso].reset_index(drop=True)

    basins = np.unique(iso_in['basin'])
    dmg_all = pd.read_csv(os.path.join(path, 'Impact/Glofris/glofris_Built-up_Area_geogunit0.csv'))
    
    file = os.path.join(path, 'Impact/complete_dependence_Built-up_Area_geogunit0.csv')
    if not os.path.isfile(file) or dmg_all['1000'][iso]==0: continue
    
    # No protection
    dmg_all_nopro = pd.read_csv(os.path.join(path, 'Impact/complete_dependence_Built-up_Area_geogunit0.csv'))
    annual_max_all_nopro = pd.read_csv(os.path.join(path,'Impact/Annual_damage/Geogunit0/nopro/annual_max_dmg_unit{:03d}.csv'.format(iso)))
    annual_max_nopro = annual_max_all_nopro[annual_max_all_nopro.columns[0:3]]
    dmg_unit_nopro = np.concatenate([[0,0],np.array(dmg_all_nopro.loc[dmg_all_nopro['FID']==iso])[0][2:]])
    rp = np.array([1,2,5,10,25,50,100,250,500,1000]).astype(int)
    dmg_cd_nopro = pd.DataFrame(data=dict(rp=rp,dmg=dmg_unit_nopro))
    
    # With protection
    dmg_all_pro = pd.read_csv(os.path.join(path, 'Impact/complete_dependence_Built-up_Area_geogunit0_pro.csv'))
    annual_max_all_pro = pd.read_csv(os.path.join(path,'Impact/Annual_damage/Geogunit0/pro/annual_max_dmg_unit{:03d}.csv'.format(iso)))
    annual_max_pro = annual_max_all_pro[annual_max_all_pro.columns[0:3]]
    dmg_cd_pro = pd.DataFrame(data=dict(rp= np.arange(1,1001),dmg=np.array(dmg_all_pro.loc[dmg_all_pro['FID']==iso])[0][1:]))
    
    rp_lowest_dmg = np.where(dmg_cd_pro['dmg']>0)[0][0].min()
    if not rp_lowest_dmg in rp.astype(int):
        rp_nearest_idx = np.argmin(abs(rp_lowest_dmg-rp))
        if rp[rp_nearest_idx]<rp_lowest_dmg: 
            rp_insert = np.arange(rp_lowest_dmg, rp[rp_nearest_idx+1])
        else:
            rp_insert = np.arange(rp_lowest_dmg, rp[rp_nearest_idx])
        rp_new = rp
        rp_new = np.sort(np.concatenate([rp,rp_insert]))
        dmg_cd_pro = dmg_cd_pro.loc[dmg_cd_pro['rp'].isin(rp_new),:]
      
    else:
        dmg_cd_pro = dmg_cd_pro.loc[dmg_cd_pro['rp'].isin(rp),:]
    
    # Bootstrapping
    rp_bootstrap = np.arange(1,1001)
    yr = 1000
    n = 100
    # event_bootstrap = pd.DataFrame(index=range(n),columns=rp_bootstrap)
    max_pro_bootstrap = pd.DataFrame(index=range(n),columns=rp_bootstrap)
    max_nopro_bootstrap = pd.DataFrame(index=range(n),columns=rp_bootstrap)
    
    for i in range(n):
    
        # No protection
        max_dmg_nopro = np.array(annual_max_nopro['dmg'].sample(n=int(len(annual_max_nopro)/(10000/yr)),replace=True))
        max_rp_nopro = get_return_periods(max_dmg_nopro, yr=1000)
        max_annual_nopro = pd.DataFrame(data=dict(dmg=max_dmg_nopro,rp=max_rp_nopro))
        max_annual_nopro = max_annual_nopro.sort_values(by='rp')
        max_nopro_bootstrap.iloc[i,:] = np.interp(rp_bootstrap,max_annual_nopro['rp'],max_annual_nopro['dmg'])
        
        # With protection
        max_dmg_pro = np.array(annual_max_pro['dmg'].sample(n=int(len(annual_max_pro)/(10000/yr)),replace=True))
        max_rp_pro = get_return_periods(max_dmg_pro, yr=1000)
        max_annual_pro = pd.DataFrame(data=dict(dmg=max_dmg_pro,rp=max_rp_pro))
        max_annual_pro = max_annual_pro.sort_values(by='rp')
        max_pro_bootstrap.iloc[i,:] = np.interp(rp_bootstrap,max_annual_pro['rp'],max_annual_pro['dmg'])
    
    
    # Medians and confidence intervals
    max_nopro_mean = max_nopro_bootstrap.mean().astype(float)
    max_nopro_ci_lower = max_nopro_bootstrap.quantile(0.05, numeric_only=False).astype(float)
    max_nopro_ci_upper = max_nopro_bootstrap.quantile(0.95, numeric_only=False).astype(float)
    
    max_pro_mean = max_pro_bootstrap.mean().astype(float)
    max_pro_ci_lower = max_pro_bootstrap.quantile(0.05, numeric_only=False).astype(float)
    max_pro_ci_upper = max_pro_bootstrap.quantile(0.95, numeric_only=False).astype(float)
    
    ### plotting

    ax[idx] = plt.subplot(gs[gs1,gs2])
    ax[idx].text(-0.07,1.03,label,fontsize=12,transform=ax[idx].transAxes) 
    
    ax[idx].plot(dmg_cd_nopro['rp'],dmg_cd_nopro['dmg'], color = c1, linewidth=2, label='RP-based without protection')
    ax[idx].plot(rp_bootstrap,max_nopro_mean,color=c2,linewidth=2,label='Event-based without protection')
    ax[idx].fill_between(rp_bootstrap,max_nopro_ci_lower,max_nopro_ci_upper,facecolor=c2,alpha=0.3)
    
    ax[idx].plot(dmg_cd_pro['rp'],dmg_cd_pro['dmg'], color = c3, alpha=0.75, linewidth=2, linestyle = "--",label='RP-based with protection')
    ax[idx].plot(rp_bootstrap,max_pro_mean,color=c4,linewidth=2,linestyle = "--",label='Event-based with protection')
    ax[idx].fill_between(rp_bootstrap,max_pro_ci_lower,max_pro_ci_upper,facecolor=c4,alpha=0.3)
    
    # ax[idx].plot(dmg_cd_nopro['rp'],dmg_cd_nopro['dmg']/scale, color = c1, linewidth=2, label='RP-based without protection',path_effects=[pe.Stroke(linewidth=5, foreground='w'), pe.Normal()])
    # ax[idx].plot(rp_bootstrap,max_nopro_mean/scale,color=c2,linewidth=2,label='Event-based without protection')
    # ax[idx].fill_between(rp_bootstrap,max_nopro_ci_lower/scale,max_nopro_ci_upper/scale,facecolor=c2,alpha=0.3)
    
    # ax[idx].plot(dmg_cd_pro['rp'],dmg_cd_pro['dmg']/scale, color = c3, alpha=0.75, linewidth=2, linestyle = "--",label='RP-based with protection',path_effects=[pe.Stroke(linewidth=5, foreground='w'), pe.Normal()])
    # ax[idx].plot(rp_bootstrap,max_pro_mean/scale,color=c4,linewidth=2,linestyle = "--",label='Event-based with protection')
    # ax[idx].fill_between(rp_bootstrap,max_pro_ci_lower/scale,max_pro_ci_upper/scale,facecolor=c4,alpha=0.3)

    
    ax[idx].set_xlim([1, 1000])
    ax[idx].set_ylim(ymin=0)
    ax[idx].set_xscale('log')
    ax[idx].grid(color='gray',linestyle = 'solid', alpha = 0.2)
    ax[idx].set_xticks([1, 2, 5, 10, 25, 50, 100, 250, 500, 1000])
    # ax[idx].set_xticklabels([])
    ax[idx].set_xticklabels([1, 2, 5, 10, 25, 50, 100, 250, 500, 1000])
    
    ax[idx].set_xlabel('Return period (years)',fontsize=12)
    ax[idx].set_ylabel('Flood economic damage (USD2005)',fontsize=12)
    
    if idx==0:
        ax[idx].legend(loc ='upper left',handlelength=2.6)
    
    # if gs2==0:
    #     ax[idx].set_ylabel('Economic damage (USD2005)')
    
    # if gs1==1:
        
    #     ax[idx].set_xticklabels([1, 2, 5, 10, 25, 50, 100, 250, 500, 1000])
        # ax[idx].set_xlabel('Return Period (year)')
    
    # ax[idx].set_title(label+' '+country,fontsize=14)
    ax[idx].set_title(country,fontsize=12)
    
# fig.text(0.5, 0.05, 'Return period (years)',fontsize=12, ha='center')
# fig.text(0.08, 0.5, 'Flood economic damage (USD2005)',fontsize=12, va='center', rotation='vertical')
outputdir = r'C:\Users\hli490\OneDrive - Vrije Universiteit Amsterdam\Desktop\Paper2\Figures\Risk_curve'
plt.savefig(os.path.join(outputdir,'impact_curve_supplement.png'),format='png',dpi=400,bbox_inches='tight')
# plt.savefig(os.path.join(outputdir,'impact_curve_seleted_countries_new.png'),format='png',dpi=225)
# plt.savefig(os.path.join(outputdir,'impact_curve_{:s}'.format(country)+'_{:03d}.png'.format(iso)),format='png',dpi=600,bbox_inches='tight')

t = max_nopro_mean[np.array([5,10,25,50,100,250,500,1000])]


ead_cd_nopro = np.trapz(dmg_cd_nopro['dmg'][::-1],1/dmg_cd_nopro['rp'][::-1])
ead_cd_nopro = np.trapz(dmg_cd_nopro['dmg'][2:10][::-1],1/dmg_cd_nopro['rp'][2:10][::-1])
28066960310.43103/40020373723.62424


ead_cd_pro = np.trapz(dmg_cd_pro['dmg'][::-1],1/dmg_cd_pro['rp'][::-1])
7063215888.921698
dmg_cd_pro = dmg_cd_pro.reset_index(drop=True)
np.trapz(dmg_cd_pro['dmg'][15:21][::-1],1/dmg_cd_pro['rp'][15:21][::-1])
6868524613.906331/7063215888.921698



# ax.plot(T,ci_lower,color=c3,linewidth=1,linestyle='dotted',label='Modelled dependence (5th-95th)')
# ax.plot(T,ci_upper,color=c3,linewidth=1,linestyle='dotted')

# for i in range(n):
    
#     # Event damage
#     event_dmg = np.array(event_dmg_all['dmg'].sample(n=int(len(event_dmg_all)/(10000/yr)),replace=True))
#     event_rp = get_return_periods(event_dmg, yr=1000)
#     event = pd.DataFrame(data=dict(dmg=event_dmg,rp=event_rp))
#     event = event.sort_values(by='rp')
#     event_bootstrap.iloc[i,:] = np.interp(rp_bootstrap,event['rp'],event['dmg'])

#     # Annual total damages
#     total_dmg = np.array(annual_total_all['dmg'].sample(n=int(len(annual_total_all)/(10000/yr)),replace=True))
#     total_rp = get_return_periods(total_dmg, yr=1000)
#     total = pd.DataFrame(data=dict(dmg=total_dmg,rp=total_rp))
#     total = total.sort_values(by='rp')
#     total_bootstrap.iloc[i,:] = np.interp(rp_bootstrap,total['rp'],total['dmg'])
    
#     # Annual max damages
#     max_dmg = np.array(annual_max_all['dmg'].sample(n=int(len(annual_max_all)/(10000/yr)),replace=True))
#     max_rp = get_return_periods(max_dmg, yr=1000)
#     max_annual = pd.DataFrame(data=dict(dmg=max_dmg,rp=max_rp))
#     max_annual = max_annual.sort_values(by='rp')
#     max_bootstrap.iloc[i,:] = np.interp(rp_bootstrap,max_annual['rp'],max_annual['dmg'])

# # Medians and confidence intervals
# event_median = event_bootstrap.median().astype(float)
# event_ci_lower = event_bootstrap.quantile(0.05, numeric_only=False).astype(float)
# event_ci_upper = event_bootstrap.quantile(0.95, numeric_only=False).astype(float)

# total_median = total_bootstrap.median().astype(float)
# total_ci_lower = total_bootstrap.quantile(0.05, numeric_only=False).astype(float)
# total_ci_upper = total_bootstrap.quantile(0.95, numeric_only=False).astype(float)

# max_median = max_bootstrap.median().astype(float)
# max_ci_lower = max_bootstrap.quantile(0.05, numeric_only=False).astype(float)
# max_ci_upper = max_bootstrap.quantile(0.95, numeric_only=False).astype(float)
