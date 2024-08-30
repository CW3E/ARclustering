#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  2 00:06:51 2023

@author: zhiqiyang
"""



import matplotlib.pyplot as plt
import glob
from scipy import stats
import pandas as pd
import os 
import sys
from datetime import datetime, timedelta
from matplotlib.pyplot import savefig
import netCDF4 as nc
import numpy as np
from mpl_toolkits.basemap import Basemap
import numpy.ma as ma
import matplotlib.cm as cm
from numpy import arange
import matplotlib.patches as patches


flalo=8
s=4
col=(0.85,0.85,0.85)
linew=0.2
fsc=10
fs=12
fscoast=0.9
cmap=cm.viridis
cmap2=cm.PuOr_r
 
obj=nc.Dataset('/Users/zhiqiyang/Downloads/globalARcatalog_MERRA2_1980-2021_v3.0all.nc')
lat=obj.variables['lat'][240:280]  #30-50N
lon=obj.variables['lon'][368:400]  #30-50N

year=obj.variables['year'][0,:,0,:]
mon=obj.variables['month'][0,:,0,:]
day=obj.variables['day'][0,:,0,:]

#[55,56,57,58,59,60,61,62,63,64,65,66,67,68]
cslev=[56,57,58,59,60,61,62,63,64,65,66,67,68,69]
 
ivt_clus_ndj=np.load('ivtdir_uniAR/kivtdir_uniARclust1wk_daily1982-2021NDJ_lat240280lon368400.npy')
ivt_clus_jfm=np.load('ivtdir_uniAR/kivtdir_uniARclust1wk_daily1982-2021JFM_lat240280lon368400.npy')
ivt_clus_ndjfm=np.load('ivtdir_uniAR/kivtdir_uniARclust1wk_daily1982-2021NDJFM_lat240280lon368400.npy')

ivt_all_ndj=np.load('ivtdir_uniAR/kivtdir_uniAR_daily1982-2021NDJ_lat240280lon368400.npy')
ivt_all_jfm=np.load('ivtdir_uniAR/kivtdir_uniAR_daily1982-2021JFM_lat240280lon368400.npy')
ivt_all_ndjfm=np.load('ivtdir_uniAR/kivtdir_uniAR_daily1982-2021NDJFM_lat240280lon368400.npy')

#=============clust============
#a1=240
#a2=280

#b1=368
#b2=400

ivtdir_dailyave=np.copy(ivt_clus_ndjfm)
ivtdir_dailyave=np.nanmean(ivtdir_dailyave,0)

'''
fig=plt.figure(1)#
ax = fig.add_axes([-0.72,0.9,0.8,0.8])#  125 59     0.1,0.1
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=30,urcrnrlat=47,llcrnrlon=230,urcrnrlon=249,resolution='l',ax=ax)
m.drawcoastlines(linewidth=fscoast)
m.drawcountries(linewidth=fscoast)
m.drawstates(linewidth=fscoast)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 
#
cs=m.contourf(x1,y1,ivtdir_dailyave,cslev,cmap=cmap,extend='both')
plt.title('Clustered Unique ARs life-cycle_dir NDJFM',fontsize=fs)

parallels = np.arange(20.,60,5.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[1,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(235.,250.,5.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)
#cbar1=fig.colorbar(cs,cax=plt.axes([0.72, 0.9, 0.03, 0.8])) #51
#cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
#cbar1.set_ticks(cslev)
#cbar1.set_ticklabels((cslev))
'''



#----------------------
ivtdir_dailyave=np.copy(ivt_clus_ndj)
ivtdir_dailyave=np.nanmean(ivtdir_dailyave,0)


fig=plt.figure(1)#
ax = fig.add_axes([0,0.9,0.8,0.8])#  125 59     0.1,0.1
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=30,urcrnrlat=47,llcrnrlon=230,urcrnrlon=249,resolution='l',ax=ax)
m.drawcoastlines(linewidth=fscoast)
m.drawcountries(linewidth=fscoast)
m.drawstates(linewidth=fscoast)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 
#
cs=m.contourf(x1,y1,ivtdir_dailyave,cslev,cmap=cmap,extend='both')
plt.title('Clustered Unique ARs life-cycle_dir NDJ',fontsize=fs)

parallels = np.arange(20.,60,5.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[1,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(235.,250.,5.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)
#cbar1=fig.colorbar(cs,cax=plt.axes([0.72, 0.9, 0.03, 0.8])) #51
#cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
#cbar1.set_ticks(cslev)
#cbar1.set_ticklabels((cslev))

#=====

ivtdir_daily=np.copy(ivt_clus_jfm)
ivtdir_dailyave=np.nanmean(ivtdir_daily,0)

fig=plt.figure(1)#
ax = fig.add_axes([0.72,0.9,0.8,0.8])#  125 59     0.1,0.1
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=30,urcrnrlat=47,llcrnrlon=230,urcrnrlon=249,resolution='l',ax=ax)
m.drawcoastlines(linewidth=fscoast)
m.drawcountries(linewidth=fscoast)
m.drawstates(linewidth=fscoast)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 
#
cs=m.contourf(x1,y1,ivtdir_dailyave,cslev,cmap=cmap,extend='both')
plt.title('Clustered Unique ARs life-cycle_dir JFM ',fontsize=fs)


parallels = np.arange(20.,60,5.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(235.,250.,5.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)
#cbar1=fig.colorbar(cs,cax=plt.axes([1.47, 0.9, 0.03, 0.8])) #51
#cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
#cbar1.set_ticks(cslev)
#cbar1.set_ticklabels((cslev))

#=========------all-----
ivtdir_dailyave=np.copy(ivt_all_ndjfm)
ivtdir_dailyave=np.nanmean(ivtdir_dailyave,0)

'''
fig=plt.figure(1)#
ax = fig.add_axes([-0.72,0,0.8,0.8])#  125 59     0.1,0.1
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=30,urcrnrlat=47,llcrnrlon=230,urcrnrlon=249,resolution='l',ax=ax)
m.drawcoastlines(linewidth=fscoast)
m.drawcountries(linewidth=fscoast)
m.drawstates(linewidth=fscoast)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 
#
cs=m.contourf(x1,y1,ivtdir_dailyave,cslev,cmap=cmap,extend='both')
plt.title('Unique ARs life-cycle_dir NDJFM',fontsize=fs)

parallels = np.arange(20.,60,5.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[1,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(235.,250.,5.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)
cbar1=fig.colorbar(cs,cax=plt.axes([0, 0, 0.03, 1.7])) #51 1.7
cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
cbar1.set_ticks(cslev)
cbar1.set_ticklabels((cslev))
'''






#-------------
ivtdir_dailyave=np.copy(ivt_all_ndj)
ivtdir_dailyave=np.nanmean(ivtdir_dailyave,0)


fig=plt.figure(1)#
ax = fig.add_axes([0,0,0.8,0.8])#  125 59     0.1,0.1
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=30,urcrnrlat=47,llcrnrlon=230,urcrnrlon=249,resolution='l',ax=ax)
m.drawcoastlines(linewidth=fscoast)
m.drawcountries(linewidth=fscoast)
m.drawstates(linewidth=fscoast)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 
#
cs=m.contourf(x1,y1,ivtdir_dailyave,cslev,cmap=cmap,extend='both')
plt.title('Unique ARs life-cycle_dir NDJ',fontsize=fs)

parallels = np.arange(20.,60,5.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[1,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(235.,250.,5.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)
cbar1=fig.colorbar(cs,cax=plt.axes([0.72, 0, 0.03, 1.7])) #51
cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
cbar1.set_ticks(cslev)
cbar1.set_ticklabels((cslev))

#=====

ivtdir_daily=np.copy(ivt_all_jfm)
ivtdir_dailyave=np.nanmean(ivtdir_daily,0)

fig=plt.figure(1)#
ax = fig.add_axes([0.72,0,0.8,0.8])#  125 59     0.1,0.1
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=30,urcrnrlat=47,llcrnrlon=230,urcrnrlon=249,resolution='l',ax=ax)
m.drawcoastlines(linewidth=fscoast)
m.drawcountries(linewidth=fscoast)
m.drawstates(linewidth=fscoast)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 
#
cs=m.contourf(x1,y1,ivtdir_dailyave,cslev,cmap=cmap,extend='both')
plt.title('Unique ARs life-cycle_dir JFM ',fontsize=fs)

parallels = np.arange(20.,60,5.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(235.,250.,5.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)
cbar1=fig.colorbar(cs,cax=plt.axes([1.44, 0, 0.03, 1.7])) #51
cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
cbar1.set_ticks(cslev)
cbar1.set_ticklabels((cslev))

#========================

#======diff===
cslev=[-8,-6,-4,-2,-1,-0.5,-0.1,0,0.1,0.5,1,2,4,6,8]


clus_avendj=np.nanmean(ivt_clus_ndj,0)
clus_avejfm=np.nanmean(ivt_clus_jfm,0)
clus_avendjfm=np.nanmean(ivt_clus_ndjfm,0)

all_avendj=np.nanmean(ivt_all_ndj,0)
all_avejfm=np.nanmean(ivt_all_jfm,0)
all_avendjfm=np.nanmean(ivt_all_ndjfm,0)


diff_clus_season=clus_avendj-clus_avejfm
diff_all_season=all_avendj-all_avejfm
diff_clus_all_ndj=clus_avendj-all_avendj
diff_clus_all_jfm=clus_avejfm-all_avejfm
diff_clus_all_ndjfm=clus_avendjfm-all_avendjfm


#------diff_clus_season---
fig=plt.figure(1)#
ax = fig.add_axes([1.44,0.9,0.8,0.8])#  125 59     0.1,0.1
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=30,urcrnrlat=47,llcrnrlon=230,urcrnrlon=249,resolution='l',ax=ax)
m.drawcoastlines(linewidth=fscoast)
m.drawcountries(linewidth=fscoast)
m.drawstates(linewidth=fscoast)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 
#
cs=m.contourf(x1,y1,diff_clus_season,cslev,cmap=cmap2,extend='both')
plt.title('Clustered Unique ARs NDJ-JFM',fontsize=fs)

for i in arange(0,np.shape(lat)[0]):
    for j in arange(0,np.shape(lon)[0]):
        
        lin1=ivt_clus_jfm[:,i,j]
        lin2=ivt_clus_ndj[:,i,j]
        
        lin1=lin1[~np.isnan(lin1)]
        lin2=lin2[~np.isnan(lin2)]
        
        t_statistic, sig = stats.ttest_ind(lin1, lin2)
        
        if (sig < 0.05 ):
            m.plot(x1[i,j],y1[i,j],'+',color='k',markersize=s+2,zorder=5)
            del sig

parallels = np.arange(20.,60,5.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(235.,250.,5.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)


#cbar1=fig.colorbar(cs,cax=plt.axes([2.22, 0.9, 0.03, 0.8])) #51
#cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
#cbar1.set_ticks(cslev)
#cbar1.set_ticklabels((cslev))
#=========
#------diff_all_season---
fig=plt.figure(1)#
ax = fig.add_axes([1.44,0,0.8,0.8])#  125 59     0.1,0.1
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=30,urcrnrlat=47,llcrnrlon=230,urcrnrlon=249,resolution='l',ax=ax)
m.drawcoastlines(linewidth=fscoast)
m.drawcountries(linewidth=fscoast)
m.drawstates(linewidth=fscoast)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 
#
cs=m.contourf(x1,y1,diff_all_season,cslev,cmap=cmap2,extend='both')
plt.title('Unique ARs NDJ-JFM',fontsize=fs)

for i in arange(0,np.shape(lat)[0]):
    for j in arange(0,np.shape(lon)[0]):
        
        lin1=ivt_all_jfm[:,i,j]
        lin2=ivt_all_ndj[:,i,j]
        
        lin1=lin1[~np.isnan(lin1)]
        lin2=lin2[~np.isnan(lin2)]
        
        t_statistic, sig = stats.ttest_ind(lin1, lin2)
        
        if (sig < 0.05 ):
            m.plot(x1[i,j],y1[i,j],'+',color='k',markersize=s,zorder=5)
            del sig

parallels = np.arange(20.,60,5.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(235.,250.,5.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,1],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)


cbar1=fig.colorbar(cs,cax=plt.axes([2.16, 0, 0.03, 1.7])) #51
cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
cbar1.set_ticks(cslev)
cbar1.set_ticklabels((cslev))

#---------diff_clus_all_ndjfm----------
'''
fig=plt.figure(1)#
ax = fig.add_axes([-0.72,-0.9,0.8,0.8])#  125 59     0.1,0.1
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=30,urcrnrlat=47,llcrnrlon=230,urcrnrlon=249,resolution='l',ax=ax)
m.drawcoastlines(linewidth=fscoast)
m.drawcountries(linewidth=fscoast)
m.drawstates(linewidth=fscoast)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 
#
cs=m.contourf(x1,y1,diff_clus_all_ndjfm,cslev,cmap=cmap2,extend='both')
plt.title('Clustered_Unique_ARs-Unique_ARs NDJFM',fontsize=fs)

for i in arange(0,np.shape(lat)[0]):
    for j in arange(0,np.shape(lon)[0]):
        
        lin1=ivt_all_ndjfm[:,i,j]
        lin2=ivt_clus_ndjfm[:,i,j]
        
        lin1=lin1[~np.isnan(lin1)]
        lin2=lin2[~np.isnan(lin2)]
        
        t_statistic, sig = stats.ttest_ind(lin1, lin2)
        
        if (sig < 0.05 ):
            m.plot(x1[i,j],y1[i,j],'+',color='k',markersize=s,zorder=5)
            del sig

parallels = np.arange(20.,60,5.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[1,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(235.,250.,5.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,1],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)


cbar1=fig.colorbar(cs,cax=plt.axes([0, -0.9, 0.03, 0.8])) #51
cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
cbar1.set_ticks(cslev)
cbar1.set_ticklabels((cslev))
'''



#---------diff_clus_all_ndj----------
fig=plt.figure(1)#
ax = fig.add_axes([0,-0.9,0.8,0.8])#  125 59     0.1,0.1
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=30,urcrnrlat=47,llcrnrlon=230,urcrnrlon=249,resolution='l',ax=ax)
m.drawcoastlines(linewidth=fscoast)
m.drawcountries(linewidth=fscoast)
m.drawstates(linewidth=fscoast)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 
#
cs=m.contourf(x1,y1,diff_clus_all_ndj,cslev,cmap=cmap2,extend='both')
plt.title('Clustered_Unique_ARs-Unique_ARs AR NDJ',fontsize=fs)

for i in arange(0,np.shape(lat)[0]):
    for j in arange(0,np.shape(lon)[0]):
        
        lin1=ivt_all_ndj[:,i,j]
        lin2=ivt_clus_ndj[:,i,j]
        
        lin1=lin1[~np.isnan(lin1)]
        lin2=lin2[~np.isnan(lin2)]
        
        t_statistic, sig = stats.ttest_ind(lin1, lin2)
        
        if (sig < 0.05 ):
            m.plot(x1[i,j],y1[i,j],'+',color='k',markersize=s,zorder=5)
            del sig

parallels = np.arange(20.,60,5.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[1,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(235.,250.,5.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,1],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)


cbar1=fig.colorbar(cs,cax=plt.axes([0.72, -0.9, 0.03, 0.8])) #51
cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
cbar1.set_ticks(cslev)
cbar1.set_ticklabels((cslev))


#--------diff_clus_all_jfm-----------
fig=plt.figure(1)#
ax = fig.add_axes([0.72,-0.9,0.8,0.8])#  125 59     0.1,0.1
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=30,urcrnrlat=47,llcrnrlon=230,urcrnrlon=249,resolution='l',ax=ax)
m.drawcoastlines(linewidth=fscoast)
m.drawcountries(linewidth=fscoast)
m.drawstates(linewidth=fscoast)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 
#
cs=m.contourf(x1,y1,diff_clus_all_jfm,cslev,cmap=cmap2,extend='both')
plt.title('Clustered_Unique_ARs-Unique_ARs JFM',fontsize=fs)

for i in arange(0,np.shape(lat)[0]):
    for j in arange(0,np.shape(lon)[0]):
        
        lin1=ivt_all_jfm[:,i,j]
        lin2=ivt_clus_jfm[:,i,j]
        
        lin1=lin1[~np.isnan(lin1)]
        lin2=lin2[~np.isnan(lin2)]
        
        t_statistic, sig = stats.ttest_ind(lin1, lin2)
        
        if (sig < 0.05 ):
            m.plot(x1[i,j],y1[i,j],'+',color='k',markersize=s,zorder=5)
            del sig

parallels = np.arange(20.,60,5.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(235.,250.,5.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,1],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)


cbar1=fig.colorbar(cs,cax=plt.axes([1.44, -0.9, 0.03, 0.8])) 
cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
cbar1.set_ticks(cslev)
cbar1.set_ticklabels((cslev))

#=========
'''
fig=plt.figure(1)
ax = fig.add_axes([1.44,-0.9,0.8,0.8])#  125 59     0.1,0.1


ax.axis('off')
arrow_length = 0.1

start_point = np.array([0.1, 0.4])
end_point = np.array([arrow_length * np.cos(np.radians(60)), arrow_length * np.sin(np.radians(60))])

arrow = ax.arrow(start_point[0], start_point[1], end_point[0], end_point[1], head_width=0.01, head_length=0.02, fc='black', ec='black', linewidth=2)

ax.text(0.13, 0.18, "IVT orientation is defined as \nclockwise degrees from North", fontsize=12, ha='center')
ax.plot([0, 0.2],[0, 0.4], linestyle='--', color='w')
ax.plot([0.1, 0.1], [0.4, 0.8], linestyle='--', color='gray')

ax.text(0.11, 0.6, "60 degrees", fontsize=12, va='center')

'''
#=========
savefig('KivtdirClimatology_daily_uniAR_control_1982-2021_3fig-lifecycle-response.pdf',bbox_inches = 'tight')
   
   
plt.show()



