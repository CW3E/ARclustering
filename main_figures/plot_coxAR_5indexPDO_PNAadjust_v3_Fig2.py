#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 22:00:05 2024

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



col=(0.85,0.85,0.85)
linew=0.2
fs=14
fsc=15
flalo=11
s=0.3

#ppt
cmap=cm.seismic
cmap2=cm.seismic #,
cslev=[-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8]


mjolead_day=0

obj=nc.Dataset('globalARcatalog_MERRA2_1980-2021_v3.0all.nc')

lat=obj.variables['lat'][220:300]  #30-50N
lon=obj.variables['lon'][200:415]  #230-245


arcoxRMM2=np.loadtxt('PNAadjustMonMJO'+str(mjolead_day)+'/ARregion_coeffRMM21982-2021indp_nov2mar_lead'+str(mjolead_day)+'day.txt')
pvalRMM2=np.loadtxt('PNAadjustMonMJO'+str(mjolead_day)+'/ARregion_pvalRMM21982-2021indp_nov2mar_lead'+str(mjolead_day)+'day.txt')
arcoxRMM2[pvalRMM2==-9999]=np.nan


arcoxRMM1=np.loadtxt('PNAadjustMonMJO'+str(mjolead_day)+'/ARregion_coeffRMM11982-2021indp_nov2mar_lead'+str(mjolead_day)+'day.txt')
pvalRMM1=np.loadtxt('PNAadjustMonMJO'+str(mjolead_day)+'/ARregion_pvalRMM11982-2021indp_nov2mar_lead'+str(mjolead_day)+'day.txt')
arcoxRMM1[pvalRMM1==-9999]=np.nan

arcoxPNA=np.loadtxt('PNAadjustMonMJO'+str(mjolead_day)+'/ARregion_coeffPNAadjust1982-2021indp_nov2mar.txt')
pvalPNA=np.loadtxt('PNAadjustMonMJO'+str(mjolead_day)+'/ARregion_pvalPNAadjust1982-2021indp_nov2mar.txt')
arcoxPNA[pvalPNA==-9999]=np.nan

arcoxAO=np.loadtxt('PNAadjustMonMJO'+str(mjolead_day)+'/ARregion_coeffAO1982-2021indp_nov2mar.txt')
pvalAO=np.loadtxt('PNAadjustMonMJO'+str(mjolead_day)+'/ARregion_pvalAO1982-2021indp_nov2mar.txt')
arcoxAO[pvalAO==-9999]=np.nan


arcoxENSO=np.loadtxt('PNAadjustMonMJO'+str(mjolead_day)+'/ARregion_coeffENSO1982-2021indp_nov2mar.txt')
pvalENSO=np.loadtxt('PNAadjustMonMJO'+str(mjolead_day)+'/ARregion_pvalENSO1982-2021indp_nov2mar.txt')
arcoxENSO[pvalENSO==-9999]=np.nan

arcoxQBO=np.loadtxt('PNAadjustMonMJO'+str(mjolead_day)+'/ARregion_coeffQBO1982-2021indp_nov2mar.txt')
pvalQBO=np.loadtxt('PNAadjustMonMJO'+str(mjolead_day)+'/ARregion_pvalQBO1982-2021indp_nov2mar.txt')
arcoxQBO[pvalQBO==-9999]=np.nan


arcoxPDO=np.loadtxt('ARregion_coeffPDO1940-2018indp_nov2mar.txt')
pvalPDO=np.loadtxt('ARregion_pvalPDO1940-2018indp_nov2mar.txt')
arcoxPDO[pvalPDO==-9999]=np.nan




#-------------------
fig=plt.figure(1)#
ax = fig.add_axes([0.1,0.1,0.8,0.8])#  125 59
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
m.drawstates(linewidth=0.5)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 

cs=m.contourf(x1,y1,arcoxPNA,cslev,cmap=cm.seismic,extend='both')#!!!!!!
plt.title('PNA_noENSOrelated(NDJFM)',fontsize=fs)

parallels = np.arange(20.,60,10.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[1,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(120.,260.,20.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,1],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)
#cbar1=fig.colorbar(cs,cax=plt.axes([0.92, 0.3, 0.01, 0.4])) #51
#cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
#cbar1.set_ticks(cslev)
#cbar1.set_ticklabels((cslev))



ax = fig.add_axes([0.1,2.68,0.8,0.8])#  125 59
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
m.drawstates(linewidth=0.5)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 

cs=m.contourf(x1,y1,arcoxAO,cslev,cmap=cm.seismic,extend='both')#!!!!!!
plt.title('AO(NDJFM)',fontsize=fs)

parallels = np.arange(20.,60,10.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[1,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(120.,260.,20.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)
#cbar1=fig.colorbar(cs,cax=plt.axes([1.92, 0.3, 0.01, 0.4])) #51
#cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
#cbar1.set_ticks(cslev)
#cbar1.set_ticklabels((cslev))




ax = fig.add_axes([0.1,1.39,0.8,0.8])#  125 59
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
m.drawstates(linewidth=0.5)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 

cs=m.contourf(x1,y1,arcoxRMM1,cslev,cmap=cm.seismic,extend='both')#!!!!!!
plt.title('RMM1(NDJFM)',fontsize=fs)

parallels = np.arange(20.,60,10.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[1,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(120.,260.,20.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)
#cbar1=fig.colorbar(cs,cax=plt.axes([0.92, 0.8, 0.01, 0.4])) #51
#cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
#cbar1.set_ticks(cslev)
#cbar1.set_ticklabels((cslev))



ax = fig.add_axes([0.1,0.96,0.8,0.8])#  125 59
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
m.drawstates(linewidth=0.5)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 

cs=m.contourf(x1,y1,arcoxRMM2,cslev,cmap=cm.seismic,extend='both')#!!!!!!
plt.title('RMM2(NDJFM)',fontsize=fs)

parallels = np.arange(20.,60,10.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[1,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(120.,260.,20.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)
#cbar1=fig.colorbar(cs,cax=plt.axes([1.92, 0.8, 0.01, 0.4])) #51
#cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
#cbar1.set_ticks(cslev)
#cbar1.set_ticklabels((cslev))




ax = fig.add_axes([0.1,1.82,0.8,0.8])#  125 59
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
m.drawstates(linewidth=0.5)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 

cs=m.contourf(x1,y1,arcoxENSO,cslev,cmap=cm.seismic,extend='both')#!!!!!!
plt.title('ENSO(NDJFM)',fontsize=fs)

parallels = np.arange(20.,60,10.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[1,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(120.,260.,20.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)
#cbar1=fig.colorbar(cs,cax=plt.axes([0.92, 1.3, 0.01, 0.4])) #51
#cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
#cbar1.set_ticks(cslev)
#cbar1.set_ticklabels((cslev))




ax = fig.add_axes([0.1,2.25,0.8,0.8])#  125 59
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
m.drawstates(linewidth=0.5)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 

cs=m.contourf(x1,y1,arcoxQBO,cslev,cmap=cm.seismic,extend='both')#!!!!!!
plt.title('QBO(NDJFM)',fontsize=fs)

parallels = np.arange(20.,60,10.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[1,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(120.,260.,20.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)



#=======
obj=nc.Dataset('era5_ar_gwv3_1940-2022.nc')

latpdo=obj.variables['lat'][120:281]  #60-20N
lonpdo=obj.variables['lon'][480:1041]  #120-260

ax = fig.add_axes([0.1,0.53,0.8,0.8])#  125 59
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
m.drawstates(linewidth=0.5)

x,y=np.meshgrid(lonpdo,latpdo)
x1,y1=m(x,y) 

cs=m.contourf(x1,y1,arcoxPDO,cslev,cmap=cm.seismic,extend='both')#!!!!!!
plt.title('PDO(NDJFM)',fontsize=fs)

parallels = np.arange(20.,60,10.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[1,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(120.,260.,20.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)

'''
cbar1=fig.colorbar(cs,cax=plt.axes([1.76, 0.1, 0.02, 1.25])) #51
cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
cbar1.set_ticks(cslev)
cbar1.set_ticklabels((cslev))
'''
del arcoxRMM2
del arcoxRMM1
del arcoxPNA
del arcoxAO
del arcoxENSO
del arcoxQBO
del arcoxPDO
#==========ndj=====


obj=nc.Dataset('globalARcatalog_MERRA2_1980-2021_v3.0all.nc')

lat=obj.variables['lat'][220:300]  #30-50N
lon=obj.variables['lon'][200:415]  #230-245


arcoxRMM2=np.loadtxt('ARregion_coeffRMM21982-2021indp_NDJ.txt')
pvalRMM2=np.loadtxt('ARregion_pvalRMM21982-2021indp_NDJ.txt')
arcoxRMM2[pvalRMM2==-9999]=np.nan


arcoxRMM1=np.loadtxt('ARregion_coeffRMM11982-2021indp_NDJ.txt')
pvalRMM1=np.loadtxt('ARregion_pvalRMM11982-2021indp_NDJ.txt')
arcoxRMM1[pvalRMM1==-9999]=np.nan

arcoxPNA=np.loadtxt('ARregion_coeffPNA1982-2021indp_NDJ.txt')
pvalPNA=np.loadtxt('ARregion_pvalPNA1982-2021indp_NDJ.txt')
arcoxPNA[pvalPNA==-9999]=np.nan

arcoxAO=np.loadtxt('ARregion_coeffAO1982-2021indp_NDJ.txt')
pvalAO=np.loadtxt('ARregion_pvalAO1982-2021indp_NDJ.txt')
arcoxAO[pvalAO==-9999]=np.nan


arcoxENSO=np.loadtxt('ARregion_coeffENSO1982-2021indp_NDJ.txt')
pvalENSO=np.loadtxt('ARregion_pvalENSO1982-2021indp_NDJ.txt')
arcoxENSO[pvalENSO==-9999]=np.nan

arcoxQBO=np.loadtxt('ARregion_coeffQBO1982-2021indp_NDJ.txt')
pvalQBO=np.loadtxt('ARregion_pvalQBO1982-2021indp_NDJ.txt')
arcoxQBO[pvalQBO==-9999]=np.nan

arcoxPDO=np.loadtxt('ARregion_coeffPDO1940-2018indp_NDJ.txt')
pvalPDO=np.loadtxt('ARregion_pvalPDO1940-2018indp_NDJ.txt')
arcoxPDO[pvalPDO==-9999]=np.nan






#-------------------
fig=plt.figure(1)#
ax = fig.add_axes([0.93,0.1,0.8,0.8])#  125 59
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
m.drawstates(linewidth=0.5)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 

cs=m.contourf(x1,y1,arcoxPNA,cslev,cmap=cm.seismic,extend='both')#!!!!!!
plt.title('PNA_noENSOrelated(NDJ)',fontsize=fs)

parallels = np.arange(20.,60,10.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(120.,260.,20.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,1],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)
#cbar1=fig.colorbar(cs,cax=plt.axes([0.92, 0.3, 0.01, 0.4])) #51
#cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
#cbar1.set_ticks(cslev)
#cbar1.set_ticklabels((cslev))



ax = fig.add_axes([0.93,2.68,0.8,0.8])#  125 59
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
m.drawstates(linewidth=0.5)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 

cs=m.contourf(x1,y1,arcoxAO,cslev,cmap=cm.seismic,extend='both')#!!!!!!
plt.title('AO(NDJ)',fontsize=fs)

parallels = np.arange(20.,60,10.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(120.,260.,20.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)
#cbar1=fig.colorbar(cs,cax=plt.axes([1.92, 0.3, 0.01, 0.4])) #51
#cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
#cbar1.set_ticks(cslev)
#cbar1.set_ticklabels((cslev))




ax = fig.add_axes([0.93,1.39,0.8,0.8])#  125 59
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
m.drawstates(linewidth=0.5)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 

cs=m.contourf(x1,y1,arcoxRMM1,cslev,cmap=cm.seismic,extend='both')#!!!!!!
plt.title('RMM1(NDJ)',fontsize=fs)

parallels = np.arange(20.,60,10.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(120.,260.,20.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)
#cbar1=fig.colorbar(cs,cax=plt.axes([0.92, 0.8, 0.01, 0.4])) #51
#cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
#cbar1.set_ticks(cslev)
#cbar1.set_ticklabels((cslev))



ax = fig.add_axes([0.93,0.96,0.8,0.8])#  125 59
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
m.drawstates(linewidth=0.5)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 

cs=m.contourf(x1,y1,arcoxRMM2,cslev,cmap=cm.seismic,extend='both')#!!!!!!
plt.title('RMM2(NDJ)',fontsize=fs)

parallels = np.arange(20.,60,10.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(120.,260.,20.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)
#cbar1=fig.colorbar(cs,cax=plt.axes([1.92, 0.8, 0.01, 0.4])) #51
#cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
#cbar1.set_ticks(cslev)
#cbar1.set_ticklabels((cslev))




ax = fig.add_axes([0.93,1.82,0.8,0.8])#  125 59
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
m.drawstates(linewidth=0.5)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 

cs=m.contourf(x1,y1,arcoxENSO,cslev,cmap=cm.seismic,extend='both')#!!!!!!
plt.title('ENSO(NDJ)',fontsize=fs)

parallels = np.arange(20.,60,10.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(120.,260.,20.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)
#cbar1=fig.colorbar(cs,cax=plt.axes([0.92, 1.3, 0.01, 0.4])) #51
#cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
#cbar1.set_ticks(cslev)
#cbar1.set_ticklabels((cslev))




ax = fig.add_axes([0.93,2.25,0.8,0.8])#  125 59
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
m.drawstates(linewidth=0.5)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 

cs=m.contourf(x1,y1,arcoxQBO,cslev,cmap=cm.seismic,extend='both')#!!!!!!
plt.title('QBO(NDJ)',fontsize=fs)

parallels = np.arange(20.,60,10.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(120.,260.,20.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)



#=======
obj=nc.Dataset('/Users/zhiqiyang/Downloads/era5_ar_gwv3_1940-2022.nc')

latpdo=obj.variables['lat'][120:281]  #60-20N
lonpdo=obj.variables['lon'][480:1041]  #120-260

ax = fig.add_axes([0.93,0.53,0.8,0.8])#  125 59
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
m.drawstates(linewidth=0.5)

x,y=np.meshgrid(lonpdo,latpdo)
x1,y1=m(x,y) 

cs=m.contourf(x1,y1,arcoxPDO,cslev,cmap=cm.seismic,extend='both')#!!!!!!
plt.title('PDO(NDJ)',fontsize=fs)

parallels = np.arange(20.,60,10.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(120.,260.,20.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)

del arcoxRMM2
del arcoxRMM1
del arcoxPNA
del arcoxAO
del arcoxENSO
del arcoxQBO
del arcoxPDO
#============

#======jfm======


obj=nc.Dataset('globalARcatalog_MERRA2_1980-2021_v3.0all.nc')

lat=obj.variables['lat'][220:300]  #30-50N
lon=obj.variables['lon'][200:415]  #230-245


arcoxRMM2=np.loadtxt('ARregion_coeffRMM21982-2021indp_JFM.txt')
pvalRMM2=np.loadtxt('ARregion_pvalRMM21982-2021indp_JFM.txt')
arcoxRMM2[pvalRMM2==-9999]=np.nan


arcoxRMM1=np.loadtxt('ARregion_coeffRMM11982-2021indp_JFM.txt')
pvalRMM1=np.loadtxt('ARregion_pvalRMM11982-2021indp_JFM.txt')
arcoxRMM1[pvalRMM1==-9999]=np.nan

arcoxPNA=np.loadtxt('ARregion_coeffPNA1982-2021indp_JFM.txt')
pvalPNA=np.loadtxt('ARregion_pvalPNA1982-2021indp_JFM.txt')
arcoxPNA[pvalPNA==-9999]=np.nan

arcoxAO=np.loadtxt('ARregion_coeffAO1982-2021indp_JFM.txt')
pvalAO=np.loadtxt('ARregion_pvalAO1982-2021indp_JFM.txt')
arcoxAO[pvalAO==-9999]=np.nan


arcoxENSO=np.loadtxt('ARregion_coeffENSO1982-2021indp_JFM.txt')
pvalENSO=np.loadtxt('ARregion_pvalENSO1982-2021indp_JFM.txt')
arcoxENSO[pvalENSO==-9999]=np.nan

arcoxQBO=np.loadtxt('ARregion_coeffQBO1982-2021indp_JFM.txt')
pvalQBO=np.loadtxt('ARregion_pvalQBO1982-2021indp_JFM.txt')
arcoxQBO[pvalQBO==-9999]=np.nan

arcoxPDO=np.loadtxt('ARregion_coeffPDO1940-2018indp_JFM.txt')
pvalPDO=np.loadtxt('ARregion_pvalPDO1940-2018indp_JFM.txt')
arcoxPDO[pvalPDO==-9999]=np.nan






#-------------------
fig=plt.figure(1)#
ax = fig.add_axes([1.76,0.1,0.8,0.8])#  125 59
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
m.drawstates(linewidth=0.5)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 

cs=m.contourf(x1,y1,arcoxPNA,cslev,cmap=cm.seismic,extend='both')#!!!!!!
plt.title('PNA_noENSOrelated(JFM)',fontsize=fs)

parallels = np.arange(20.,60,10.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(120.,260.,20.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,1],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)
#cbar1=fig.colorbar(cs,cax=plt.axes([0.92, 0.3, 0.01, 0.4])) #51
#cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
#cbar1.set_ticks(cslev)
#cbar1.set_ticklabels((cslev))



ax = fig.add_axes([1.76,2.68,0.8,0.8])#  125 59
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
m.drawstates(linewidth=0.5)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 

cs=m.contourf(x1,y1,arcoxAO,cslev,cmap=cm.seismic,extend='both')#!!!!!!
plt.title('AO(JFM)',fontsize=fs)

parallels = np.arange(20.,60,10.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(120.,260.,20.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)
#cbar1=fig.colorbar(cs,cax=plt.axes([1.92, 0.3, 0.01, 0.4])) #51
#cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
#cbar1.set_ticks(cslev)
#cbar1.set_ticklabels((cslev))




ax = fig.add_axes([1.76,1.39,0.8,0.8])#  125 59
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
m.drawstates(linewidth=0.5)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 

cs=m.contourf(x1,y1,arcoxRMM1,cslev,cmap=cm.seismic,extend='both')#!!!!!!
plt.title('RMM1(JFM)',fontsize=fs)

parallels = np.arange(20.,60,10.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(120.,260.,20.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)
#cbar1=fig.colorbar(cs,cax=plt.axes([0.92, 0.8, 0.01, 0.4])) #51
#cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
#cbar1.set_ticks(cslev)
#cbar1.set_ticklabels((cslev))



ax = fig.add_axes([1.76,0.96,0.8,0.8])#  125 59
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
m.drawstates(linewidth=0.5)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 

cs=m.contourf(x1,y1,arcoxRMM2,cslev,cmap=cm.seismic,extend='both')#!!!!!!
plt.title('RMM2(JFM)',fontsize=fs)

parallels = np.arange(20.,60,10.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(120.,260.,20.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)
#cbar1=fig.colorbar(cs,cax=plt.axes([1.92, 0.8, 0.01, 0.4])) #51
#cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
#cbar1.set_ticks(cslev)
#cbar1.set_ticklabels((cslev))




ax = fig.add_axes([1.76,1.82,0.8,0.8])#  125 59
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
m.drawstates(linewidth=0.5)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 

cs=m.contourf(x1,y1,arcoxENSO,cslev,cmap=cm.seismic,extend='both')#!!!!!!
plt.title('ENSO(JFM)',fontsize=fs)

parallels = np.arange(20.,60,10.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(120.,260.,20.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)
#cbar1=fig.colorbar(cs,cax=plt.axes([0.92, 1.3, 0.01, 0.4])) #51
#cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
#cbar1.set_ticks(cslev)
#cbar1.set_ticklabels((cslev))




ax = fig.add_axes([1.76,2.25,0.8,0.8])#  125 59
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
m.drawstates(linewidth=0.5)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 

cs=m.contourf(x1,y1,arcoxQBO,cslev,cmap=cm.seismic,extend='both')#!!!!!!
plt.title('QBO(JFM)',fontsize=fs)

parallels = np.arange(20.,60,10.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(120.,260.,20.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)


#=======
obj=nc.Dataset('era5_ar_gwv3_1940-2022.nc')

latpdo=obj.variables['lat'][120:281]  #60-20N
lonpdo=obj.variables['lon'][480:1041]  #120-260

ax = fig.add_axes([1.76,0.53,0.8,0.8])#  125 59
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
m.drawstates(linewidth=0.5)

x,y=np.meshgrid(lonpdo,latpdo)
x1,y1=m(x,y) 

cs=m.contourf(x1,y1,arcoxPDO,cslev,cmap=cm.seismic,extend='both')#!!!!!!
plt.title('PDO(JFM)',fontsize=fs)

parallels = np.arange(20.,60,10.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(120.,260.,20.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)






cbar1=fig.colorbar(cs,cax=plt.axes([2.59, 0.3, 0.02, 2.95])) #51
cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
cbar1.set_ticks(cslev)
cbar1.set_ticklabels((cslev))



savefig('Cox_regression_coefficient1982-2021indp_5indPDO_PNAadjustMon_v3_Fig2-response.pdf',bbox_inches = 'tight')










