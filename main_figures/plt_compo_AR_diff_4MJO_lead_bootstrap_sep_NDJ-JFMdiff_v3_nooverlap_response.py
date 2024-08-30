#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 19 11:48:41 2023

@author: zhiqiyang
"""


import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
from scipy import stats
import pandas as pd
import os 
from matplotlib.pyplot import savefig
from numpy import arange,array,ones,ndarray
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap


h=[0.1,0.92,0.1,0.92,0.1,0.92,0.1,0.92] #[0.5*5+0.1,0.5*4+0.1,0.5*3+0.1,0.5*2+0.1,,]
z=[0.5*3+0.1,0.5*3+0.1,0.5*2+0.1,0.5*2+0.1,0.5+0.1,0.5+0.1,0.1,0.1]

col=(0.85,0.85,0.85)
linew=0.2
fs=12
fsc=13
flalo=10
s=2.5
cmap=cm.PuOr_r


for wind in np.arange(1,2): #week window
 

   
    
    if wind==1:
        cslev=[-0.5,-0.4,-0.3,-0.2,-0.1,-0.03,0,0.03,0.1,0.2,0.3,0.4,0.5]
    elif wind==2:
        cslev=[-1,-0.8,-0.6,-0.4,-0.2,-0.1,0,0.1,0.2,0.4,0.6,0.8,1]
    elif wind==3:
        cslev=[-1,-0.8,-0.6,-0.4,-0.2,-0.1,0,0.1,0.2,0.4,0.6,0.8,1]
    elif wind==4:
        cslev=[-1,-0.8,-0.6,-0.4,-0.2,-0.1,0,0.1,0.2,0.4,0.6,0.8,1]
    elif wind==5:
        cslev=[-1,-0.8,-0.6,-0.4,-0.2,-0.1,0,0.1,0.2,0.4,0.6,0.8,1]
    elif wind==6:
        cslev=[-1,-0.8,-0.6,-0.4,-0.2,-0.1,0,0.1,0.2,0.4,0.6,0.8,1]
    elif wind==7:
        cslev=[-1,-0.8,-0.6,-0.4,-0.2,-0.1,0,0.1,0.2,0.4,0.6,0.8,1]
    elif wind==8:
        cslev=[-1,-0.8,-0.6,-0.4,-0.2,-0.1,0,0.1,0.2,0.4,0.6,0.8,1]
    
    obj=nc.Dataset('/Users/zhiqiyang/Downloads/globalARcatalog_MERRA2_1980-2021_v3.0all.nc')
    
    lat=obj.variables['lat'][220:300]  #30-50N
    lon=obj.variables['lon'][200:415]  #230-245
    
    #======================
    month='NDJ'
    month2='NDJ'
        
    compMJO1ndj=np.zeros((80,215))
    compMJO2ndj=np.zeros((80,215))
    compMJO3ndj=np.zeros((80,215))
    compMJO4ndj=np.zeros((80,215))
    
    compMJO1ndj_all=np.zeros((80,215,512))
    compMJO2ndj_all=np.zeros((80,215,512))
    compMJO3ndj_all=np.zeros((80,215,512))
    compMJO4ndj_all=np.zeros((80,215,512))
  
    
    
    #index--winter
    mjo=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/mjo8phase_amp1982-2021_'+month+'_nowk1lead5.txt')
    mjoph=mjo[:,0]
    amp=mjo[:,1]
    
    m1=np.copy(mjoph)
    m2=np.copy(mjoph)
    m3=np.copy(mjoph)
    m4=np.copy(mjoph)
    
    
    
    m1[(m1 >= 2) & (m1 <= 7)] = np.nan
    m1[amp<1]= np.nan
    
    m2[(m2 != 2) & (m2 != 3)] = np.nan
    m2[amp<1]= np.nan
    
    m3[(m3 != 4) & (m3 != 5)] = np.nan
    m3[amp<1]= np.nan
    
    m4[(m4 != 6) & (m4 != 7)] = np.nan
    m4[amp<1]= np.nan
    
    
    arwinter=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
    climatndj=np.nanmean(arwinter,2)

    #--1982-2021 ndjfm
    for i in np.arange(220,300):
        
                
        for j in np.arange(0,215):
            
            e1=arwinter[i-220,j,1:513] #no week1 MJO
 
            even1=np.copy(e1)
            even2=np.copy(even1)
            even3=np.copy(even1)
            even4=np.copy(even1)
            
            even1[np.isnan(m1[0:np.shape(even1)[0]])]=np.nan
            even2[np.isnan(m2[0:np.shape(even1)[0]])]=np.nan
            even3[np.isnan(m3[0:np.shape(even1)[0]])]=np.nan
            even4[np.isnan(m4[0:np.shape(even1)[0]])]=np.nan
            
                        
            compMJO1ndj[i-220,j]=np.nanmean(even1)
            compMJO2ndj[i-220,j]=np.nanmean(even2)
            compMJO3ndj[i-220,j]=np.nanmean(even3)
            compMJO4ndj[i-220,j]=np.nanmean(even4)
            
            
            compMJO1ndj_all[i-220,j,:]=even1
            compMJO2ndj_all[i-220,j,:]=even2
            compMJO3ndj_all[i-220,j,:]=even3
            compMJO4ndj_all[i-220,j,:]=even4
            
    
 
#================================================================   
    
    #=======
    month='JFM'
    month2='JFM'
    
    compMJO1jfm=np.zeros((80,215))
    compMJO2jfm=np.zeros((80,215))
    compMJO3jfm=np.zeros((80,215))
    compMJO4jfm=np.zeros((80,215))
    
    compMJO1jfm_all=np.zeros((80,215,502))
    compMJO2jfm_all=np.zeros((80,215,502))
    compMJO3jfm_all=np.zeros((80,215,502))
    compMJO4jfm_all=np.zeros((80,215,502))
    
 
    
    #index--winter
    mjo=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/mjo8phase_amp1982-2021_'+month+'_nowk1lead5.txt')
    mjoph=mjo[:,0]
    amp=mjo[:,1]
    
    m1=np.copy(mjoph)
    m2=np.copy(mjoph)
    m3=np.copy(mjoph)
    m4=np.copy(mjoph)
    
    
    
    m1[(m1 >= 2) & (m1 <= 7)] = np.nan
    m1[amp<1]= np.nan
    
    m2[(m2 != 2) & (m2 != 3)] = np.nan
    m2[amp<1]= np.nan
    
    m3[(m3 != 4) & (m3 != 5)] = np.nan
    m3[amp<1]= np.nan
    
    m4[(m4 != 6) & (m4 != 7)] = np.nan
    m4[amp<1]= np.nan
    
    
    arwinter=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
    climatjfm=np.nanmean(arwinter,2)

    #--1982-2021 ndjfm
    for i in np.arange(220,300):
        
                
        for j in np.arange(0,215):
            
            e1=arwinter[i-220,j,1:503] #no week1 MJO
 
            even1=np.copy(e1)
            even2=np.copy(even1)
            even3=np.copy(even1)
            even4=np.copy(even1)
            
            even1[np.isnan(m1[0:np.shape(even1)[0]])]=np.nan
            even2[np.isnan(m2[0:np.shape(even1)[0]])]=np.nan
            even3[np.isnan(m3[0:np.shape(even1)[0]])]=np.nan
            even4[np.isnan(m4[0:np.shape(even1)[0]])]=np.nan
            
        
            
            compMJO1jfm[i-220,j]=np.nanmean(even1)
            compMJO2jfm[i-220,j]=np.nanmean(even2)
            compMJO3jfm[i-220,j]=np.nanmean(even3)
            compMJO4jfm[i-220,j]=np.nanmean(even4)
         
            
            compMJO1jfm_all[i-220,j,:]=even1
            compMJO2jfm_all[i-220,j,:]=even2
            compMJO3jfm_all[i-220,j,:]=even3
            compMJO4jfm_all[i-220,j,:]=even4
            
 
#================================================================   
    
    
    
    #---------------
    fig=plt.figure(1)#
    ax = fig.add_axes([h[0],z[0],0.8,0.8])#  125 59     0.1,0.1
    m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)
    m.drawstates(linewidth=0.5)
    
    x,y=np.meshgrid(lon,lat)
    x1,y1=m(x,y) 
    
    cs=m.contourf(x1,y1,compMJO1ndj-climatndj-(compMJO1jfm-climatjfm),cslev,cmap=cmap,extend='both')
    plt.title('MJO1&8_'+str(wind)+'week ARs NDJ-JFM',fontsize=fs)
    
    for i in arange(0,np.shape(lat)[0],2):
        for j in arange(0,np.shape(lon)[0],2):
            
            lin1=compMJO1jfm_all[i,j,:]
            lin2=compMJO1ndj_all[i,j,:]
            
            lin1=lin1[~np.isnan(lin1)]
            lin2=lin2[~np.isnan(lin2)]
            
            repeated_lin1 = climatjfm[i, j]
            repeated_lin2 = climatndj[i, j]
            
            lin1=lin1-repeated_lin1
            lin2=lin2-repeated_lin2
            
            
            
            t_statistic, sig = stats.ttest_ind(lin1, lin2)
            
            if (sig < 0.05 ):
                m.plot(x1[i,j],y1[i,j],'+',color='k',markersize=s,zorder=5)
                del sig
                                            
                    
    
    parallels = np.arange(20.,60,10.)
    m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[1,0,0,0],fontsize=flalo,zorder=2)
    meridians = np.arange(120.,260.,20.)
    m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,1],fontsize=flalo,zorder=2)
    
    ax.spines['bottom'].set_linewidth(linew)
    ax.spines['top'].set_linewidth(linew)
    ax.spines['left'].set_linewidth(linew)
    ax.spines['right'].set_linewidth(linew)
    
  
    
    
    ax = fig.add_axes([h[1],z[1],0.8,0.8])#  125 59   1.02,0.1
    m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)
    m.drawstates(linewidth=0.5)
    
    x,y=np.meshgrid(lon,lat)
    x1,y1=m(x,y) 
    
    cs=m.contourf(x1,y1,compMJO2ndj-climatndj-(compMJO2jfm-climatjfm),cslev,cmap=cmap,extend='both')
    plt.title('MJO2&3_'+str(wind)+'week ARs NDJ-JFM',fontsize=fs)
    
    for i in arange(0,np.shape(lat)[0],2):
        for j in arange(0,np.shape(lon)[0],2):
            
            lin1=compMJO2jfm_all[i,j,:]
            lin2=compMJO2ndj_all[i,j,:]
            
            lin1=lin1[~np.isnan(lin1)]
            lin2=lin2[~np.isnan(lin2)]
            
            repeated_lin1 = climatjfm[i, j]
            repeated_lin2 = climatndj[i, j]
            
            lin1=lin1-repeated_lin1
            lin2=lin2-repeated_lin2
            
            
            
            t_statistic, sig = stats.ttest_ind(lin1, lin2)
            
            if (sig < 0.05 ):
                m.plot(x1[i,j],y1[i,j],'+',color='k',markersize=s,zorder=5)
                del sig                 
    
    parallels = np.arange(20.,60,10.)
    m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
    meridians = np.arange(120.,260.,20.)
    m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,1],fontsize=flalo,zorder=2)
    
    ax.spines['bottom'].set_linewidth(linew)
    ax.spines['top'].set_linewidth(linew)
    ax.spines['left'].set_linewidth(linew)
    ax.spines['right'].set_linewidth(linew)
    #cbar1=fig.colorbar(cs,cax=plt.axes([1.75, 0.3, 0.03, 2.9])) #51
    #cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
    #cbar1.set_ticks(cslev)
    #cbar1.set_ticklabels((cslev))
    
 
    ax = fig.add_axes([h[2],z[2],0.8,0.8])#  125 59     0.1,0.1
    m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)
    m.drawstates(linewidth=0.5)
    
    x,y=np.meshgrid(lon,lat)
    x1,y1=m(x,y) 
    
    cs=m.contourf(x1,y1,compMJO3ndj-climatndj-(compMJO3jfm-climatjfm),cslev,cmap=cmap,extend='both')
    plt.title('MJO4&5_'+str(wind)+'week ARs NDJ-JFM',fontsize=fs)
    
    for i in arange(0,np.shape(lat)[0],2):
        for j in arange(0,np.shape(lon)[0],2):
            
            lin1=compMJO3jfm_all[i,j,:]
            lin2=compMJO3ndj_all[i,j,:]
            
            lin1=lin1[~np.isnan(lin1)]
            lin2=lin2[~np.isnan(lin2)]
            
            repeated_lin1 = climatjfm[i, j]
            repeated_lin2 = climatndj[i, j]
            
            lin1=lin1-repeated_lin1
            lin2=lin2-repeated_lin2
            
            
            
            t_statistic, sig = stats.ttest_ind(lin1, lin2)
            
            if (sig < 0.05 ):
                m.plot(x1[i,j],y1[i,j],'+',color='k',markersize=s,zorder=5)
                del sig 
    
    
    parallels = np.arange(20.,60,10.)
    m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[1,0,0,0],fontsize=flalo,zorder=2)
    meridians = np.arange(120.,260.,20.)
    m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,1],fontsize=flalo,zorder=2)
    
    ax.spines['bottom'].set_linewidth(linew)
    ax.spines['top'].set_linewidth(linew)
    ax.spines['left'].set_linewidth(linew)
    ax.spines['right'].set_linewidth(linew)
    
    
    
    
    
    ax = fig.add_axes([h[3],z[3],0.8,0.8])#  125 59     0.1,0.1
    m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)
    m.drawstates(linewidth=0.5)
    
    x,y=np.meshgrid(lon,lat)
    x1,y1=m(x,y) 
    
    cs=m.contourf(x1,y1,compMJO4ndj-climatndj-(compMJO4jfm-climatjfm),cslev,cmap=cmap,extend='both')
    plt.title('MJO6&7_'+str(wind)+'week ARs NDJ-JFM',fontsize=fs)
    
    for i in arange(0,np.shape(lat)[0],2):
        for j in arange(0,np.shape(lon)[0],2):
            
            lin1=compMJO4jfm_all[i,j,:]
            lin2=compMJO4ndj_all[i,j,:]
            
            lin1=lin1[~np.isnan(lin1)]
            lin2=lin2[~np.isnan(lin2)]
            
            repeated_lin1 = climatjfm[i, j]
            repeated_lin2 = climatndj[i, j]
            
            lin1=lin1-repeated_lin1
            lin2=lin2-repeated_lin2
            
            
            
            t_statistic, sig = stats.ttest_ind(lin1, lin2)
            
            if (sig < 0.05 ):
                m.plot(x1[i,j],y1[i,j],'+',color='k',markersize=s,zorder=5)
                del sig
    
      
    parallels = np.arange(20.,60,10.)
    m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
    meridians = np.arange(120.,260.,20.)
    m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,1],fontsize=flalo,zorder=2)
    
    ax.spines['bottom'].set_linewidth(linew)
    ax.spines['top'].set_linewidth(linew)
    ax.spines['left'].set_linewidth(linew)
    ax.spines['right'].set_linewidth(linew)
    
    
    
    ax.spines['bottom'].set_linewidth(linew)
    ax.spines['top'].set_linewidth(linew)
    ax.spines['left'].set_linewidth(linew)
    ax.spines['right'].set_linewidth(linew)
    
    cbar1=fig.colorbar(cs,cax=plt.axes([1.75, 1.3, 0.03, 0.9])) #51
    cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
    cbar1.set_ticks(cslev)
    cbar1.set_ticklabels((cslev))
    
  
    
    savefig('CompoMJO4window_'+str(wind)+'week_1982-2021indp_diff_'+month+'_lead5_bootstrap_NDJ-JFM_v3nooverlap-response.pdf',bbox_inches = 'tight')
    
    plt.show()
    
  
    
    
    












