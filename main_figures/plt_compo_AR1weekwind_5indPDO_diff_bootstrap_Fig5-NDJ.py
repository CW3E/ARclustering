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




col=(0.85,0.85,0.85)
linew=0.2
fs=16
fsc=17
s=1.5
flalo=12
cmap=cm.seismic
#'PNA'AO'
namesss=['PNA','AO','QBO','ENSO','RMM1','RMM2','PDO']

allcomphi=np.zeros((80,215,7))
allcomplo=np.zeros((80,215,7))

month='NDJ' #'JFM' nov2mar NDJ
month2='NDJ'  #NDJFM NDJ 
month3='NDJ'  #NDJFM NDJ 

h=[0.5*6+0.1,0.5*5+0.1,0.5*4+0.1,0.5*3+0.1,0.5*2+0.1,0.5+0.1,0.1]
h2=[0.5*6+0.1,0.5*5+0.1,0.5*4+0.1,0.5*3+0.1,0.5*2+0.1,0.5+0.1,0.1]

cslev=[-0.8,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,-0.05,-0.01,0,0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.8]

for wind in np.arange(1,2): #week window
  
    
    obj=nc.Dataset('globalARcatalog_MERRA2_1980-2021_v3.0all.nc')
    
    lat=obj.variables['lat'][220:300]  #30-50N
    lon=obj.variables['lon'][200:415]  #230-245
    
    climat=np.load('all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
    climat=np.nanmean(climat,2)

    
    #===========
    name='PNA'
    ind5_hi=np.loadtxt('pna1982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    ind5_lo=np.loadtxt('pna1982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    
    ind5_hi[ind5_hi<0]=np.nan
    ind5_lo[ind5_lo>0]=np.nan    
    
    comp_ind5_hi=np.load('all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
    comp_ind5_lo=np.load('all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
        
    comp_ind5_hi[:,:,np.isnan(ind5_hi)]=np.nan
    comp_ind5_lo[:,:,np.isnan(ind5_lo)]=np.nan
    
    comp_ind5_hi=np.nanmean(comp_ind5_hi,2)
    comp_ind5_lo=np.nanmean(comp_ind5_lo,2)
    
    allcomphi[:,:,0]=comp_ind5_hi
    allcomplo[:,:,0]=comp_ind5_lo
    #===========

    #===========
    name='AO'
    ind5_hi=np.loadtxt('ao1982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    ind5_lo=np.loadtxt('ao1982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    
    ind5_hi[ind5_hi<0]=np.nan
    ind5_lo[ind5_lo>0]=np.nan    
    
    comp_ind5_hi=np.load('all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
    comp_ind5_lo=np.load('all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
        
    comp_ind5_hi[:,:,np.isnan(ind5_hi)]=np.nan
    comp_ind5_lo[:,:,np.isnan(ind5_lo)]=np.nan
    
    comp_ind5_hi=np.nanmean(comp_ind5_hi,2)
    comp_ind5_lo=np.nanmean(comp_ind5_lo,2)
    
    allcomphi[:,:,1]=comp_ind5_hi
    allcomplo[:,:,1]=comp_ind5_lo
    #===========
    name='QBO'
    ind5_hi=np.loadtxt('QBO1982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    ind5_lo=np.loadtxt('QBO1982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    
    ind5_hi[ind5_hi<0]=np.nan
    ind5_lo[ind5_lo>0]=np.nan    
    
    comp_ind5_hi=np.load('all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
    comp_ind5_lo=np.load('all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
        
    comp_ind5_hi[:,:,np.isnan(ind5_hi)]=np.nan
    comp_ind5_lo[:,:,np.isnan(ind5_lo)]=np.nan
    
    comp_ind5_hi=np.nanmean(comp_ind5_hi,2)
    comp_ind5_lo=np.nanmean(comp_ind5_lo,2)
    
    allcomphi[:,:,2]=comp_ind5_hi
    allcomplo[:,:,2]=comp_ind5_lo
    
    #===========
    name='ENSO'
    ind5_hi=np.loadtxt('nino1982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    ind5_lo=np.loadtxt('nino1982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    
    ind5_hi[ind5_hi<0]=np.nan
    ind5_lo[ind5_lo>0]=np.nan    
    
    comp_ind5_hi=np.load('all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
    comp_ind5_lo=np.load('all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
        
    comp_ind5_hi[:,:,np.isnan(ind5_hi)]=np.nan
    comp_ind5_lo[:,:,np.isnan(ind5_lo)]=np.nan
    
    comp_ind5_hi=np.nanmean(comp_ind5_hi,2)
    comp_ind5_lo=np.nanmean(comp_ind5_lo,2)
    
    allcomphi[:,:,3]=comp_ind5_hi
    allcomplo[:,:,3]=comp_ind5_lo

    #===========
    name='RMM1'
    ind5_hi=np.loadtxt('rmm11982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    ind5_lo=np.loadtxt('rmm11982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    
    ind5_hi[ind5_hi<0]=np.nan
    ind5_lo[ind5_lo>0]=np.nan    
    
    comp_ind5_hi=np.load('all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
    comp_ind5_lo=np.load('all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
        
    comp_ind5_hi[:,:,np.isnan(ind5_hi)]=np.nan
    comp_ind5_lo[:,:,np.isnan(ind5_lo)]=np.nan
    
    comp_ind5_hi=np.nanmean(comp_ind5_hi,2)
    comp_ind5_lo=np.nanmean(comp_ind5_lo,2)
    
    allcomphi[:,:,4]=comp_ind5_hi
    allcomplo[:,:,4]=comp_ind5_lo
    
    #===========
    name='RMM2'
    ind5_hi=np.loadtxt('rmm21982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    ind5_lo=np.loadtxt('rmm21982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    
    ind5_hi[ind5_hi<0]=np.nan
    ind5_lo[ind5_lo>0]=np.nan    
    
    comp_ind5_hi=np.load('all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
    comp_ind5_lo=np.load('all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
        
    comp_ind5_hi[:,:,np.isnan(ind5_hi)]=np.nan
    comp_ind5_lo[:,:,np.isnan(ind5_lo)]=np.nan
    
    comp_ind5_hi=np.nanmean(comp_ind5_hi,2)
    comp_ind5_lo=np.nanmean(comp_ind5_lo,2)
    
    allcomphi[:,:,5]=comp_ind5_hi
    allcomplo[:,:,5]=comp_ind5_lo


    #_------xzx1
    
    allclimate=np.load('comprad_AR'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')

   
    sig=np.zeros((np.shape(lat)[0],np.shape(lon)[0],6))   
    sig2=np.zeros((np.shape(lat)[0],np.shape(lon)[0],6))  
    
    
    
    for i in arange(0,np.shape(lat)[0]):
       for j in arange(0,np.shape(lon)[0]): 
           
               amax=np.percentile(allclimate[i,j,:],95)
               amin=np.percentile(allclimate[i,j,:],5)
               
               for nam_n in arange(0,6):
               
                   if ((allcomphi[i,j,nam_n]<amax) and (allcomphi[i,j,nam_n]>amin)):
                       sig[i,j,nam_n]=np.nan
                   elif  (allcomphi[i,j,nam_n]>=amax) :
                       sig[i,j,nam_n]=1
                       
                   
                   if ((allcomplo[i,j,nam_n]<amax) and (allcomplo[i,j,nam_n]>amin)):
                       sig2[i,j,nam_n]=np.nan
                   elif  (allcomplo[i,j,nam_n]>=amax) :
                       sig2[i,j,nam_n]=1
        
#===========pdo======!!!!!

    obj=nc.Dataset('era5_ar_gwv3_1940-2022.nc')
    
    latpdo=obj.variables['lat'][120:281]  #60-20N
    lonpdo=obj.variables['lon'][480:1041]  #120-260
    
    climatpdo=np.load('ERA5all_num'+str(wind)+'week_indp_'+month+'_1940alldomain_bootstrap_nooverlap.npy')
    climatpdo=np.nanmean(climatpdo,2)
    
    
    #===========
    name='pdo'
    ind5_hi=np.loadtxt('pdo1940-2018_'+month2+'_ave1wk.txt')
    ind5_lo=np.loadtxt('pdo1940-2018_'+month2+'_ave1wk.txt')
    
    ind5_hi[ind5_hi<0]=np.nan
    ind5_lo[ind5_lo>0]=np.nan    
    
    comp_ind5_hi=np.load('ERA5all_num'+str(wind)+'week_indp_'+month+'_1940alldomain_bootstrap_nooverlap.npy')
    comp_ind5_lo=np.load('ERA5all_num'+str(wind)+'week_indp_'+month+'_1940alldomain_bootstrap_nooverlap.npy')
        
    comp_ind5_hi[:,:,np.isnan(ind5_hi)]=np.nan
    comp_ind5_lo[:,:,np.isnan(ind5_lo)]=np.nan
    
    comp_ind5_hi=np.nanmean(comp_ind5_hi,2)
    comp_ind5_lo=np.nanmean(comp_ind5_lo,2)
    
    allcomphipdo=comp_ind5_hi
    allcomplopdo=comp_ind5_lo
    #===========
    #_------xzxpdo
    
    allclimate=np.load('ERA5comprad_AR'+str(wind)+'week_indp_'+month+'_1940alldomain_bootstrap_nooverlap.npy')
    
    
    sigpdo=np.zeros((np.shape(latpdo)[0],np.shape(lonpdo)[0]))   
    sig2pdo=np.zeros((np.shape(latpdo)[0],np.shape(lonpdo)[0]))  
    


    for i in arange(0,np.shape(latpdo)[0]):
       for j in arange(0,np.shape(lonpdo)[0]): 
           
            amax=np.percentile(allclimate[i,j,:],95)
            amin=np.percentile(allclimate[i,j,:],5)
            
               
        
            if ((allcomphipdo[i,j]<amax) and (allcomphipdo[i,j]>amin)):
                sigpdo[i,j]=np.nan
            elif  (allcomphipdo[i,j]>=amax) :
                sigpdo[i,j]=1
                
            
            if ((allcomplopdo[i,j]<amax) and (allcomplopdo[i,j]>amin)):
                sig2pdo[i,j]=np.nan
            elif  (allcomplopdo[i,j]>=amax) :
                sig2pdo[i,j]=1
        


    

                  
    #---------------
    
for nam_n in arange(0,7):    
    
    fig=plt.figure(1)#
    
    if nam_n==6:
        ax = fig.add_axes([0.1,h[nam_n-1],0.8,0.8])#  125 59
        m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        m.drawstates(linewidth=0.5)
        
        x,y=np.meshgrid(lonpdo,latpdo)
        x1,y1=m(x,y) 
        #
        cs=m.contourf(x1,y1,allcomphipdo-climatpdo,cslev,cmap=cmap,extend='both')
        plt.title('Composite +'+namesss[nam_n]+'_'+str(wind)+'week ARs '+month3,fontsize=fs)
        
        for i in arange(0,np.shape(latpdo)[0],5):
            for j in arange(0,np.shape(lonpdo)[0],5):
                
                    if (sigpdo[i,j]==1 ):
                        m.plot(x1[i,j],y1[i,j],'+',color='k',markersize=s,zorder=5)
                    if (sigpdo[i,j]==0 ):   
                        m.plot(x1[i,j],y1[i,j],'+',color='k',markersize=s,zorder=5)
                       
                            
        parallels = np.arange(20.,60,10.)
        m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[1,0,0,0],fontsize=flalo,zorder=2)
        meridians = np.arange(120.,260.,20.)
        m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
        
        ax.spines['bottom'].set_linewidth(linew)
        ax.spines['top'].set_linewidth(linew)
        ax.spines['left'].set_linewidth(linew)
        ax.spines['right'].set_linewidth(linew)
        
        #---------------
        
        
        ax = fig.add_axes([0.92,h2[nam_n-1],0.8,0.8])#  125 59
        m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        m.drawstates(linewidth=0.5)
        
        x,y=np.meshgrid(lonpdo,latpdo)
        x1,y1=m(x,y) 
        #
        cs=m.contourf(x1,y1,allcomplopdo-climatpdo,cslev,cmap=cmap,extend='both')
        
        for i in arange(0,np.shape(latpdo)[0],5):
            for j in arange(0,np.shape(lonpdo)[0],5):
                
                    if (sig2pdo[i,j]==1 ):
                        m.plot(x1[i,j],y1[i,j],'+',color='k',markersize=s,zorder=5)
                    if (sig2pdo[i,j]==0 ):   
                        m.plot(x1[i,j],y1[i,j],'+',color='k',markersize=s,zorder=5)
                        
        plt.title('Composite -'+namesss[nam_n]+'_'+str(wind)+'week ARs '+month3,fontsize=fs)
        
        parallels = np.arange(20.,60,10.)
        m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
        meridians = np.arange(120.,260.,20.)
        m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
        
        ax.spines['bottom'].set_linewidth(linew)
        ax.spines['top'].set_linewidth(linew)
        ax.spines['left'].set_linewidth(linew)
        ax.spines['right'].set_linewidth(linew)
        cbar1=fig.colorbar(cs,cax=plt.axes([1.75, 0.3, 0.03, 3.4])) #51
        cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
        cbar1.set_ticks(cslev)
        cbar1.set_ticklabels((cslev))
         
    
    
    
    else:
    
        ax = fig.add_axes([0.1,h[nam_n-1],0.8,0.8])#  125 59
        m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        m.drawstates(linewidth=0.5)
        
        x,y=np.meshgrid(lon,lat)
        x1,y1=m(x,y) 
        #
        cs=m.contourf(x1,y1,allcomphi[:,:,nam_n]-climat,cslev,cmap=cmap,extend='both')
        plt.title('Composite +'+namesss[nam_n]+'_'+str(wind)+'week ARs '+month2,fontsize=fs)
        
        for i in arange(0,np.shape(lat)[0],3):
            for j in arange(0,np.shape(lon)[0],3):
                
                    if (sig[i,j,nam_n]==1 ):
                        m.plot(x1[i,j],y1[i,j],'+',color='k',markersize=s,zorder=5)
                    if (sig[i,j,nam_n]==0 ):   
                        m.plot(x1[i,j],y1[i,j],'+',color='k',markersize=s,zorder=5)
                       
                            
        parallels = np.arange(20.,60,10.)
        m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[1,0,0,0],fontsize=flalo,zorder=2)
        meridians = np.arange(120.,260.,20.)
        m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
        
        if namesss[nam_n]=='PNA':
            m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,1],fontsize=flalo,zorder=2)
        
        ax.spines['bottom'].set_linewidth(linew)
        ax.spines['top'].set_linewidth(linew)
        ax.spines['left'].set_linewidth(linew)
        ax.spines['right'].set_linewidth(linew)
        
       
        
        
        ax = fig.add_axes([0.92,h2[nam_n-1],0.8,0.8])#  125 59
        m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        m.drawstates(linewidth=0.5)
        
        x,y=np.meshgrid(lon,lat)
        x1,y1=m(x,y) 
        #
        cs=m.contourf(x1,y1,allcomplo[:,:,nam_n]-climat,cslev,cmap=cmap,extend='both')
        
        for i in arange(0,np.shape(lat)[0],3):
            for j in arange(0,np.shape(lon)[0],3):
                
                    if (sig2[i,j,nam_n]==1 ):
                        m.plot(x1[i,j],y1[i,j],'+',color='k',markersize=s,zorder=5)
                    if (sig2[i,j,nam_n]==0 ):   
                        m.plot(x1[i,j],y1[i,j],'+',color='k',markersize=s,zorder=5)
                        
        plt.title('Composite -'+namesss[nam_n]+'_'+str(wind)+'week ARs '+month2,fontsize=fs)
        
        parallels = np.arange(20.,60,10.)
        m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
        meridians = np.arange(120.,260.,20.)
        m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
        
        if namesss[nam_n]=='PNA':
            m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,1],fontsize=flalo,zorder=2)

        
        ax.spines['bottom'].set_linewidth(linew)
        ax.spines['top'].set_linewidth(linew)
        ax.spines['left'].set_linewidth(linew)
        ax.spines['right'].set_linewidth(linew)
        
        
        
savefig('Compo-window_1week_1982-2021indp_diff_'+month+'_5indPDO_bootstrap_v3_nooverlap-response.pdf',bbox_inches = 'tight')
    
plt.show()
    
     



