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
fsc=12
s=1
flalo=10

cmap=cm.seismic


for wind in np.arange(1,2): #week window
 
    
    if wind==1:
        cslev=[-1,-0.8,-0.6,-0.4,-0.2,-0.1,-0.02,0,0.02,0.1,0.2,0.4,0.6,0.8,1]
    elif wind==2:
        cslev=[0,0.5,1,1.5,2,2.5,3]
    elif wind==3:
        cslev=[0,0.5,1,1.5,2,2.5,3,3.5,4,4.5]
    elif wind==4:
        cslev=[0,1,2,3,4,5,6]
    elif wind==5:
        cslev=[0,1,2,3,4,5,6,7,8]
    elif wind==6:
        cslev=[0,2,4,6,8,10]
    
    
    
    
    obj=nc.Dataset('/Users/zhiqiyang/Downloads/globalARcatalog_MERRA2_1980-2021_v3.0all.nc')
    
    lat=obj.variables['lat'][220:300]  #30-50N
    lon=obj.variables['lon'][200:415]  #230-245
    
    
        
    compMJO1=np.zeros((80,215))
    compMJO2=np.zeros((80,215))
    compMJO3=np.zeros((80,215))
    compMJO4=np.zeros((80,215))
    
    
    
    #index--winter
    mjo=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/mjo8phase_amp1982-2021_ndjfm_nowk1lead5.txt')
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
    
    
    arwinter=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/all_num'+str(wind)+'week_indp_nov2mar_1982alldomain_bootstrap_nooverlap.npy')
    climat=np.nanmean(arwinter,2)

    #--1982-2021 ndjfm
    for i in np.arange(220,300):
        
                
        for j in np.arange(0,215):
            
            e1=arwinter[i-220,j,1:843] #no week1 MJO
 
            even1=np.copy(e1)
            even2=np.copy(even1)
            even3=np.copy(even1)
            even4=np.copy(even1)
            
            even1[np.isnan(m1[0:np.shape(even1)[0]])]=np.nan
            even2[np.isnan(m2[0:np.shape(even1)[0]])]=np.nan
            even3[np.isnan(m3[0:np.shape(even1)[0]])]=np.nan
            even4[np.isnan(m4[0:np.shape(even1)[0]])]=np.nan
            
            
            
            #evenhi=np.array(calculate_average_ar_number(e1hi,wind))
            #evenlo=np.array(calculate_average_ar_number(e2lo,wind))
            
            
            #evenhi[np.isnan(hiwindow)]=np.nan
            #evenhi=evenhi[~np.isnan(evenhi)]
            
            
        
            
            compMJO1[i-220,j]=np.nanmean(even1)
            compMJO2[i-220,j]=np.nanmean(even2)
            compMJO3[i-220,j]=np.nanmean(even3)
            compMJO4[i-220,j]=np.nanmean(even4)
            
 
#_------xzx1
allclimate=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/comprad_AR'+str(wind)+'week_indp_nov2mar_1982alldomain_bootstrap_nooverlap.npy')


mtsig1=np.zeros((np.shape(lat)[0],np.shape(lon)[0]))  
mtsig2=np.zeros((np.shape(lat)[0],np.shape(lon)[0]))  
mtsig3=np.zeros((np.shape(lat)[0],np.shape(lon)[0]))  
mtsig4=np.zeros((np.shape(lat)[0],np.shape(lon)[0]))  



for i in arange(0,np.shape(lat)[0]):
   for j in arange(0,np.shape(lon)[0]): 
       
           amax=np.percentile(allclimate[i,j,:],95)
           amin=np.percentile(allclimate[i,j,:],5)
           
           if ((compMJO1[i,j]<amax) and (compMJO1[i,j]>amin)):
               mtsig1[i,j]=np.nan
           elif  (compMJO1[i,j]>=amax) :
               mtsig1[i,j]=1
               
               
           if ((compMJO2[i,j]<amax) and (compMJO2[i,j]>amin)):
               mtsig2[i,j]=np.nan
           else :
               mtsig2[i,j]=1    
               
           if ((compMJO3[i,j]<amax) and (compMJO3[i,j]>amin)):
               mtsig3[i,j]=np.nan
           else :
               mtsig3[i,j]=1    
               
               
           if ((compMJO4[i,j]<amax) and (compMJO4[i,j]>amin)):
               mtsig4[i,j]=np.nan
           else :
               mtsig4[i,j]=1       
               
            
    

#---------------
fig=plt.figure(1)#
ax = fig.add_axes([h[0],z[0],0.8,0.8])#  125 59     0.1,0.1
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
m.drawstates(linewidth=0.5)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 

cs=m.contourf(x1,y1,compMJO1-climat,cslev,cmap=cmap,extend='both')
plt.title('Composite lead5days MJO1&8_'+str(wind)+'week ARs (NDJFM)',fontsize=fs)

for i in arange(0,np.shape(lat)[0],2):
    for j in arange(0,np.shape(lon)[0],2):
        
            if (~np.isnan(mtsig1[i,j])):
                m.plot(x1[i,j],y1[i,j],'+',color='k',markersize=s,zorder=5)


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




ax = fig.add_axes([h[1],z[1],0.8,0.8])#  125 59   1.02,0.1
m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
m.drawcoastlines(linewidth=0.5)
m.drawcountries(linewidth=0.5)
m.drawstates(linewidth=0.5)

x,y=np.meshgrid(lon,lat)
x1,y1=m(x,y) 

cs=m.contourf(x1,y1,compMJO2-climat,cslev,cmap=cmap,extend='both')
plt.title('Composite lead5days MJO2&3_'+str(wind)+'week ARs (NDJFM)',fontsize=fs)

for i in arange(0,np.shape(lat)[0],2):
    for j in arange(0,np.shape(lon)[0],2):
        
            if (mtsig2[i,j]==1 ):
                m.plot(x1[i,j],y1[i,j],'+',color='k',markersize=s,zorder=5)
            if (mtsig2[i,j]==0 ):   
                m.plot(x1[i,j],y1[i,j],'+',color='k',markersize=s,zorder=5)

parallels = np.arange(20.,60,10.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(120.,260.,20.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,1],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)
#cbar1=fig.colorbar(cs,cax=plt.axes([1.84, 0.3, 0.01, 0.4])) #51
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

cs=m.contourf(x1,y1,compMJO3-climat,cslev,cmap=cmap,extend='both')
plt.title('Composite lead5days MJO4&5_'+str(wind)+'week ARs (NDJFM)',fontsize=fs)

for i in arange(0,np.shape(lat)[0],2):
    for j in arange(0,np.shape(lon)[0],2):
        
            if (mtsig3[i,j]==1 ):
                m.plot(x1[i,j],y1[i,j],'+',color='k',markersize=s,zorder=5)
            if (mtsig3[i,j]==0 ):   
                m.plot(x1[i,j],y1[i,j],'+',color='k',markersize=s,zorder=5)


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

cs=m.contourf(x1,y1,compMJO4-climat,cslev,cmap=cmap,extend='both')
plt.title('Composite lead5days MJO6&7_'+str(wind)+'week ARs (NDJFM)',fontsize=fs)

for i in arange(0,np.shape(lat)[0],2):
    for j in arange(0,np.shape(lon)[0],2):
        
            if (mtsig4[i,j]==1 ):
                m.plot(x1[i,j],y1[i,j],'+',color='k',markersize=s,zorder=5)
            if (mtsig4[i,j]==0 ):   
                m.plot(x1[i,j],y1[i,j],'+',color='k',markersize=s,zorder=5)


parallels = np.arange(20.,60,10.)
m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
meridians = np.arange(120.,260.,20.)
m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,1],fontsize=flalo,zorder=2)

ax.spines['bottom'].set_linewidth(linew)
ax.spines['top'].set_linewidth(linew)
ax.spines['left'].set_linewidth(linew)
ax.spines['right'].set_linewidth(linew)





cbar1=fig.colorbar(cs,cax=plt.axes([1.75, 1.3, 0.03, 0.9])) #51
cbar1.ax.tick_params(labelsize=fsc, width=linew,length=1)
cbar1.set_ticks(cslev)
cbar1.set_ticklabels((cslev))




savefig('CompoMJO4comb_'+str(wind)+'week_1982-2021indp_diff_amp1lead5_nov2mar_nooverlap_v4_response.pdf',bbox_inches = 'tight')
plt.show()

  




 










