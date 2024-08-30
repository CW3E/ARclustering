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
cmap=cm.PuOr_r
#'PNA'AO'
namesss=['PNA','AO','QBO','ENSO','RMM1','RMM2','PDO']

allcomphiNDJ=np.zeros((80,215,6))
allcomploNDJ=np.zeros((80,215,6))

allcomphiJFM=np.zeros((80,215,6))
allcomploJFM=np.zeros((80,215,6))

NDJ_hi_allweek=np.zeros((80,215,513,6))
NDJ_lo_allweek=np.zeros((80,215,513,6))

JFM_hi_allweek=np.zeros((80,215,503,6))
JFM_lo_allweek=np.zeros((80,215,503,6))


h=[0.5*6+0.1,0.5*5+0.1,0.5*4+0.1,0.5*3+0.1,0.5*2+0.1,0.5+0.1,0.1]
h2=[0.5*6+0.1,0.5*5+0.1,0.5*4+0.1,0.5*3+0.1,0.5*2+0.1,0.5+0.1,0.1]

cslev=[-0.5,-0.4,-0.3,-0.2,-0.1,-0.05,-0.01,0,0.01,0.05,0.1,0.2,0.3,0.4,0.5]

for wind in np.arange(1,2): #week window
  
    month='NDJ' #'JFM' nov2mar NDJ
    month2='NDJ'  #NDJFM NDJ 
    obj=nc.Dataset('/Users/zhiqiyang/Downloads/globalARcatalog_MERRA2_1980-2021_v3.0all.nc')
    
    lat=obj.variables['lat'][220:300]  #30-50N
    lon=obj.variables['lon'][200:415]  #230-245
    
    climat=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
    climatNDJ=np.nanmean(climat,2)

    
    #===========
    name='PNA'
    ind5_hi=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/pna1982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    ind5_lo=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/pna1982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    
    ind5_hi[ind5_hi<0]=np.nan
    ind5_lo[ind5_lo>0]=np.nan    
    
    comp_ind5_hi=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
    comp_ind5_lo=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
        
    comp_ind5_hi[:,:,np.isnan(ind5_hi)]=np.nan
    comp_ind5_lo[:,:,np.isnan(ind5_lo)]=np.nan
    
    NDJ_hi_allweek[:,:,:,0]=np.copy(comp_ind5_hi)
    NDJ_lo_allweek[:,:,:,0]=np.copy(comp_ind5_lo)
    
    
    comp_ind5_hi=np.nanmean(comp_ind5_hi,2)
    comp_ind5_lo=np.nanmean(comp_ind5_lo,2)
    
    allcomphiNDJ[:,:,0]=comp_ind5_hi
    allcomploNDJ[:,:,0]=comp_ind5_lo
    #===========

    #===========
    name='AO'
    ind5_hi=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/ao1982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    ind5_lo=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/ao1982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    
    ind5_hi[ind5_hi<0]=np.nan
    ind5_lo[ind5_lo>0]=np.nan    
    
    comp_ind5_hi=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
    comp_ind5_lo=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
        
    comp_ind5_hi[:,:,np.isnan(ind5_hi)]=np.nan
    comp_ind5_lo[:,:,np.isnan(ind5_lo)]=np.nan
    
    NDJ_hi_allweek[:,:,:,1]=np.copy(comp_ind5_hi)
    NDJ_lo_allweek[:,:,:,1]=np.copy(comp_ind5_lo)
    
    
    comp_ind5_hi=np.nanmean(comp_ind5_hi,2)
    comp_ind5_lo=np.nanmean(comp_ind5_lo,2)
    
    allcomphiNDJ[:,:,1]=comp_ind5_hi
    allcomploNDJ[:,:,1]=comp_ind5_lo
    #===========
    name='QBO'
    ind5_hi=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/QBO1982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    ind5_lo=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/QBO1982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    
    ind5_hi[ind5_hi<0]=np.nan
    ind5_lo[ind5_lo>0]=np.nan    
    
    comp_ind5_hi=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
    comp_ind5_lo=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
        
    comp_ind5_hi[:,:,np.isnan(ind5_hi)]=np.nan
    comp_ind5_lo[:,:,np.isnan(ind5_lo)]=np.nan
    
    NDJ_hi_allweek[:,:,:,2]=np.copy(comp_ind5_hi)
    NDJ_lo_allweek[:,:,:,2]=np.copy(comp_ind5_lo)
    
    
    comp_ind5_hi=np.nanmean(comp_ind5_hi,2)
    comp_ind5_lo=np.nanmean(comp_ind5_lo,2)
    
    allcomphiNDJ[:,:,2]=comp_ind5_hi
    allcomploNDJ[:,:,2]=comp_ind5_lo
    
    #===========
    name='ENSO'
    ind5_hi=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/nino1982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    ind5_lo=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/nino1982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    
    ind5_hi[ind5_hi<0]=np.nan
    ind5_lo[ind5_lo>0]=np.nan    
    
    comp_ind5_hi=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
    comp_ind5_lo=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
        
    comp_ind5_hi[:,:,np.isnan(ind5_hi)]=np.nan
    comp_ind5_lo[:,:,np.isnan(ind5_lo)]=np.nan
    
    NDJ_hi_allweek[:,:,:,3]=np.copy(comp_ind5_hi)
    NDJ_lo_allweek[:,:,:,3]=np.copy(comp_ind5_lo)
    
    
    comp_ind5_hi=np.nanmean(comp_ind5_hi,2)
    comp_ind5_lo=np.nanmean(comp_ind5_lo,2)
    
    allcomphiNDJ[:,:,3]=comp_ind5_hi
    allcomploNDJ[:,:,3]=comp_ind5_lo

    #===========
    name='RMM1'
    ind5_hi=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/rmm11982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    ind5_lo=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/rmm11982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    
    ind5_hi[ind5_hi<0]=np.nan
    ind5_lo[ind5_lo>0]=np.nan    
    
    comp_ind5_hi=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
    comp_ind5_lo=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
        
    comp_ind5_hi[:,:,np.isnan(ind5_hi)]=np.nan
    comp_ind5_lo[:,:,np.isnan(ind5_lo)]=np.nan
    
    NDJ_hi_allweek[:,:,:,4]=np.copy(comp_ind5_hi)
    NDJ_lo_allweek[:,:,:,4]=np.copy(comp_ind5_lo)
    
    
    comp_ind5_hi=np.nanmean(comp_ind5_hi,2)
    comp_ind5_lo=np.nanmean(comp_ind5_lo,2)
    
    allcomphiNDJ[:,:,4]=comp_ind5_hi
    allcomploNDJ[:,:,4]=comp_ind5_lo
    
    #===========
    name='RMM2'
    ind5_hi=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/rmm21982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    ind5_lo=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/rmm21982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    
    ind5_hi[ind5_hi<0]=np.nan
    ind5_lo[ind5_lo>0]=np.nan    
    
    comp_ind5_hi=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
    comp_ind5_lo=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
        
    comp_ind5_hi[:,:,np.isnan(ind5_hi)]=np.nan
    comp_ind5_lo[:,:,np.isnan(ind5_lo)]=np.nan
    
    NDJ_hi_allweek[:,:,:,5]=np.copy(comp_ind5_hi)
    NDJ_lo_allweek[:,:,:,5]=np.copy(comp_ind5_lo)
    
    
    comp_ind5_hi=np.nanmean(comp_ind5_hi,2)
    comp_ind5_lo=np.nanmean(comp_ind5_lo,2)
    
    allcomphiNDJ[:,:,5]=comp_ind5_hi
    allcomploNDJ[:,:,5]=comp_ind5_lo


        
    
#========================================
    month='JFM' #'JFM' nov2mar NDJ
    month2='JFM'  #NDJFM NDJ 
    
    climat=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
    climatJFM=np.nanmean(climat,2)

    
    #===========
    name='PNA'
    ind5_hi=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/pna1982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    ind5_lo=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/pna1982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    
    ind5_hi[ind5_hi<0]=np.nan
    ind5_lo[ind5_lo>0]=np.nan    
    
    comp_ind5_hi=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
    comp_ind5_lo=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
        
    comp_ind5_hi[:,:,np.isnan(ind5_hi)]=np.nan
    comp_ind5_lo[:,:,np.isnan(ind5_lo)]=np.nan
    
    JFM_hi_allweek[:,:,:,0]=np.copy(comp_ind5_hi)
    JFM_lo_allweek[:,:,:,0]=np.copy(comp_ind5_lo)
    
    
    comp_ind5_hi=np.nanmean(comp_ind5_hi,2)
    comp_ind5_lo=np.nanmean(comp_ind5_lo,2)
    
    allcomphiJFM[:,:,0]=comp_ind5_hi
    allcomploJFM[:,:,0]=comp_ind5_lo
    #===========

    #===========
    name='AO'
    ind5_hi=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/ao1982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    ind5_lo=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/ao1982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    
    ind5_hi[ind5_hi<0]=np.nan
    ind5_lo[ind5_lo>0]=np.nan    
    
    comp_ind5_hi=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
    comp_ind5_lo=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
        
    comp_ind5_hi[:,:,np.isnan(ind5_hi)]=np.nan
    comp_ind5_lo[:,:,np.isnan(ind5_lo)]=np.nan
    
    JFM_hi_allweek[:,:,:,1]=np.copy(comp_ind5_hi)
    JFM_lo_allweek[:,:,:,1]=np.copy(comp_ind5_lo)
    
    
    comp_ind5_hi=np.nanmean(comp_ind5_hi,2)
    comp_ind5_lo=np.nanmean(comp_ind5_lo,2)
    
    allcomphiJFM[:,:,1]=comp_ind5_hi
    allcomploJFM[:,:,1]=comp_ind5_lo
    #===========
    name='QBO'
    ind5_hi=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/QBO1982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    ind5_lo=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/QBO1982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    
    ind5_hi[ind5_hi<0]=np.nan
    ind5_lo[ind5_lo>0]=np.nan    
    
    comp_ind5_hi=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
    comp_ind5_lo=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
        
    comp_ind5_hi[:,:,np.isnan(ind5_hi)]=np.nan
    comp_ind5_lo[:,:,np.isnan(ind5_lo)]=np.nan
    
    JFM_hi_allweek[:,:,:,2]=np.copy(comp_ind5_hi)
    JFM_lo_allweek[:,:,:,2]=np.copy(comp_ind5_lo)
    
    
    comp_ind5_hi=np.nanmean(comp_ind5_hi,2)
    comp_ind5_lo=np.nanmean(comp_ind5_lo,2)
    
    allcomphiJFM[:,:,2]=comp_ind5_hi
    allcomploJFM[:,:,2]=comp_ind5_lo
    
    #===========
    name='ENSO'
    ind5_hi=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/nino1982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    ind5_lo=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/nino1982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    
    ind5_hi[ind5_hi<0]=np.nan
    ind5_lo[ind5_lo>0]=np.nan    
    
    comp_ind5_hi=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
    comp_ind5_lo=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
        
    comp_ind5_hi[:,:,np.isnan(ind5_hi)]=np.nan
    comp_ind5_lo[:,:,np.isnan(ind5_lo)]=np.nan
    
    JFM_hi_allweek[:,:,:,3]=np.copy(comp_ind5_hi)
    JFM_lo_allweek[:,:,:,3]=np.copy(comp_ind5_lo)
    
    
    
    comp_ind5_hi=np.nanmean(comp_ind5_hi,2)
    comp_ind5_lo=np.nanmean(comp_ind5_lo,2)
    
    allcomphiJFM[:,:,3]=comp_ind5_hi
    allcomploJFM[:,:,3]=comp_ind5_lo

    #===========
    name='RMM1'
    ind5_hi=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/rmm11982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    ind5_lo=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/rmm11982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    
    ind5_hi[ind5_hi<0]=np.nan
    ind5_lo[ind5_lo>0]=np.nan    
    
    comp_ind5_hi=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
    comp_ind5_lo=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
        
    comp_ind5_hi[:,:,np.isnan(ind5_hi)]=np.nan
    comp_ind5_lo[:,:,np.isnan(ind5_lo)]=np.nan
    
    JFM_hi_allweek[:,:,:,4]=np.copy(comp_ind5_hi)
    JFM_lo_allweek[:,:,:,4]=np.copy(comp_ind5_lo)
    
    
    
    comp_ind5_hi=np.nanmean(comp_ind5_hi,2)
    comp_ind5_lo=np.nanmean(comp_ind5_lo,2)
    
    allcomphiJFM[:,:,4]=comp_ind5_hi
    allcomploJFM[:,:,4]=comp_ind5_lo
    
    #===========
    name='RMM2'
    ind5_hi=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/rmm21982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    ind5_lo=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/rmm21982-2021_'+month2+'_ave'+str(wind)+'wk.txt')
    
    ind5_hi[ind5_hi<0]=np.nan
    ind5_lo[ind5_lo>0]=np.nan    
    
    comp_ind5_hi=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
    comp_ind5_lo=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/all_num'+str(wind)+'week_indp_'+month+'_1982alldomain_bootstrap_nooverlap.npy')
        
    comp_ind5_hi[:,:,np.isnan(ind5_hi)]=np.nan
    comp_ind5_lo[:,:,np.isnan(ind5_lo)]=np.nan
    
    JFM_hi_allweek[:,:,:,5]=np.copy(comp_ind5_hi)
    JFM_lo_allweek[:,:,:,5]=np.copy(comp_ind5_lo)
    
    
    
    comp_ind5_hi=np.nanmean(comp_ind5_hi,2)
    comp_ind5_lo=np.nanmean(comp_ind5_lo,2)
    
    allcomphiJFM[:,:,5]=comp_ind5_hi
    allcomploJFM[:,:,5]=comp_ind5_lo


 
diffhi=allcomphiNDJ-allcomphiJFM
difflo=allcomploNDJ-allcomploJFM

del ind5_hi
del ind5_lo
#========PDO================================
obj=nc.Dataset('/Users/zhiqiyang/Downloads/era5_ar_gwv3_1940-2022.nc')

latpdo=obj.variables['lat'][120:281]  #60-20N
lonpdo=obj.variables['lon'][480:1041]  #120-260


#===========
name='PDO'
month='NDJ'
climatpdo_NDJ=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/ERA5all_num'+str(wind)+'week_indp_'+month+'_1940alldomain_bootstrap_nooverlap.npy')
climatpdo_NDJ=np.nanmean(climatpdo_NDJ,2)

ind5_hi=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/pdo1940-2018_'+month+'_ave1wk.txt')
ind5_lo=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/pdo1940-2018_'+month+'_ave1wk.txt')

ind5_hi[ind5_hi<0]=np.nan
ind5_lo[ind5_lo>0]=np.nan    

comp_ind5_hi=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/ERA5all_num'+str(wind)+'week_indp_'+month+'_1940alldomain_bootstrap_nooverlap.npy')
comp_ind5_lo=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/ERA5all_num'+str(wind)+'week_indp_'+month+'_1940alldomain_bootstrap_nooverlap.npy')
    
comp_ind5_hi[:,:,np.isnan(ind5_hi)]=np.nan
comp_ind5_lo[:,:,np.isnan(ind5_lo)]=np.nan

NDJ_hi_allpdo=np.copy(comp_ind5_hi)
NDJ_lo_allpdo=np.copy(comp_ind5_lo)


comp_ind5_hi=np.nanmean(comp_ind5_hi,2)
comp_ind5_lo=np.nanmean(comp_ind5_lo,2)

allcomphipdo_NDJ=comp_ind5_hi
allcomplopdo_NDJ=comp_ind5_lo

del ind5_hi
del ind5_lo
#----------
month='JFM'
climatpdo_JFM=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/ERA5all_num'+str(wind)+'week_indp_'+month+'_1940alldomain_bootstrap_nooverlap.npy')
climatpdo_JFM=np.nanmean(climatpdo_JFM,2)

ind5_hi=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/pdo1940-2018_'+month2+'_ave1wk.txt')
ind5_lo=np.loadtxt('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/pdo1940-2018_'+month2+'_ave1wk.txt')

ind5_hi[ind5_hi<0]=np.nan
ind5_lo[ind5_lo>0]=np.nan    

comp_ind5_hi=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/ERA5all_num'+str(wind)+'week_indp_'+month+'_1940alldomain_bootstrap_nooverlap.npy')
comp_ind5_lo=np.load('/Users/zhiqiyang/Documents/A-ucsd-postdoc/nooverlap/ERA5all_num'+str(wind)+'week_indp_'+month+'_1940alldomain_bootstrap_nooverlap.npy')
    
comp_ind5_hi[:,:,np.isnan(ind5_hi)]=np.nan
comp_ind5_lo[:,:,np.isnan(ind5_lo)]=np.nan

JFM_hi_allpdo=np.copy(comp_ind5_hi)
JFM_lo_allpdo=np.copy(comp_ind5_lo)

comp_ind5_hi=np.nanmean(comp_ind5_hi,2)
comp_ind5_lo=np.nanmean(comp_ind5_lo,2)

allcomphipdo_JFM=comp_ind5_hi
allcomplopdo_JFM=comp_ind5_lo



#========================================
    
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
        cs=m.contourf(x1,y1,allcomphipdo_NDJ-allcomphipdo_JFM-climatpdo_NDJ+climatpdo_JFM,cslev,cmap=cmap,extend='both')
        plt.title('Composite +'+namesss[nam_n]+'_'+str(wind)+'week ARs (NDJ-JFM)',fontsize=fs)
        
        for i in arange(0,np.shape(latpdo)[0],5):
            for j in arange(0,np.shape(lonpdo)[0],5):
                
                lin1=JFM_hi_allpdo[i,j,:]
                lin2=NDJ_hi_allpdo[i,j,:]
                
                lin1=lin1[~np.isnan(lin1)]
                lin2=lin2[~np.isnan(lin2)]
                
                repeated_lin1 = climatpdo_JFM[i, j]
                repeated_lin2 = climatpdo_NDJ[i, j]
                
                lin1=lin1-repeated_lin1
                lin2=lin2-repeated_lin2
                
                
                
                t_statistic, sig = stats.ttest_ind(lin1, lin2)
                
                if (sig < 0.1 ):
                    m.plot(x1[i,j],y1[i,j],'+',color='k',markersize=s,zorder=5)
                    del sig
                            
        parallels = np.arange(20.,60,10.)
        m.drawparallels(parallels,color=col,linewidth=linew,dashes=[1,0.1],labels=[1,0,0,0],fontsize=flalo,zorder=2)
        meridians = np.arange(120.,260.,20.)
        m.drawmeridians(meridians,color=col,linewidth=linew,dashes=[1,0.1],labels=[0,0,0,0],fontsize=flalo,zorder=2)
        
        ax.spines['bottom'].set_linewidth(linew)
        ax.spines['top'].set_linewidth(linew)
        ax.spines['left'].set_linewidth(linew)
        ax.spines['right'].set_linewidth(linew)
        
       
        
        
        ax = fig.add_axes([0.92,h2[nam_n-1],0.8,0.8])#  125 59
        m = Basemap(area_thresh=100000.,projection='cyl',llcrnrlat=20,urcrnrlat=56,llcrnrlon=133,urcrnrlon=255,resolution='l',ax=ax)
        m.drawcoastlines(linewidth=0.5)
        m.drawcountries(linewidth=0.5)
        m.drawstates(linewidth=0.5)
        
        x,y=np.meshgrid(lonpdo,latpdo)
        x1,y1=m(x,y) 
        #
        cs=m.contourf(x1,y1,allcomplopdo_NDJ-allcomplopdo_JFM-climatpdo_NDJ+climatpdo_JFM,cslev,cmap=cmap,extend='both')
        
        for i in arange(0,np.shape(latpdo)[0],5):
            for j in arange(0,np.shape(lonpdo)[0],5):
                
                lin1=JFM_lo_allpdo[i,j,:]
                lin2=NDJ_lo_allpdo[i,j,:]
                
                lin1=lin1[~np.isnan(lin1)]
                lin2=lin2[~np.isnan(lin2)]
                
                repeated_lin1 = climatpdo_JFM[i, j]
                repeated_lin2 = climatpdo_NDJ[i, j]
                
                lin1=lin1-repeated_lin1
                lin2=lin2-repeated_lin2
                
                t_statistic, sig2 = stats.ttest_ind(lin1, lin2)
                
                if (sig2 < 0.1 ):
                    m.plot(x1[i,j],y1[i,j],'+',color='k',markersize=s,zorder=5)
                    del sig2
        
       
                        
        plt.title('Composite -'+namesss[nam_n]+'_'+str(wind)+'week ARs (NDJ-JFM)',fontsize=fs)
        
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
        cs=m.contourf(x1,y1,diffhi[:,:,nam_n]-climatNDJ+climatJFM,cslev,cmap=cmap,extend='both')
        plt.title('Composite +'+namesss[nam_n]+'_'+str(wind)+'week ARs (NDJ-JFM)',fontsize=fs)
        
        for i in arange(0,np.shape(lat)[0],2):
            for j in arange(0,np.shape(lon)[0],2):
                
                lin1=JFM_hi_allweek[i,j,:,nam_n]
                lin2=NDJ_hi_allweek[i,j,:,nam_n]
                
                lin1=lin1[~np.isnan(lin1)]
                lin2=lin2[~np.isnan(lin2)]
                
                repeated_lin1 = climatJFM[i, j]
                repeated_lin2 = climatNDJ[i, j]
                
                lin1=lin1-repeated_lin1
                lin2=lin2-repeated_lin2
                
                
                
                t_statistic, sig = stats.ttest_ind(lin1, lin2)
                
                if (sig < 0.1 ):
                    m.plot(x1[i,j],y1[i,j],'+',color='k',markersize=s,zorder=5)
                    del sig
                            
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
        cs=m.contourf(x1,y1,difflo[:,:,nam_n]-climatNDJ+climatJFM,cslev,cmap=cmap,extend='both')
        
        for i in arange(0,np.shape(lat)[0],2):
            for j in arange(0,np.shape(lon)[0],2):
                
                lin1=JFM_lo_allweek[i,j,:,nam_n]
                lin2=NDJ_lo_allweek[i,j,:,nam_n]
                
                lin1=lin1[~np.isnan(lin1)]
                lin2=lin2[~np.isnan(lin2)]
                
                repeated_lin1 = climatJFM[i, j]
                repeated_lin2 = climatNDJ[i, j]
                
                lin1=lin1-repeated_lin1
                lin2=lin2-repeated_lin2
                
                t_statistic, sig2 = stats.ttest_ind(lin1, lin2)
                
                if (sig2 < 0.1 ):
                    m.plot(x1[i,j],y1[i,j],'+',color='k',markersize=s,zorder=5)
                    del sig2
        
       
                        
        plt.title('Composite -'+namesss[nam_n]+'_'+str(wind)+'week ARs (NDJ-JFM)',fontsize=fs)
        
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
        
        
        
savefig('Compo-window_1week_1982-2021indp_PDOdiff_NDJ-JFM_5ind_v3_nooverlap_response.pdf',bbox_inches = 'tight')
    
plt.show()
    
  
    
    
    
'''
for i in np.arange(0,1000):
    plt.contourf(allclimate[:,:,i])
    plt.show()

'''




