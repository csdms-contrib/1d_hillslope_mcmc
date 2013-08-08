# -*- coding: utf-8 -*-
"""
Created on Mon Mar 08 13:50:11 2010

@author: Administrator
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from itertools import izip
argmax = lambda array: max(izip(array, xrange(len(array))))[1]

# parameters for histograms
n_hist_bins =100

# bounds of the probability windows
#CDF_interp_points = [0.025,0.05,0.25,0.75,0.95,0.975]
# these are points in the CDF that indicate 2sigma
# 1 sigma, and the mean value of the parameters
CDF_interp_points = [0.023,0.159,0.5,0.841,0.977]

# get the directory path
root=os.getcwd()

#load the chain and rate files
#paramData = np.loadtxt(root+'/test_chain2.txt', unpack=False)
#paramData2 = np.loadtxt(root+'/fixed_chain.chain', unpack=False)
paramData2 = np.loadtxt(root+'/chain2.chain', unpack=False)
dimData = paramData2.shape

n_params = 3


n_cols = dimData[1]
n_runs = dimData[0]
#print "n_cols: "+str(n_cols)+" and n times: "+str(n_runs)

paramData = paramData2[500:n_runs,:]
dimData = paramData.shape
n_cols = dimData[1]
n_runs = dimData[0]
zero_vec = np.zeros(n_runs)
#print "n_cols: "+str(n_cols)+" and n times: "+str(n_runs)

# now for each time slice create a normalized histogram
#weighted by the likliehoods
like_weights = paramData[:,8]
run = paramData[:,0]
tpeak = paramData[:,4]
Upeak = paramData[:,5]
Uwidth = paramData[:,6]


#print "like weights: "
#for i in range (0,len(like_weights)):
#    print str(like_weights[i])

#####
###
### All the plotting stuff below is used to just visually
### look for the burn in period
###
fig = plt.figure(1,figsize=(16,10))
ax = fig.add_subplot(311)
ax.fill_between(run, zero_vec, tpeak, 
                 where=None, alpha=0.1, facecolor = 'red',edgecolor='black',
                 linewidth=2)
title_str = 'T*_{peak} chain, n='+str(n_runs)
plt.title(title_str,size=30)

ax = fig.add_subplot(312)
ax.fill_between(run, zero_vec, Upeak, 
                 where=None, alpha=0.1, facecolor = 'green',edgecolor='black',
                 linewidth=2)
title_str = 'U*_{peak} chain, n='+str(n_runs)
plt.title(title_str,size=30)

ax = fig.add_subplot(313)
ax.fill_between(run, zero_vec, Uwidth, 
                 where=None, alpha=0.1, facecolor = 'blue',edgecolor='black',
                 linewidth=2)
title_str = 'U*_{width} chain, n='+str(n_runs)
plt.title(title_str,size=30)
plt.savefig('param_chains.png')
###
###
#####

# create a histogram from the data
# this data is normalized, but it is the
# probability *density* so the integral = 1
tPHist,tPbins=np.histogram(tpeak,bins=n_hist_bins,
	  range=(0.25,0.4),normed=True,weights=like_weights)
UPHist,UPbins=np.histogram(Upeak,bins=n_hist_bins,
	  range=(15,25),normed=True,weights=like_weights)
UwHist,Uwbins=np.histogram(Uwidth,bins=n_hist_bins,
	  range=(0.25,0.4),normed=True,weights=like_weights)	  

nbins = len(tPbins)
#print "bins "+str(len(tPbins))
#for i in range(0,len(tPbins)):
#    print str(tPbins[i])+" "+str(UPbins[i])+" "+str(Uwbins[i])

#print "histogram: "   
#for i in range(0,len(tPHist)):
#    print str(tPHist[i])+" "+str(UPHist[i])+" "+str(UwHist[i])
    
# convert into probability
tPt_density = sum(tPHist)
tP_prob = np.divide(tPHist,tPt_density)
UPt_density = sum(UPHist)
UP_prob = np.divide(UPHist,UPt_density)
Uwt_density = sum(UwHist)
Uw_prob = np.divide(UwHist,Uwt_density)

#print "probability"
#for i in range(0,len(tP_prob)):
#    print str(tP_prob[i])+ " "+ str(UP_prob[i])+" "+str(Uw_prob[i])


# now get the cumulative density
sz_DHist = tPHist.size
left_bar_tP = np.zeros(sz_DHist)
right_bar_tP = np.zeros(sz_DHist)
midpoint_bar_tP =  np.zeros(sz_DHist)
Data_cdf_tP = np.zeros(sz_DHist)
last_val_tP = 0
for i in range(0,sz_DHist):
    left_bar_tP[i]=tPbins[i]
    right_bar_tP[i]=tPbins[i+1]
    midpoint_bar_tP[i] = (left_bar_tP[i]+right_bar_tP[i])/2
    Data_cdf_tP[i]=tP_prob[i]+last_val_tP
    last_val_tP = Data_cdf_tP[i]
tP_bounds = np.interp(CDF_interp_points,Data_cdf_tP,midpoint_bar_tP)     

# now get the cumulative density
sz_DHist = UPHist.size
left_bar_UP = np.zeros(sz_DHist)
right_bar_UP = np.zeros(sz_DHist)
midpoint_bar_UP =  np.zeros(sz_DHist)
Data_cdf_UP = np.zeros(sz_DHist)
last_val_UP = 0
for i in range(0,sz_DHist):
    left_bar_UP[i]=UPbins[i]
    right_bar_UP[i]=UPbins[i+1]
    midpoint_bar_UP[i] = (left_bar_UP[i]+right_bar_UP[i])/2
    Data_cdf_UP[i]=UP_prob[i]+last_val_UP
    last_val_UP = Data_cdf_UP[i]
UP_bounds = np.interp(CDF_interp_points,Data_cdf_UP,midpoint_bar_UP)     
    
# now get the cumulative density
sz_DHist = UwHist.size
left_bar_Uw = np.zeros(sz_DHist)
right_bar_Uw = np.zeros(sz_DHist)
midpoint_bar_Uw =  np.zeros(sz_DHist)
Data_cdf_Uw = np.zeros(sz_DHist)
last_val_Uw = 0
for i in range(0,sz_DHist):
    left_bar_Uw[i]=Uwbins[i]
    right_bar_Uw[i]=Uwbins[i+1]
    midpoint_bar_Uw[i] = (left_bar_Uw[i]+right_bar_Uw[i])/2
    Data_cdf_Uw[i]=Uw_prob[i]+last_val_Uw
    last_val_Uw = Data_cdf_Uw[i]   
Uw_bounds = np.interp(CDF_interp_points,Data_cdf_Uw,midpoint_bar_Uw)   

n_bounds = len(Uw_bounds)
for i in range(0,n_bounds):
    print str(tP_bounds[i])+" "+str(UP_bounds[i])+" "+str(Uw_bounds[i])+" "
