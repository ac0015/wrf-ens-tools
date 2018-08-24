"""
Created on Fri Aug 10 15:54:23 2018

@author: aucolema
"""

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from os import chdir
from datetime import datetime, timedelta

ncfiles = ['/lustre/research/bancell/aucolema/HWT2016runs/2016051300/fss.nc']
# For second plot
sens_var = ['300_hPa_GPH', '300_hPa_T', '300_hPa_U-Wind', '300_hPa_V-Wind',
             '500_hPa_GPH', '500_hPa_T', '500_hPa_U-Wind', '500_hPa_V-Wind',
             '700_hPa_GPH', '700_hPa_T', '700_hPa_U-Wind', '700_hPa_V-Wind',
             '850_hPa_GPH', '850_hPa_T', '850_hPa_U-Wind', '850_hPa_V-Wind',
             '925_hPa_T', '925_hPa_U-Wind', '925_hPa_V-Wind',
             'SLP', '2m_Temp', '2m_Dewpt', '10m_U-Wind', '10m_V-Wind']

# Plot FSS for subset size/pperf sigma
fullensnum = 42
sigma = 2
figpath = 'fss_scatter_sigma{}.png'.format(sigma)
a = 0.3
plt.figure(1, figsize=(10,10))
for file in ncfiles:
    dat = Dataset(file)
    sig = dat.variables['Prac_Perf_Sigma'][:]
    inds = np.where(sig == sigma)
    rtime = dat.variables['Response_Time'][inds]
    subsize = dat.variables['Subset_Size'][inds]       
    fss_all = dat.variables['Subset_FSS_Total'][inds]
    fss_rbox = dat.variables['Subset_FSS_Rbox'][inds]
    fullens_fss_all = dat.variables['Full_Ens_FSS_Total'][inds]
    fullens_fss_rbox = dat.variables['Full_Ens_FSS_Rbox'][inds]
    init = dat.variables['Run_Init'][inds]
    sub_full_mean = []
    sub_rbox_mean = []
    for size in np.unique(subsize[subsize<fullensnum]):
        sub_inds = np.where(subsize == size)
        sub_full_mean.append(np.mean(fss_all[sub_inds]))
        sub_rbox_mean.append(np.mean(fss_rbox[sub_inds]))
    sub_full_mean, sub_rbox_mean = np.array(sub_full_mean), np.array(sub_rbox_mean)
    plt.scatter(subsize[subsize<fullensnum], fss_all[subsize<fullensnum], c='lightgreen', edgecolor='k', alpha=a,
            label='Total Domain FSS')
    plt.scatter(np.unique(subsize[subsize<fullensnum]), sub_full_mean,
            c='firebrick', edgecolor='k', s=60., alpha=1, label="Total Domain Mean FSS")
    plt.scatter(fullensnum*np.ones_like(subsize), fullens_fss_all, alpha=a, c='lightgreen', edgecolor='k')
    plt.scatter(subsize[subsize<fullensnum], np.array(fss_rbox[subsize<fullensnum]), c='purple', edgecolor='k', alpha=a,
            label='Response Box FSS')
    plt.scatter(np.unique(subsize[subsize<fullensnum]), sub_rbox_mean,
            c='yellow', edgecolor='k', s=60., alpha=1, label="Response Box Mean FSS")

    plt.scatter(fullensnum*np.ones_like(subsize), np.array(fullens_fss_rbox), alpha=a, c='purple', edgecolor='k')
plt.legend()
plt.title(r'Fractional Skill Scores of Subsets from {} UTC run valid f{} with $\sigma$ = {}'.format(str(init[0]), str(rtime[0]), str(sigma)))
plt.xlabel('Subset Size')
plt.ylabel('FSS')
plt.grid()
plt.tight_layout()
plt.savefig(figpath)
plt.close()
dat.close()

# Plot FSS for each type of subset with subset size
if len(sens_var) > 6:
    figpath = 'fss_subset_size_sigma{}_{}.png'.format(sigma, 'all')
else:
    figpath = 'fss_subset_size_sigma{}_{}.png'.format(sigma, '-'.join(str(var) for var in sens_var))
plt.figure(2, figsize=(10,10))
for file in ncfiles:
    dat = Dataset(file)
    sig = dat.variables['Prac_Perf_Sigma'][:]
    inds = np.where(sig == sigma)
    del(sig)
    uhthresh = dat.variables['Response_Thresh'][inds]
    sig = dat.variables['Prac_Perf_Sigma'][inds]
    #print(np.min(uhthresh[inds]), np.max(uhthresh[inds]), np.min(sig[inds]), np.max(sig[inds]))
    #print(np.shape(inds))
    subsize = dat.variables['Subset_Size'][inds]
    methods = dat.variables['Subset_Method'][inds]
    percents = dat.variables['Sens_Threshold'][inds]
    sensvars = dat.variables['Sens_Vars'][inds]
    fens_all = dat.variables['Full_Ens_FSS_Total'][inds]
    fens_rbox = dat.variables['Full_Ens_FSS_Rbox'][inds]
    sens_mask = []
    for var in sensvars:
        #print(list(var), list(sens_var))
        if list(var) != list(sens_var):
           sens_mask.append(False)
        else:
           sens_mask.append(True)
    senstimes = dat.variables['Sens_Time'][inds]
    nbrs = dat.variables['Neighborhood'][inds]
    for method in np.unique(methods):
        for percent in np.unique(percents):
            for stime in np.unique(senstimes):
                for nbr in np.unique(nbrs):
                    for uh_thresh in np.unique(uhthresh):
                        x = subsize[sens_mask]
                        #print(np.shape(fss_all))
                        fss_all_tmp = fss_all[sens_mask]
                        fss_rbox_tmp = fss_rbox[sens_mask]
                        sub_inds = np.where((methods[sens_mask] == method) & (percents[sens_mask] == percent) \
                                            & (senstimes[sens_mask] == stime) & (nbrs[sens_mask] == nbr) \
                                            & (uhthresh[sens_mask] == uh_thresh))
                        x = x[sub_inds]
                        y_all = fss_all_tmp[sub_inds]
                        y_rbox = fss_rbox_tmp[sub_inds]
                        if (method == 'percent') and (percent == 0.):
                            print('We have a winner!')
                            color_all = 'darkgray'
                            color_rbox = 'grey'
                            lw = 1
                            zorder = 10
                        else:
                            color_all = 'lightgreen'
                            color_rbox = 'purple'
                            lw = 1
                            zorder = 8
                        plt.plot(x[x<fullensnum], y_all[x<fullensnum], 'd', c=color_all, ls='-', 
                                 linewidth=lw, zorder=zorder)
                        plt.plot(x[x<fullensnum], y_rbox[x<fullensnum], 'd', c=color_rbox, ls='-', 
                                 linewidth=lw, zorder=zorder)
                        print(method, percent, stime, nbr, uhthresh[sub_inds], sig[sub_inds], y_all, y_rbox)
                        del sub_inds
                        del x
                        del y_all
                        del y_rbox
                    
    plt.plot(np.unique(subsize[subsize<fullensnum]), 
             np.ones_like(np.unique(subsize[subsize<fullensnum]))*np.mean(fens_all), ls='--',
                             linewidth=2., c='maroon', label='Full Ensemble Mean Total Domain FSS')
    plt.plot(np.unique(subsize[subsize<fullensnum]),  
             np.ones_like(np.unique(subsize[subsize<fullensnum]))*np.mean(fens_rbox), ls='--',
                             linewidth=2., c='blue', label='Full Ensemble Mean Response Box FSS')
plt.legend()
if len(list(sens_var)) > 6:
    sensstr = 'All'
else:
    sensstr = ', '.join(str(var).replace('_', ' ') for var in list(sens_var))
plt.title(r'Fractional Skill Scores of Subsets from {} UTC run valid f{} with $\sigma$ = {}'.format(str(init[0]), 
          str(rtime[0]), str(sigma)) + \
         "\nSensitivity Variables: {}".format(sensstr),fontsize=11)
plt.xlabel('Subset Size')
plt.ylabel('FSS')
plt.grid()
#plt.tight_layout()
plt.savefig(figpath)
dat.close()
        
             
