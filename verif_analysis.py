"""
Created on Fri Aug 10 15:54:23 2018

@author: aucolema
"""

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from os import chdir
from datetime import datetime, timedelta

ncfiles = ['/lustre/research/bancell/aucolema/HWT2016runs/2016050900/stats_reduced.nc']

# Plot FSS for subset size/pperf sigma
fullensnum = 42
sigma = 1
stime = 6
figpath = 'fss_scatter_sigma{}_stime{}.png'.format(sigma, stime)
a = 0.3
runs = []
plt.figure(1, figsize=(10,10))
n = 0
for file in ncfiles:
    dat = Dataset(file)
    sig = dat.variables['Prac_Perf_Sigma'][:]
    senstime = dat.variables['Sens_Time'][:]
    inds = np.where((sig == sigma) & (senstime == stime))
    rtime = dat.variables['Response_Time'][inds]
    subsize = dat.variables['Subset_Size'][inds]       
    fss_all = dat.variables['Subset_FSS_Total'][inds]
    fss_rbox = dat.variables['Subset_FSS_Rbox'][inds]
    fullens_fss_all = dat.variables['Full_Ens_FSS_Total'][inds]
    fullens_fss_rbox = dat.variables['Full_Ens_FSS_Rbox'][inds]
    init = dat.variables['Run_Init'][inds]
    runs.append(str(init[0]))
    sub_full_mean = []
    sub_rbox_mean = []
    for size in np.unique(subsize[subsize<fullensnum]):
        sub_inds = np.where(subsize == size)
        sub_full_mean.append(np.mean(fss_all[sub_inds]))
        sub_rbox_mean.append(np.mean(fss_rbox[sub_inds]))
    sub_full_mean, sub_rbox_mean = np.array(sub_full_mean), np.array(sub_rbox_mean)
    if n > 0:
        totlabel = None
        rboxlabel = None
    else:
        totlabel = 'Total Domain FSS'
        rboxlabel = 'Response Box FSS'
    plt.scatter(subsize[subsize<fullensnum], fss_all[subsize<fullensnum], c='lightgreen', edgecolor='k', alpha=a,
            label=totlabel)
#    plt.scatter(np.unique(subsize[subsize<fullensnum]), sub_full_mean,
#            c='firebrick', edgecolor='k', s=60., alpha=1, label="Total Domain Mean FSS", zorder=13)
    plt.scatter(fullensnum*np.ones_like(subsize), fullens_fss_all, alpha=a, c='lightgreen', edgecolor='k')
    plt.scatter(subsize[subsize<fullensnum], np.array(fss_rbox[subsize<fullensnum]), c='purple', edgecolor='k', alpha=a,
            label=rboxlabel)
#    plt.scatter(np.unique(subsize[subsize<fullensnum]), sub_rbox_mean,
#            c='yellow', edgecolor='k', s=60., alpha=1, label="Response Box Mean FSS", zorder=13)

    plt.scatter(fullensnum*np.ones_like(subsize), np.array(fullens_fss_rbox), alpha=a, c='purple', edgecolor='k')
    n += 1
plt.legend()
print(runs)
plt.title(r'Fractional Skill Scores of Subsets from {} UTC run(s) with $\sigma$ = {}'.format(', '.join(str(i) for i in runs), 
          str(sigma)))
plt.xlabel('Subset Size')
plt.ylabel('FSS')
plt.grid()
plt.tight_layout()
plt.savefig(figpath)
plt.close()
dat.close()

# Plot FSS for each type of subset with subset size
uhthresholds = [40]
nbrhds = [30]
figpath = 'fss_subset_size_sigma{}_nbr{}_uh{}.png'.format(sigma, nbrhds[0], uhthresholds[0])
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
    fss_all = dat.variables['Subset_FSS_Total'][inds]
    fss_rbox = dat.variables['Subset_FSS_Rbox'][inds]
    sens_vars = [list(x) for x in set(tuple(x) for x in sensvars[:1])]
    colors_rbox = ['purple', 'darkorange', 'maroon', 'purple', 'midnightblue', 
                  'darkolivegreen', 'seagreen', 'red']
    colors_all = ['green', 'coral', 'indianred', 'mediumorchid', 'royalblue',
                   'olive', 'mediumseagreen', 'orangered']
    senstimes = dat.variables['Sens_Time'][inds]
    nbrs = dat.variables['Neighborhood'][inds]
    i = 0
    for sens_var in sens_vars:
        print(sens_var)
        color_all = colors_all[i]
        color_rbox = colors_rbox[i]
        lw = 1
        if ('500_hPa_GPH' in sens_var) and ('SLP' in sens_var):
            label_all = "FSS Tot for All Sens Vars"
            label_rbox = "FSS Rbox for All Sens Vars"
        else:
            label_all="FSS Tot for Sens Vars: {}".format(' '.join(str(var) for var in sens_var))
            label_rbox="FSS RBox for Sens Vars: {}".format(' '.join(str(var) for var in sens_var))
        plt.plot(0, 0, color=colors_all[i], label=label_all)
        plt.plot(0, 0, color=colors_rbox[i], label=label_rbox)
        sens_mask = []
        for k in range(len(sensvars)):
            if list(sensvars[k]) == list(sens_var):
                sens_mask.append(True)
            else:
                sens_mask.append(False)
        for method in np.unique(methods):
            for percent in np.unique(percents):
                for stime in np.unique(senstimes):
                    for nbr in nbrhds:
                        for uh_thresh in uhthresholds:
                            x = subsize[sens_mask]
                            #print(np.shape(fss_all))
                            fss_all_tmp = fss_all[sens_mask]
                            fss_rbox_tmp = fss_rbox[sens_mask]
                            sub_inds = np.where((methods[sens_mask] == method) & (percents[sens_mask] == percent) \
                                                & (senstimes[sens_mask] == stime) & (nbrs[sens_mask] == nbr) \
                                                & (uhthresh[sens_mask] == uh_thresh))
                            x = x[sub_inds]
                            bugmask = (subsize[sens_mask][sub_inds] > 4) & (fss_rbox_tmp[sub_inds] < 0.001)
                            y_all = fss_all_tmp[sub_inds]
                            y_rbox = fss_rbox_tmp[sub_inds]
                            x_masked = np.ma.masked_array(x, mask=bugmask)
                            y_all_masked = np.ma.masked_array(y_all, mask=bugmask)
                            y_rbox_masked = np.ma.masked_array(y_rbox, mask=bugmask)
                            if (method == 'percent') and (percent == 0.):
                                print('We have a winner!')
                                color_all = 'gray'
                                color_rbox = 'k'
                                lw = 3.
                                zorder = 10
                                alpha = 0.9
                            else:
                                color_all = colors_all[i]
                                color_rbox = colors_rbox[i]
                                zorder = 8
                                lw = 1.
                                alpha = 0.5
                            plt.plot(x_masked[x_masked<fullensnum], y_all_masked[x_masked<fullensnum], 'd', c=color_all, ls='-', 
                                     linewidth=lw, zorder=zorder, alpha=alpha)
                            plt.plot(x_masked[x_masked<fullensnum], y_rbox_masked[x_masked<fullensnum], 'd', c=color_rbox, ls='-', 
                                     linewidth=lw, zorder=zorder, alpha=alpha)
                            print(method, percent, stime, nbr, uhthresh[sub_inds], sig[sub_inds], y_all, y_rbox)
                            del sub_inds
                            del x
                            del y_all
                            del y_rbox
        i += 1 
    x = np.unique(subsize[subsize<fullensnum])         
    fens_all_masked = fens_all[(uhthresh == uhthresholds[0]) & (nbrs == nbrhds[0])]
    fens_rbox_masked = fens_rbox[(uhthresh == uhthresholds[0]) & (nbrs == nbrhds[0])]
    print(np.unique(fens_all_masked), np.unique(fens_rbox_masked))         
    plt.plot(x, np.ones_like(x)*np.unique(fens_all_masked)[0], ls='--',
                             linewidth=3., c='navy', zorder=10, label='Full Ensemble Total Domain FSS')
    plt.plot(x, np.ones_like(x)*np.unique(fens_rbox_masked)[0], ls='--',
                             linewidth=3., c='maroon', zorder=10, label='Full Ensemble Response Box FSS')
l = plt.legend(fontsize=5)
l.set_zorder(13)
if len(list(sens_var)) > 6:
    sensstr = 'All'
else:
    sensstr = ', '.join(str(var).replace('_', ' ') for var in list(sens_var))
plt.title(r'Fractional Skill Scores of Subsets from {} UTC run valid f{} with $\sigma$ = {}'.format(str(init[0]), 
          str(rtime[0]), str(sigma)) + '\n' + \
         r"UH Threshold = {} m$^2$/s$^2$; Neighborhood = {} km".format(uhthresholds[0], nbrhds[0]),fontsize=11)
plt.xlabel('Subset Size')
plt.ylabel('FSS')
plt.grid()
#plt.tight_layout()
plt.savefig(figpath)
dat.close()
        
             
