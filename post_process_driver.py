#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 19:21:54 2018

@author: aucolema
"""

from datetime import datetime, timedelta
from sens import Sens
from subset import Subset
import os
import sys

################### EXPERIMENT SETUP ####################################
direc = str(sys.argv[1])
date = str(sys.argv[2])
yr, mo, day, hr = int(date[:4]), int(date[4:6]), int(date[6:8]), int(date[8:])
print(yr, mo, day, hr)
init = datetime(yr, mo, day, hr)
sens_times = [6, 12]
subset_sizes = [1, 2, 4, 6, 10, 15, 20, 25, 30, 35, 40]
subset_methods = ['weight', 'percent']
# Each sig level and then all sens vars
sensvars = [['300_hPa_GPH', '300_hPa_T', '300_hPa_U-Wind', '300_hPa_V-Wind'], 
            ['500_hPa_GPH', '500_hPa_T', '500_hPa_U-Wind', '500_hPa_V-Wind'],
            ['700_hPa_GPH', '700_hPa_T', '700_hPa_U-Wind', '700_hPa_V-Wind'],
            ['850_hPa_GPH', '850_hPa_T', '850_hPa_U-Wind', '850_hPa_V-Wind'],
            ['925_hPa_T', '925_hPa_U-Wind', '925_hPa_V-Wind'],
            ['SLP', '2m_Temp', '2m_Dewpt', '10m_U-Wind', '10m_V-Wind'], 
            ['300_hPa_GPH', '300_hPa_T', '300_hPa_U-Wind', '300_hPa_V-Wind',
             '500_hPa_GPH', '500_hPa_T', '500_hPa_U-Wind', '500_hPa_V-Wind',
             '700_hPa_GPH', '700_hPa_T', '700_hPa_U-Wind', '700_hPa_V-Wind',
             '850_hPa_GPH', '850_hPa_T', '850_hPa_U-Wind', '850_hPa_V-Wind',
             '925_hPa_T', '925_hPa_U-Wind', '925_hPa_V-Wind',
             'SLP', '2m_Temp', '2m_Dewpt', '10m_U-Wind', '10m_V-Wind']]
percents = [0., 10., 20., 30., 40., 50., 60., 70., 80., 90.]
uh_thresh = [25., 40., 100.]
nbrs = [30., 45., 60.]
sigma = [0, 1, 2]
# TO-DO
analysis = ["RAP"]

########################################################################

import post_process as post

# Move response box and time choices to correct directory
try:
    os.popen('mv /home/aucolema/subsetGUI.txt ' + direc + '.')
except:
    print('subsetGUI.txt not found. Assuming already in {}'.format(direc))

for senstime in sens_times:    
    # Create Sens object    
    S = Sens(infile=True, gui=True, run=init, 
             senstime=senstime)
    print('Using sens obj:', str(S))
    S.runAll()
    # Calculate Full Ensemble probs and Practically Perfect probs
    rtime = S.getRTime()
    rdate = init + timedelta(hours=rtime)
    pperfpaths = []
    for sig in sigma:
        outpath = direc + 'sigma{}_pperf.nc'.format(sig)
        pperfpaths.append(outpath)
        post.storePracPerf(init, [rtime], 
                           outpath, sigma=sig)
    for nbr in nbrs:
        sub = Subset(S, nbrhd=nbr)
        sub.calcProbs(sub._fullens)
        outpath = direc + 'reliability_ob_nbr{}.nc'.format(int(nbr))
        post.storeNearestNeighbor(init, [rtime], 
                                  outpath, nbrhd=nbr)
    sub.interpRAP()
    del sub
    statspath = direc + 'stats.nc'
    # Get stats for different subset combos
    for subsize in subset_sizes:
        for method in subset_methods:
            for percent in percents:
                for varlist in sensvars:
                    for nbr in nbrs:
                        for thresh in uh_thresh:
                            # Create Subset object
                            sub = Subset(S, subset_size=subsize, subset_method=method,
                                         percent=percent, sensvars=varlist, nbrhd=nbr, thresh=thresh)
                            print('Using subset obj:', str(sub))
                            sub.calcSubset()
                            sub.calcProbs(sub.getSubMembers())
                            #sub.plotSixPanels()
                            reliabilityobpath = direc + 'reliability_ob_nbr{}.nc'.format(int(nbr))
                            for obpath in pperfpaths:
                                sub.storeUHStats(outpath=statspath,pperfpath=obpath, 
                                                 reliabilityobpath=reliabilityobpath)
                            del sub
                
                
                
                
                
                    
                                   
