#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 19:21:54 2018

@author: aucolema
"""

from datetime import datetime, timedelta
from wrf_ens_tools.sensitivity import Sens
from wrf_ens_tools.sensitivity import Subset
import os
import sys

################### EXPERIMENT SETUP ####################################
direc = str(sys.argv[1])
date = str(sys.argv[2])
sixhour=True
yr, mo, day, hr = int(date[:4]), int(date[4:6]), int(date[6:8]), int(date[8:])
print(yr, mo, day, hr)
init = datetime(yr, mo, day, hr)
sens_times = [6]
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
nbrs = [30.0]
pperfnbrs = [80.0]
rfunctions = ["UH Coverage"]
analysis = ["WRF", "RAP"]
analysis_fhr = 0
########################################################################

from wrf_ens_tools.post import post_process as post

# Move response box and time choices to correct directory
try:
    os.popen('mv /home/aucolema/subsetGUI.txt ' + direc + '.')
except:
    print('subsetGUI.txt not found. Assuming already in {}'.format(direc))

for rfunc in rfunctions:
    for senstime in sens_times:
        for thresh in uh_thresh:
            # Create Sens object
            S = Sens(infile=True, gui=True, run=init,
                     senstime=senstime, sixhr=sixhour,
                     rfuncstr=rfunc, rthresh=thresh,
                     ensbasepath=direc)
            print('Using sens obj:', str(S))
            #S.runAll()
            # Calculate Full Ensemble probs and Practically Perfect probs
            rtime = S.getRTime()
            rdate = init + timedelta(hours=rtime)
            pperfpaths = []
            # Calculate full ensemble probs and reliability obs
            # Using neighborhood for reliability obs
            # TO-DO:Convert reliability to allow for six hour time frames
            for nbr in nbrs:
                sub = Subset(S, nbrhd=nbr, thresh=thresh)
                sub.calcProbs(sub._fullens)
                outpath = direc + 'reliability_ob_nbr{}.nc'.format(int(nbr))
                print("Currently in:", os.getcwd())
                os.chdir(direc)
                print("Changed to:", os.getcwd())
                if os.path.exists(outpath):
                   os.remove(outpath)
                post.storeNearestNeighbor(init, [rtime],
                                          str(outpath), sixhour=sixhour,
                                          wrfrefpath="/lustre/scratch/aucolema/2016052600/wrfoutREFd2")
                print("YAY I'm done with reliability ob stuff")
            for pperfnbr in pperfnbrs:
                # Calculate practically perfect w different distance sigma vals
                if sixhour:
                    outpath = direc + 'sixhr_sigma_nbr{}_pperf.nc'.format(pperfnbr)
                else:
                    outpath = direc + 'onehr_sigma_nbr{}_pperf.nc'.format(pperfnbr)
                pperfpaths.append(outpath)
                if os.path.exists(outpath):
                   os.remove(outpath)
                post.storePracPerfSPCGrid(init, [rtime],
                                   str(outpath), nbrhd=pperfnbr,
                                   dx=80., # Using 80-km SPC grid for AMS results
                                   sixhour=sixhour,
                                   wrfrefpath='/lustre/scratch/aucolema/2016052600/wrfoutREFd2')
                                   #dx=sub.getHorizGridSpacingD2()/1000.,

            del sub
            # Move onto analysis interpolation (if needed)
            for anl in analysis:
                if anl == "WRF":
                    sensdate = S.getRunInit() + timedelta(hours=S.getSensTime())
                    TTUanalysisbasepath = "/lustre/research/bancell/SE2016/"
                    dirdate = sensdate.strftime('%Y%m%d%H') + '/'
                    anlfile = "anTTU{}.analysis".format(analysis_fhr)
                    fullpath = TTUanalysisbasepath + dirdate + anlfile
                    print("TTU Analysis directory:", fullpath)
                elif anl == "RAP":
                    fullpath = None
                sub = Subset(S, nbrhd=nbr, analysis_type=anl,
                                 wrfanalysis_to_post_path=fullpath)
                if os.path.exists(sub.getAnalysis()) == False:
                    if anl == "WRF":
                       sub.processWRFAnalysis(True)
                    elif anl == "RAP":
                        sub.interpRAP()
            del sub
            if sixhour:
                statspath = direc + 'brians_subsetting_way_sixhr_stats_sig1_nbr30.nc'
            else:
                statspath = direc + 'onehr_stats_sig1_nbr30.nc'
            # Get stats for different subset combos
            for anl in analysis:
                for subsize in subset_sizes:
                    for method in subset_methods:
                        for percent in percents:
                            for varlist in sensvars:
                                for nbr in nbrs:
                                    # Create Subset object
                                    sub = Subset(S, subset_size=subsize, subset_method=method,
                                                 percent=percent, sensvars=varlist,
                                                 nbrhd=nbr, thresh=thresh,
                                                 analysis_type=anl)
                                    print('Using subset obj:', str(sub))
                                    sub.calcSubset()
                                    sub.calcProbs(sub.getSubMembers())
                                    #sub.plotSixPanels()
                                    reliabilityobpath = direc + 'reliability_ob_nbr{}.nc'.format(int(nbr))
                                    for obpath in pperfpaths:
                                        sub.storeUHStatsNetCDF(outpath=statspath,
                                                            pperfpath=obpath,
                                                            reliabilityobpath=reliabilityobpath)
                                    del sub
