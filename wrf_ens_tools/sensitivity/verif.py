
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import os
from datetime import timedelta
from wrf_ens_tools.calc import FSS, ReliabilityTotal
import cartopy.crs as ccrs
import cartopy.feature as cfeat
# import cmocean
from wrf_ens_tools.sensitivity import Subset
from wrf_ens_tools.plots import nclcmaps


def storeEnsStats(ensprobpath, obpath, reliabilityobpath,
                  outpath, runinit, fhr,
                  probvar='updraft_helicity', probthresh=25., nbrhd=0.,
                  subset=False, rboxpath=None, **subsetparms):
    '''
    Stores ensemble verification stats using FSS and Reliability
    calculations from calc.py. Currently only supports updraft
    helicity and uh thresholds based on formatting described in
    probpath input.

    Inputs
    ------
    probpath ------ path to netCDF file containing
                    probability values. IMPORTANT
                    NOTE - prob file is expected to
                    be organized like in probcalcSUBSET.f
                     Variable P_HYD[0,i,:,:]:
                      i         Var
                     [0]   Refl > 40 dBZ probs
                     [1]   UH > 25 m2/s2 probs
                     [2]   UH > 40 m2/s2 probs
                     [3]   UH > 100 m2/s2 probs
                     [4]   Wind Speed > 40 mph probs
    obspath ------- path to netCDF file containing
                    observation verification values.
                    Currently, practically perfect
                    for verifying UH is the only
                    verification type supported.
    runinit ------- datetime obj describing model init
                    time.
    fhr ----------- integer specifying forecast hour to
                    store stats for.
    var ----------- string describing variable to verify.
                    Only supports 'updraft_helicity'
                    option as of right now.
    probthresh ---- float describing threshold of
                    variable to use when
                    pulling probs. Choices
                    for UH are 25, 40, and 100 m2/s2.
                    Choices for Reflectivity and
                    Win Speed are 40 (dbz) and
                    40 (mph) respectively.
    rboxpath ------ optional path to sensitivity calc
                    input file. Only needed if using
                    FSS to verify subsets. Otherwise
                    leave as None.
    nbrhd --------- optional float to store the neighborhood
                    distance (in km) of the probabilities.
    **subsetparms - accepts a variety of subset-related keyword
                    args if function is being used to verify
                    subsets. (See subset.py or Subset class
                    for more detail)
                     Options:
                     sens - sensitivity object
                     subset_size - size of desired subset
                     subset_method - accepts 'weight' or 'percent'
                     percent - sensitivity threshold float
                     sensvalfile - relative path to stored
                         sensitivity values
                     nbrhd - float specifying desired neighborhood
                         in km
                     thresh - response threshold (i.e. UH > 25 m2/s2 probs,
                         25 is the threshold)
                     sensvars - list of sensitivity variables to use for
                         subset. See subset.py or esens_subsetting.py for
                         info on naming conventions.
                     analysis_type - supports "RAP" or "WRF" at the moment
    '''
    if subset:
        sub = Subset(**subsetparms)
        sub.storeUHStats(outpath, obpath, reliabilityobpath)
    else:
        if os.path.exists(outpath):
            statsout = Dataset(outpath, 'a')
        else:
            statsout = Dataset(outpath, 'w')
            statsout.createDimension('Times',  None)
            statsout.createDimension('vars', None)
            statsout.createDimension('rbox', 4)
            statsout.createDimension('rel', 3)
            statsout.createDimension('bins', 10)
            statsout.createVariable('Run_Init', str, ('Times'))
            statsout.createVariable('Neighborhood', float, ('Times'))
            statsout.createVariable('Full_Ens_Reliability_Total',
                                    float, ('Times', 'rel', 'bins'))
            statsout.createVariable('Full_Ens_FSS_Total', float, ('Times'))

            ###########################################################
            # In non-subsetting terms, the response_func simply
            #   refers to the variable for which the probs are
            #   valid. Response_thresh specifies what threshold
            #   the probs are for (e.g. UH > 25 m2/s2 - 25 is the
            #   threshold). Response_time describes the forecast
            #   hour for which the probabilities are valid. For
            #   instance, a response time of 23 means the 1-hr probs
            #   are valid from f22-f23.
            ###########################################################
            statsout.createVariable('Response_Func', str, ('Times'))
            statsout.createVariable('Response_Thresh', float, ('Times'))
            statsout.createVariable('Response_Time', int, ('Times'))
            ###########################################################
            # The practically perfect sigma describes the sigma used
            #    in the gaussian filter for creating "probabilistic"
            #    storm reports. The sigma value used operationally
            #    to verify operational SPC convective outlooks is 2.
            #    See Sobash et al 2011 for more info.
            ###########################################################
            statsout.createVariable('Prac_Perf_Sigma', int, ('Times'))
        #### Continue onto calculations ###
        if os.path.exists(obpath):
            rtimedate = runinit + timedelta(hours=fhr)
            ## Start with FSS ##
            fens_fss = FSS(ensprobpath, obpath, rtimedate, var='updraft_helicity',
                          thresh=probthresh, rboxpath=rboxpath)
            ## Reliability next... ##
            fens_reliability = ReliabilityTotal(ensprobpath, runinit, fhr,
                                            obpath=None, var='updraft_helicity',
                                            thresh=probthresh, rboxpath=rboxpath,
                                            sixhr=False)
            # If not subsetting, the reliability values for the response box will all be missing vals
            prob_bins, f_fcstfreq_tot, ob_hr_tot = fens_reliability
            f_fss_tot, f_fss_rbox, sig = fens_fss
            # Start pulling variables
            init = statsout.variables['Run_Init']
            rfunc = statsout.variables['Response_Func']
            resptime = statsout.variables['Response_Time']
            fensfsstot = statsout.variables['Full_Ens_FSS_Total']
            nbr = statsout.variables['Neighborhood']
            fens_rel_tot = statsout.variables['Full_Ens_Reliability_Total']
            respthresh = statsout.variables['Response_Thresh']
            pperfsig = statsout.variables['Prac_Perf_Sigma']
            n = len(init)
            init[n] = runinit
            rfunc[n] = probvar
            resptime[n] = fhr
            respthresh[n] = probthresh
            fensfsstot[n] = f_fss_tot
            pperfsig[n] = sig
            nbr[n] = nbrhd
            fens_rel_tot[n,:,:] = np.atleast_2d(np.vstack((prob_bins,
                                                f_fcstfreq_tot, ob_hr_tot)))[:]
            statsout.close()
            return

# Plot FSS against a variable from the netCDF file
def plotUHFSS(ncfilepaths, xvarname,
              sigmas=[0,1,2], uhthresholds=[25., 40., 100.], nbrs=[30, 45, 60],
              senstimes=[6,12], normalize=False,
              onlyplot=None, subset=False):
    '''
    Plots FSS from a netCDF file containing output either from the
    storeEnsStats() function or the storeUHStats() from the subset
    library.
    '''
    for sigma in sigmas:
        for uhthresh in uhthresholds:
            for nbr in nbrs:
                for stime in senstimes:
                    figpath = "fss_by_{}_sig{}_uh{}_nbr{}_stimef{}.png".format(xvarname,
                                      sigma, uhthresh, nbr, stime)
                    plt.figure(1, figsize=(12,10))
                    ncruns = []
                    for path in ncfilepaths:
                        print("NC File Path {}; Sigma {}; Neighborhood {}; UH Thresh {}; Sens Time: {}".format(path,
                              sigma, nbr, uhthresh, stime))
                        # Open dataset
                        fssdat = Dataset(path)
                        # Pull all values of sigma, response thresholds,
                        #  sens times, and nbrhds
                        ncsig = fssdat.variables['Prac_Perf_Sigma'][:]
                        ncprobthresh = fssdat.variables['Response_Thresh'][:]
                        ncnbrhd = fssdat.variables['Neighborhood'][:]
                        ncsenstimes = fssdat.variables['Sens_Time'][:]
                        # Find indices for particular combo of sigma,
                        #  threshold, sens times, and neighborhood values
                        inds = np.where((ncsig == sigma) & (ncprobthresh == uhthresh) & (ncnbrhd == nbr) & (ncsenstimes == stime))
                        # Run initialization
                        ncrun = fssdat.variables['Run_Init'][inds]
                        ncruns.append(ncrun)
                        # Pull FSS values
                        f_tot_fss = fssdat.variables['Full_Ens_FSS_Total'][inds]
                        # Define x axis increments with xvarname
                        xvar = fssdat.variables[xvarname][inds]
                        x = np.unique(xvar)
                        # Pull many more extra variables if plotting subsets
                        if subset:
                            f_rbox_fss = fssdat.variables['Full_Ens_FSS_Rbox'][inds]
                            s_tot_fss = fssdat.variables['Subset_FSS_Total'][inds]
                            s_rbox_fss = fssdat.variables['Subset_FSS_Rbox'][inds]
                            # Pull more subset specific vars
                            ncsensvars = fssdat.variables['Sens_Vars'][inds]
                            ncanalysis = fssdat.variables['Analysis'][inds]
                            ncmethod = fssdat.variables['Subset_Method'][inds]
                            if xvarname != 'Sens_Threshold':
                                # If xvar isn't sens_threshold, it'll be subset size
                                ncyvar = fssdat.variables['Sens_Threshold'][inds]
                            else:
                                # Will need subset size too
                                ncyvar = fssdat.variables['Subset_Size'][inds]
                            # List of unique sensitivity variable lists
                            sensvars =  [list(x) for x in set(tuple(x) for x in ncsensvars)]
                            nvars = len(sensvars)
                            # Get color lists based on number of different sensitivity variables
                            allcmap = nclcmaps.cmap('grads_rainbow').colors
                            rboxcmap = nclcmaps.cmap('grads_rainbow').colors
                            cmapincr = np.linspace(0, len(allcmap)-1, nvars, dtype=int)
                            totcolors = [allcmap[x] for x in cmapincr]
                            rboxcolors = [rboxcmap[x] for x in cmapincr]
                            i = 0
                            for sensvar in sensvars:
                                if ('500_hPa_GPH' in sensvar) and ('SLP' in sensvar):
                                    varlabel = 'All'
                                else:
                                    varlabel = ' '.join(var for var in sensvar)
                                # For legend label purposes
                                plt.plot(0, 0, color=totcolors[i],
                                         label="FSS Total Domain with Sensitivity Variables: {}".format(varlabel))
                                plt.plot(0, 0, color=rboxcolors[i],
                                        label="FSS Response Box with Sensitivity Variables: {}".format(varlabel))
                                # Need to create mask here because shape of sensvars is different
                                #   from shape of other vars (can't broadcast together to create mask)
                                sens_mask = []
                                for k in range(len(ncsensvars[:,0])):
                                    if list(ncsensvars[k]) == list(sensvar):
                                        sens_mask.append(True)
                                    else:
                                        sens_mask.append(False)
                                # Move onto rest of variables
                                for method in np.unique(ncmethod):
                                    for var in np.unique(ncyvar):
                                        ###### Set colors here ######
                                        if xvarname == 'Subset_Size':
                                            # Means we can plot non-sensitivity-based subsets by subset size
                                            if (method == "percent") and (var == 0.):
                                                totcolor = 'grey'
                                                rboxcolor = 'black'
                                                lw = 3
                                                a = 1.
                                                totcolor = totcolors[i]
                                                rboxcolor = rboxcolors[i]
                                                lw = 2
                                                a = 0.5
                                                zorder=10
                                            else:
                                                totcolor = totcolors[i]
                                                rboxcolor = rboxcolors[i]
                                                lw = 2
                                                a = 0.5
                                                zorder=3
                                        else:
                                            totcolor = totcolors[i]
                                            rboxcolor = rboxcolors[i]
                                            lw=2
                                            a=0.5
#                                            if (method == 'percent'):
#                                                lw = 2
#                                                a = 0.6
#                                            else:
#                                                lw=2
#                                                a=0.0
                                        ### Done with color setting #####
                                        for analysis in np.unique(ncanalysis):
                                            xtmp = xvar[sens_mask]
                                            sub_fss_all_tmp = s_tot_fss[sens_mask]
                                            sub_fss_rbox_tmp = s_rbox_fss[sens_mask]
                                            subinds = np.where((ncmethod[sens_mask] == method) & (ncyvar[sens_mask] == var) & (ncanalysis[sens_mask] == analysis))
                                            if normalize:
                                                print('Normalizing by full ensemble, so FSS vals are full-ens-relative')
                                                #y_tot_domain = sub_fss_all_tmp[subinds] - f_tot_fss[0]
                                                y_rbox = sub_fss_rbox_tmp[subinds] - f_rbox_fss[0]
                                            else:
                                                y_tot_domain = sub_fss_all_tmp[subinds]
                                                y_rbox = sub_fss_rbox_tmp[subinds]
                                            # Careful with this if-statement
                                            # Only use when you know zeroes are
                                            #  due to probcalc.f bug and not
                                            #  actual zero FSS values. Otherwise
                                            #  leave commented out.
                                            if (uhthresh < 100) & (xvarname == 'Subset_Size'):
                                                totmask = (xtmp[subinds] > 4) & (y_tot_domain < 0.001)
                                                rboxmask = (xtmp[subinds] > 4) & (y_rbox < 0.001)
                                                y_tot_domain = np.ma.masked_array(y_tot_domain, mask=totmask)
                                                y_rbox = np.ma.masked_array(y_rbox, mask=rboxmask)
                                            elif (xvarname == 'Subset_Size'):
                                                totmask = (xtmp[subinds] > 6) & (y_tot_domain < 0.001)
                                                rboxmask = (xtmp[subinds] > 6) & (y_rbox < 0.001)
                                                y_tot_domain = np.ma.masked_array(y_tot_domain, mask=totmask)
                                                y_rbox = np.ma.masked_array(y_rbox, mask=rboxmask)
                                            # Plot single line representing subset of a unique parameter set by subset size
                                            plt.plot(xtmp[subinds], y_tot_domain, color=totcolor, lw=lw, alpha=a,zorder=zorder)
                                            plt.plot(xtmp[subinds], y_rbox, color=rboxcolor, lw=lw, alpha=a,zorder=zorder)
                                i += 1
                            fssdat.close()
                        plt.grid()
                        print(np.unique(f_rbox_fss))
                        if normalize == False:
                            plt.plot(x, np.ones_like(x)*np.unique(f_rbox_fss)[0], color='maroon',
                                    linewidth=4, linestyle='--', label="Full Ensemble FSS Response Box", zorder=13)
                            plt.plot(x, np.ones_like(x)*np.unique(f_tot_fss)[0], color='midnightblue',
                                linewidth=4, linestyle='--', label="Full Ensemble FSS Total Domain", zorder=13)
                            plt.ylim(0, 1)
                    plt.xlabel(xvarname.replace('_',' '), fontsize=12)
                    if normalize:
                        plt.ylabel('Fractional Skill Score Difference (Subset FSS - Full Ensemble FSS', fontsize=12)
                        plt.title(r"FSS for Runs Initiated {}".format(', '.join(str(run) for run in np.unique(ncruns[0][:]))) + '\n' + r"Practically Perfect $\sigma$ = {}; Neighborhood = {} km; UH Threshold = {} m$^2$/s$^2$; Sens Time = f{}".format(sigma,
                                  nbr, uhthresh, stime), fontsize=14)
                    else:
                        plt.ylabel('Fractional Skill Score', fontsize=12)
                        plt.title(r"FSS for Runs Initiated {}".format(', '.join(str(run) for run in np.unique(ncruns[0][:]))) + '\n' + r"Practically Perfect $\sigma$ = {}; Neighborhood = {} km; UH Threshold = {} m$^2$/s$^2$; Sens Time = f{}".format(sigma,
                                  nbr, uhthresh, stime), fontsize=14)
                    plt.xlim(x.min(), x.max())
                    l = plt.legend(fontsize=10, loc=9, bbox_to_anchor=(0.5,0.), borderaxespad=4.)
                    l.set_zorder = 12
                    plt.savefig(figpath, bbox_extra_artists=(l,), bbox_inches='tight')
                    plt.close()
    return

def storeUHStatsCSV(outpath, probpath, pperfpath, reliabilityobpath,
                    verif_fhr, nbrhd, probthresh, probvar='P_HYD'):
    """
    Stores UH verification stats to csv file.

    Inputs
    ------
    outpath ----------- absolute output filepath as a string
    pperfpath --------- absolute filepath for practically perfect
                        gridded data
    reliabilityobpath - absolute filepath for reliability observational
                        point data on ensemble grid
    verif_fhr -------- forecast hour to verify (e.g. 24 to verify
                        UH between forecast hours 23 and 24)
    nbrhd ------------- neighborhood distance in km
    probthresh -------- probability threshold (e.g. UH > 25 m2/s2)
    probvar ----------- name of netCDF variable in which
                        probabilities are stored (as string).
                        If using probability calculations from
                        this library, probabilities will by default
                        be stored in 'P_HYD'

    Outputs
    -------
    returns NULL and stores UH verification stats to outpath as csv
    """
    probdat = Dataset(probpath)
    ######## THIS WILL NOT WORK ON THE FIRST TRY #######################
    ######## Get correct format string when quanah returns #############
    runinit = datetime(probdat.START_DATE)
    veriftime = runinit + timedelta(hours=verif_fhrs)
    var = 'updraft_helicity'
    if os.path.exists(pperfpath):
        fens_fss = FSS(probpath, pperfpath, veriftime, var=probvar,
                  thresh=probthresh, rboxpath=S.getDir()+'esens.in')
        fens_reliability = ReliabilityTotal(probpath, runinit, verif_fhr,
                                       obpath=reliabilityobpath, var='updraft_helicity',
                                       thresh=probthresh, rboxpath=None,
                                       sixhr=False, nbrhd=nbrhd)
        # rbox parameters will be zeros since we didn't provide any response box
        #   information.
        prob_bins, f_fcstfreq_tot, f_ob_hr_tot = fens_reliability
        f_fss_tot, f_fss_rbox, sig = fens_fss
    else:
        raise FileNotFoundError('Please run storePracPerf() or set correct obpath.')
    if os.path.exists(outpath):
        mode = 'a'
    else:
        mode = 'w'
    # Open file with specified mode
    with open(outpath, mode=mode) as outfile:
        if mode == 'w':
            # If new file, add header row
            cols = ['Run_Init', 'Response_Func',
                    'Response_Time', 'Full_Ens_FSS_Total',
                    'Prac_Perf_Sigma', 'Neighborhood',
                    'Response_Thresh', 'Full_Ens_Reliability_Total']
            outfile_writer = csv.writer(outfile, delimiter=',')
            outfile_writer.writerow(cols)
        # Add row valid for current ensemble
        entry = [runinit, var, verif_fhr, f_fss_tot, sig,
                self._nbr, probthresh, (prob_bins, f_fcstfreq_tot, f_ob_hr_tot)]
        outfile_writer = csv.writer(outfile, delimiter=',')
        outfile_writer.writerow(entry)

def plotHistUHFSS(statspaths, outpath, xvarname='Subset_Size'):
    '''
    Plots 3D histogram of FSS by some chosen variable.
    '''
    pass

def plotReliabilityfromCSV(statspaths, outpath, subset=False, subgroupby="Sens_Vars",
                    sigmas=[0,1,2], uhthresholds=[25., 40., 100.], nbrs=[30, 45, 60],
                    senstimes=[6,12], subsize=10):
    '''
    Plot reliability diagram from csv files.
    '''
    # TO-DO: Add plotting functions from csv
    pass

def plotReliabilityfromNETCDF(statspaths, outpath, subset=False, subgroupby="Sens_Vars",
                    sigmas=[0,1,2], uhthresholds=[25., 40., 100.], nbrs=[30, 45, 60],
                    senstimes=[6,12], subsize=10):
    '''
    Plot reliability diagrams.
    '''
    for path in statspaths:
        stats = Dataset(path)
        probbins = []
        rel = []
        colors = []
        binfreq = []
        fensrel = stats.variables['Full_Ens_Reliability_Total'][:]
        bins, fcstfreqtmp, obhitratetmp = fensrel[:,0], fensrel[:,1], fensrel[:,2]
        obhitrate = np.ma.masked_array(obhitratetmp, mask=9e9)
        fcstfreq = np.ma.masked_array(fcstfreqtmp, mask=9e9)
        # If subset
        if subset:
            ncsubsizes = stats.variables['Subset_Size'][:]
            ncsig = stats.variables['Prac_Perf_Sigma'][:]
            ncrthresh = stats.variables['Response_Thresh'][:]
            ncnbrhds = stats.variables['Neighborhood'][:]
            ncstimes = stats.variables['Sens_Time'][:]
            for sigma in sigmas:
                for uhthresh in uhthresholds:
                    for nbr in nbrs:
                        for stime in senstimes:
                            plt.figure(figsize=(10,10)) # Initialize figure
                            probbins = []
                            rel = []
                            colors = []
                            binfreq = []
                            # Get indices for subset of relibility array
                            inds = np.where((ncsubsizes == subsize) & \
                                            (ncsig == sigma) & (ncrthresh == uhthresh) & \
                                            (ncnbrhds == nbr) & (ncstimes == stime))
                            # Grab full ensemble reliability for given indices
                            # Full ensemble reliability for total
                            fensrel = stats.variables['Full_Ens_Reliability_Total'][inds]
                            bins, fcstfreqtmp, obhitratetmp = fensrel[:,0], fensrel[:,1], fensrel[:,2]
                            obhitrate = np.ma.masked_array(obhitratetmp, mask=(obhitratetmp==9e9))
                            fcstfreq = np.ma.masked_array(fcstfreqtmp, mask=(obhitratetmp==9e9))
                            print(bins[0], obhitrate[0])
                            plt.plot(bins[0], obhitrate[0], color='navy', linestyle='-.',
                                     linewidth=3, label='Full Ens Total Domain Reliability', zorder=3)
                            # Full ensemble reliability for response box
                            fensrelrbox = stats.variables['Full_Ens_Reliability_Rbox'][inds]
                            bins, fcstfreq, obhitrate = fensrelrbox[:,0], fensrelrbox[:,1], fensrelrbox[:,2]
                            obhitrate = np.ma.masked_array(obhitratetmp, mask=(obhitratetmp==9e9))
                            fcstfreq = np.ma.masked_array(fcstfreqtmp, mask=(fcstfreqtmp==9e9))
                            print(bins[0], obhitrate[0])
                            plt.plot(bins[0], obhitrate[0], linewidth=3,
                                     linestyle='-.', color='maroon', zorder=4,
                                     label='Full Ens Response Box Reliability')
                            # Going to average reliability by this variable
                            groupby = stats.variables[subgroupby][inds]
                            # Finding the set of unique values in most
                            # variables is easy, but when working with
                            # lists of lists associated with sensitivity
                            # variabes, we need to do a little extra work.
                            if subgroupby == 'Sens_Vars':
                                unique =  [list(x) for x in set(tuple(x) for x in groupby)]
                            else:
                                unique = np.unique(groupby)
                            nvars = len(unique)
                            # Get color lists based on number of different sensitivity variables
                            allcmap = nclcmaps.cmap('grads_rainbow').colors
                            rboxcmap = nclcmaps.cmap('grads_rainbow').colors
                            cmapincr = np.linspace(0, len(allcmap)-1, nvars, dtype=int)
                            totcolors = [allcmap[x] for x in cmapincr]
                            rboxcolors = [rboxcmap[x] for x in cmapincr]
                            i = 0
                            # Average reliability for each unique groupby variable value
                            for groupbyval in unique:
                                if subgroupby == 'Sens_Vars':
                                    # Format label for legend
                                    if ('500_hPa_GPH' in groupbyval) and ('SLP' in groupbyval):
                                        varlabel = 'All'
                                    else:
                                        varlabel = ' '.join(var for var in groupbyval)
                                    # Iterate through lists of lists
                                    mask = []
                                    for k in range(len(groupby[:,0])):
                                        if list(groupby[k]) == list(groupbyval):
                                            mask.append(True)
                                        else:
                                            mask.append(False)
                                    plt.plot(10, 0, linewidth=2, color=totcolors[i],
                                             linestyle='--',
                                             label='Rel for Total Domain Sens Vars: {}'.format(varlabel))
                                    plt.plot(10, 0, linewidth=2, color=rboxcolors[i],
                                             label='Rel for Response Box Sens Vars: {}'.format(varlabel))
                                else:
                                    mask = (groupby == groupbyval)
                                    varlabel = str(groupbyval)
                                # Subset reliability for total domain
                                subrel = stats.variables['Subset_Reliability_Total'][inds][mask]
                                bins, fcstfreqtmp, obhitratetmp = subrel[:,0], subrel[:,1], subrel[:,2]
                                probbins.append(bins[0])
                                obhitrate = np.ma.masked_array(obhitratetmp, mask=(obhitratetmp==9e9))
                                fcstfreq = np.ma.masked_array(fcstfreqtmp, mask=(fcstfreqtmp==9e9))
                                rel.append(np.mean(obhitrate, axis=0))
                                binfreq.append(np.mean(fcstfreq, axis=0))
                                colors.append(totcolors[i])
                                # Subset reliability for response box
                                subrelrbox = stats.variables['Subset_Reliability_Rbox'][inds][mask]
                                bins, fcstfreq, obhitrate = subrelrbox[:,0], subrelrbox[:,1], subrelrbox[:,2]
                                probbins.append(bins[0])
                                obhitrate = np.ma.masked_array(obhitratetmp, mask=(obhitratetmp==9e9))
                                fcstfreq = np.ma.masked_array(fcstfreqtmp, mask=(fcstfreqtmp==9e9))
                                rel.append(np.mean(obhitrate, axis=0))
                                binfreq.append(np.mean(fcstfreq, axis=0))
                                colors.append(rboxcolors[i])
                                i+=1
                        probbins, rel, binfreq = np.array(probbins[:]), np.array(rel[:]), np.array(binfreq)
                        mask = (rel == 9e9)
                        rel = np.ma.masked_array(rel, mask=mask)
                        for i in range(len(rel)-1):
                            plt.plot(probbins[i], rel[i], color=colors[i], linewidth=2, alpha=0.8)
                        plt.xticks(probbins[0], ['0-10', '10-20', '20-30', '30-40', '40-50',
                                   '50-60', '60-70', '70-80', '80-90', '90-100'])
                        l = plt.legend(fontsize=10, loc=9, bbox_to_anchor=(0.5,0.), borderaxespad=4.)
                        l.set_zorder(12)
                        figpath = outpath + '/rel_rthresh{}_subsize{}_sig{}_nbr{}_stime{}.png'.format(uhthresh,
                                                         subsize, sigma, nbr, stime)
                        plt.title('Reliability Diagram of Subsets and Full Ensemble \n\
                                  for Sens Time: {}, Subset Size: {}, Response Threshold: {}, Neighborhood: {}'.format(stime, subsize, uhthresh, nbr))
                        plt.tight_layout()
                        plt.savefig(figpath)
                        plt.close()
        else:
            probbins.append(bins)
            rel.append(obhitrate)
            colors.append('blue')
            binfreq.append(fcstfreq)
            plt.figure(figsize=(10,10))
            for i in range(len(rel)):
                for j in range(len(rel[0,:,0])):
                    plt.plot(probbins[0,0,:], rel[i,j,:], color=colors[i], alpha=0.5)
            plt.legend()
            plt.savefig(outpath)
            plt.close()
        return

def plotPracPerf(pperfpath, sixhour=False,
                 wrfoutref="/lustre/research/bancell/aucolema/HWT2016runs/2016050800/wrfoutREFd2"):
    """
    Plots practically perfect probabilities given
    the desired practically perfect filepath and
    store it to the designated outfile path.
    """
    pperfnc = Dataset(pperfpath)
    init = pperfnc.START_DATE
    initdate = datetime(int(init[:4]), int(init[5:7]),
                        int(init[8:10]), int(init[11:13]))
    fhrs = pperfnc.variables['fhr'][:]
    sigma = pperfnc.variables['sigma'][:]
    pperf = pperfnc.variables['practically_perfect']
    ref = Dataset(wrfoutref)
    lons = ref.variables['XLONG'][0]
    lats = ref.variables['XLAT'][0]
    print(np.shape(pperf))
    print(np.shape(lons))

    for t in range(len(fhrs)):
        fig = plt.figure(figsize=(10, 10))
        time = initdate + timedelta(hours=int(fhrs[t]))
        # Build map
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.LambertConformal())
        ax.set_extent([-108., -85., 28., 45.])
        state_borders = cfeat.NaturalEarthFeature(category='cultural',
                   name='admin_1_states_provinces_lakes', scale='50m', facecolor='None')
        ax.add_feature(state_borders, linestyle="-", edgecolor='dimgray')
        ax.add_feature(cfeat.BORDERS, edgecolor='dimgray')
        ax.add_feature(cfeat.COASTLINE, edgecolor='dimgray')
        cflevels = np.arange(10, 110, 10)
        cf = ax.contourf(lons, lats, pperf[t], cflevels, transform=ccrs.PlateCarree(),
                         cmap=nclcmaps.cmap("perc2_9lev"))
        #plt.clabel(cs)
        plt.colorbar(cf, ax=ax, fraction=0.039, pad=0.01, orientation='vertical', label='Probability')
        if sixhour:
            plt.title(r'$\sigma = $' + str(sigma[0]) + ' Practically Perfect valid: ' + \
                      str(time-timedelta(hours=6)) + ' to ' + str(time))
            plt.savefig('sixhr_pperf_f{}.png'.format(fhrs[t]))
        else:
            plt.title(r'$\sigma = $' + str(sigma[0]) + ' Practically Perfect valid: ' + \
                      str(time-timedelta(hours=1)) + ' to ' + str(time))
            plt.savefig('onehr_pperf_f{}.png'.format(fhrs[t]))
