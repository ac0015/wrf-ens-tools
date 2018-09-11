
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import os
from os import chdir
from datetime import datetime, timedelta
from calc import FSS, Reliability
from subset import Subset
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import nclcmaps

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
            statsout.createVariable('Full_Ens_Reliability_Total', float, ('Times', 'rel', 'bins'))
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
            fens_reliability = Reliability(ensprobpath, runinit, fhr,
                                            obpath=None, var='updraft_helicity',
                                            thresh=probthresh, rboxpath=rboxpath,
                                            sixhr=False)
            # If not subsetting, the reliability values for the response box will all be missing vals
            prob_bins, f_fcstfreq_tot, ob_hr_tot, f_fcstfreq_rbox, ob_hr_rbox = fens_reliability
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
            fens_rel_tot[n,:,:] = np.atleast_2d(np.vstack((prob_bins, f_fcstfreq_tot, ob_hr_tot)))[:]
            statsout.close()
            return

# Plot FSS against a variable from the netCDF file
def plotUHFSS(ncfilepaths, xvarname,
              sigmas=[0,1,2], uhthresholds=[25., 40., 100.], nbrs=[30, 45, 60], 
              onlyplot=None, subset=False):
    '''
    Plots FSS from a netCDF file containing output either from the
    storeEnsStats() function or the storeUHStats() from the subset
    library.
    '''
    for sigma in sigmas:
        for uhthresh in uhthresholds:
            for nbr in nbrs:
                figpath = "fss_by_{}_sig{}_uh{}_nbr{}.png".format(xvarname, sigma, uhthresh, nbr)
                plt.figure(1, figsize=(12,10))
                ncruns = []
                for path in ncfilepaths:
                    print("NC File Path {}; Sigma {}; Neighborhood {}; UH Thresh {}".format(path, sigma, nbr, uhthresh))
                    # Open dataset
                    fssdat = Dataset(path)
                    # Pull all values of sigma, response threshold, and nbrhds
                    ncsig = fssdat.variables['Prac_Perf_Sigma'][:]
                    ncprobthresh = fssdat.variables['Response_Thresh'][:]
                    ncnbrhd = fssdat.variables['Neighborhood'][:]
                    # Find indices for particular combo of sigma, threshold, and neighborhood values
                    inds = np.where((ncsig == sigma) & (ncprobthresh == uhthresh) & (ncnbrhd == nbr))
                    # Run init
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
                        ncsenstimes = fssdat.variables['Sens_Time'][inds]
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
                                            #totcolor = 'grey'
                                            #rboxcolor = 'black'
                                            totcolor = totcolors[i]
                                            rboxcolor = rboxcolors[i]
                                            lw = 2
                                            a = 0.5
                                        else:
                                            totcolor = totcolors[i]
                                            rboxcolor = rboxcolors[i]
                                            lw = 2
                                            a = 0.5
                                    else:
                                        totcolor = totcolors[i]
                                        rboxcolor = rboxcolors[i]
                                        lw = 2
                                        a = 0.5
                                    ### Done with color setting ######
                                    for stime in np.unique(ncsenstimes):
                                        for analysis in np.unique(ncanalysis):
                                            xtmp = xvar[sens_mask]
                                            sub_fss_all_tmp = s_tot_fss[sens_mask]
                                            sub_fss_rbox_tmp = s_rbox_fss[sens_mask]
                                            subinds = np.where((ncmethod[sens_mask] == method) & (ncyvar[sens_mask] == var) &                                                               (ncanalysis[sens_mask] == analysis) & (ncsenstimes[sens_mask] == stime))
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
                                            #print(method, var, analysis, stime, y_tot_domain, y_rbox)
                                            #print(xtmp[subinds])
                                            plt.plot(xtmp[subinds], y_tot_domain, color=totcolor, lw=lw, alpha=a)
                                            plt.plot(xtmp[subinds], y_rbox, color=rboxcolor, lw=lw, alpha=a)
                            i += 1
                        fssdat.close()
                        #print(np.unique(f_rbox_fss))
                        plt.grid()
                        plt.plot(x, np.ones_like(x)*np.unique(f_rbox_fss)[0], color='maroon', 
                                 linewidth=4, linestyle='--', label="Full Ensemble FSS Response Box")
                    #print(np.unique(f_tot_fss))
                    plt.plot(x, np.ones_like(x)*np.unique(f_tot_fss)[0], color='midnightblue',
                            linewidth=4, linestyle='--', label="Full Ensemble FSS Total Domain")
                plt.xlabel(xvarname.replace('_',' '))
                plt.ylabel('Fractional Skill Score')
                plt.ylim(0, 1)
                l = plt.legend(fontsize=10, loc=9, bbox_to_anchor=(0.5,0.), borderaxespad=4.)
                l.set_zorder = 12
                plt.title(r"FSS for Runs Initiated {}".format(', '.join(str(run) for run in np.unique(ncruns[0][:]))) + '\n' +                           r"Practically Perfect $\sigma$ = {}; Neighborhood = {} km; UH Threshold = {} m$^2$/s$^2$".format(sigma, nbr, uhthresh))
                plt.savefig(figpath, bbox_extra_artists=(l,), bbox_inches='tight')
                plt.close()
    return
    


# In[6]:


#plotUHFSS(ncfilepaths=ncfiles, xvarname='Subset_Size', subset=True)


# In[4]:


def plotReliability(statspaths, outpath, subset=False):
    '''
    Plot reliability diagrams.
    '''
    plt.figure()
    for path in statspaths:
        stats = Dataset(path)
        print(stats)
        probbins = []
        rel = []
        colors = []
        binfreq = []
        fensrel = stats.variables['Full_Ens_Reliability_Total'][:]
        bins, fcstfreq, obhitrate = fensrel[:,0], fensrel[:,1], fensrel[:,2]
        probbins.append(bins)
        rel.append(obhitrate)
        colors.append('blue')
        binfreq.append(fcstfreq)
        if subset:
            subrel = stats.variables['Subset_Reliability_Total'][:]
            bins, fcstfreqtmp, obhitratetmp = subrel[:,0], subrel[:,1], subrel[:,2]
            obhitrate = np.ma.masked_array(obhitratetmp, mask=9e9)
            fcstfreq = np.ma.masked_array(fcstfreqtmp, mask=9e9)
            probbins.append(bins)
            rel.append(obhitrate)
            colors.append('orange')
            binfreq.append(fcstfreq)
            subrelrbox = stats.variables['Subset_Reliability_Rbox'][:]
            bins, fcstfreq, obhitrate = subrelrbox[:,0], subrelrbox[:,1], subrelrbox[:,2]
            probbins.append(bins)
            rel.append(obhitrate)
            binfreq.append(fcstfreq)
            colors.append('firebrick')
            fensrelrbox = stats.variables['Full_Ens_Reliability_Rbox'][:]
            bins, fcstfreq, obhitrate = fensrelrbox[:,0], fensrelrbox[:,1], fensrelrbox[:,2]
            probbins.append(bins)
            rel.append(obhitrate)
            binfreq.append(fcstfreq)
            colors.append('navy')
        probbins, rel, binfreq = np.array(probbins[:]), np.array(rel[:]), np.array(binfreq)
        mask = (rel == 9e9)
        rel = np.ma.masked_array(rel, mask=mask)
        for i in range(len(rel)):
            for j in range(len(rel[0,:,0])):
                plt.plot(probbins[0,0,:], rel[i,j,:], color=colors[i], alpha=0.5)
    plt.savefig(outpath)
    plt.close()

