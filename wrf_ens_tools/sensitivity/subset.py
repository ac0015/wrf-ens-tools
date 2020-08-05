#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 10:11:46 2018

@author: aucolema
"""
import numpy as np
import os, csv, subprocess
from datetime import timedelta
from wrf_ens_tools.sensitivity import Sens
from wrf_ens_tools.sensitivity import ensSubset
from wrf_ens_tools.post import fromDatetime, interpRAPtoWRF, subprocess_cmd
from wrf_ens_tools.post import storeReliabilityRboxFortran
from wrf_ens_tools.post import storeNearestNeighborFortran
from wrf_ens_tools.plots import plotProbs, plotDiff, plotSixPanels
from wrf_ens_tools.calc import *
from wrf_ens_tools.post import process_wrf, postTTUWRFanalysis, postIdealizedAnalysis
from netCDF4 import Dataset
import pandas as pd
# from profilehooks import profile
import xarray as xr
import time

package_dir = os.path.dirname(os.path.abspath(__file__))

#################
# Subset class
#################
class Subset:
    """
    The Subset class encapsulates all the information pertaining
    to a forecast ensemble and its subset which is gathered from
    relevant sensitivity variables and response functions and their
    respective sensitivity fields.
    """
    def __init__(self, sens=Sens(), subset_size=21, subset_method='percent',
                 percent=70., sensvalfile="SENSvals.nc", nbrhd=32.1869, thresh=25.,
                 sensvars=['500_hPa_GPH', '700_hPa_T', '850_hPa_T', 'SLP'],
                 analysis_type="RAP", wrfanalysis_to_post_path=None,
                 idealized=False, truth_member=1, semi_idealized=False):
        """
        Constructor for an instance of the Subset class.

        Inputs
        ------
        sens ---------- sensitivity obj on which the Subset object is based.
                        Uses default Sens object if none provided.
        subset_size --- integer specifying the number of members to include
                        in the subset.
        subset_method - string specifying the subsetting method. Valid options
                        are the 'point' (subset w/ max sensitivity point) method,
                        the 'weight' (aka projection) method, and the 'percent'
                        method. Defaults to percent. See Ancell 2016 for further detail.
        percent ------- sensitivity threshold as a float, specifies what percentage
                        of the sensitivity field will be used for the subset. Defaults
                        to 70 (i.e. the top 30% of sensitivity magnitudes will be used).
        sensvalfile --- string for the relative path to the netcdf file in which
                        all the values of each sensitivity variable for each member
                        are stored for its coreesponding sensitivity time. Uses the
                        path specified by default from the Sens class.
        nbrhd --------- float controlling the neighborhood distance (in km)
                        used in probability calculations. Defaults to 32.1869 km
                        or 20 miles.
        thresh -------- float controlling threshold of response function to use
                        for probabilities in verification statistics. If using
                        coverage response function, make sure this is equal to
                        the response threshold of the Sens object for an apples-
                        to-apples verification.
        sensvars ------ list of strings containing the desired sensitivity variables
                        on which to base the subset. Must comply with naming conventions
                        of keys in the RAP interpolation from interp_analysis.py.
        analysis_type - string specifying type of analysis to calculate sensitivity var
                        errors. Currently accepts "RAP" or "WRF".
        wrfanalysis_to_post_path -- optional path to half post-processed WRF analysis
                                    produced by reduced file process specific to TTU.
        """
        self._sens = sens
        self._fullens = np.arange(1,sens.getEnsnum()+1,1)
        self._subsize = subset_size
        self._methodchoices = {'point' : 1, 'weight' : 2, 'percent' : 3}
        self._percent = percent
        self._sensvars = sensvars
        self._thresh = thresh
        self._nbr = nbrhd
        self._idealized = idealized
        self._semi_idealized = semi_idealized
        # Only need this if you need to post-process an analysis file
        #  and not using RAP
        self._wrfpath_to_post_proc = wrfanalysis_to_post_path

        # If thresh and sens object response threshold don't match,
        #  issue warning
        if self._thresh != self._sens.getResponseThreshold():
            print("WARNING - Subset obj response threshold and Sensitivity obj\n",
                  "don't match. Only makes a difference if working with coverage\n",
                  "response functions.")
            print("Sens obj response threshold: ", self._sens.getResponseThreshold())
            print("Sub obj response threshold: ", self._thresh)

        # If method is not in methodchoices, set to percent (RMS) method
        try:
            self._method = self._methodchoices[subset_method]
        except:
            dflt_opt = list(self._methodchoices.keys())[-1]
            print("{} not an acceptable choice from {}. Using {} method.".format(str(subset_method),
                  self._methodchoices.keys(), dflt_opt))
            self._method = self._methodchoices[dflt_opt]
        self._subset = []

        # Define analysis file path
        self._sensdate = sens.getRunInit() + timedelta(hours=sens.getSensTime())
        yr, mo, day, hr = fromDatetime(self._sensdate, interp=True)
        if self._idealized:
            self._truth_member = truth_member
        elif self._semi_idealized:
            self._truth_member = truth_member
        self._analysis_type = analysis_type
        if self._analysis_type.upper() == "RAP":
            self._analysis = "{}RAP_interp_to_WRF_{}{}{}{}.nc".format(sens.getDir(),
                          yr, mo, day, hr)
        elif self._analysis_type.upper() == "WRF":
            self._analysis = "{}WRF_analysis_{}{}{}{}.nc".format(sens.getDir(),
                          yr, mo, day, hr)
        else:
            raise ValueError('{} not a vaid option. RAP and WRF are supported analyisis types.'.format(self._analysis_type))

        # Set sensitivity time sens variable values filepath and prob filepaths
        self._sensvalfile = sens.getDir() + sensvalfile
        self._fensprob = "FULLENSwrfout_nbr{}_f{}.prob".format(int(self._nbr), self._sens.getRTime())
        self._subprob = "SUBSETwrfout_nbr{}_f{}.prob".format(int(self._nbr), self._sens.getRTime())
        self._dx = Dataset(sens.getRefFileD2()).DX / 1000. # convert to km

    def __str__(self):
        return "Subset object with full ensemble of {} members, subset size of {}, \nand " \
               "using the {} subsetting method with a threshold of {}, \nneighborhood of {}, " \
               "response threshold of {} and these sensitivity variables: {}. \nUsing {} " \
               "analysis for subsetting. Idealized experiment = {} " \
               "\n\nBased on Sens object: \n {}".format(self._sens.getEnsnum(),
                                                    self._subsize, self._method, str(self._percent),
                                                    str(self._nbr), str(self._thresh),
                                                    ','.join(self.getSensVars()),
                                                    str(self._analysis),
                        str(self._idealized) if self._idealized==False else str(self._idealized)+"; Truth Member: {}".format(self._truth_member),
                                                    str(self.getSens()))

    def setSubsetMethod(self, subset_method):
        """
        Set the ensemble subsetting method given a
        string. Current choices are: 'point', 'weight',
        or 'percent'.
        """
        self._method = self._methodchoices[subset_method]
        return

    def setSubsetSize(self, subset_size):
        """
        Set the ensemble subset size with
        an integer.
        """
        self._subsize = subset_size
        return

    def setAnalysis(self, analysispath, analysistype):
        """
        Set the absolute path of the interpolated analysis file
        to be used for verification and set the analysis type.
        Options for analysis type are 'RAP' or 'WRF'. Make sure the
        analysispath is of the same type described in analysistype.
        """
        self._analysis = analysispath
        self._analysis_type = analysistype

    def getAnalysis(self):
        """
        Returns path to analysis file at senstime.
        """
        return self._analysis

    def getSubMembers(self):
        """
        Returns a list of the members in the subset.
        """
        return self._subset

    def getFullEns(self):
        """
        Returns a list of members of the full ensemble.
        """
        return self._fullens

    def getSens(self):
        """
        Returns Sens object valid for subset.
        """
        return self._sens

    def getSensVars(self):
        """
        Returns sensitivity variables being
        used for subsetting.
        """
        return self._sensvars

    def getHorizGridSpacingD2(self):
        """
        Returns horizontal grid spacing
        of inner domain in kilometers
        """
        return self._dx

    def interpRAP(self):
        """
        If using RAP analysis, must call this function to interpolate
        the analysis to our WRF grid before doing any subsetting. Returns
        NULL but will produce outfile with filepath as described by
        the analysis attribute of the Subset instance.
        """
        os.chdir(self.getSens().getDir())
        yr, mo, day, hr = fromDatetime(self._sensdate, interp=True)
        interpRAPtoWRF(yr, mo, day, hr, self.getSens().getRefFileD1())
        return

    def processWRFAnalysis(self, half_post_processed=False, rand_err=False,
                            rand_err_std_dev=0.):
        """
        If using WRF analysis, must call this function to post-process
        all necessary sensitivity variables for subsetting. Returns NULL
        but will produce netCDF outfile with filepath described by the
        analysis attribute of the Subset instance. The half_post_processed
        argument is associated with the post-processing technique implemented
        by Team Bancell and has a specific function for post-processing (if
        you don't know what this last part means then you won't need to use it).

        Inputs
        ------
        self ---------------- instance of Subset object
        half_post_processed - boolean specific to Brian's TTU WRF analysis
                              post-processing, if not working with that dataset
                              always leave False
        rand_err ------------ optional boolean specifying whether to add random
                              error to analysis. Only works with idealized
                              settings
        rand_err_std_dev ---- if rand_err is set to True, may provide the
                              standard deviation of observation error to add to
                              the analysis. Adds this to ALL variables in the
                              analysis, so only recommended for specific use
                              cases

        Outputs
        -------
        Returns NULL but produces netCDF with naming conventions described by
        the analysis attribute of the Subset object
        """
        if half_post_processed:
            postTTUWRFanalysis(self._wrfpath_to_post_proc, self._analysis)
        elif self._idealized or self._semi_idealized:
            postIdealizedAnalysis(self._sensvalfile, self._analysis,
                                self._truth_member, rand_err=rand_err,
                                rand_err_std_dev=rand_err_std_dev)
        else:
            # Use default vars for process_wrf and default naming conventions.
            process_wrf(self._wrfpath_to_post_proc, outpath=self._analysis, reduced=True)

        return

    def calcSubset(self):
        """
        Calls the ensSubset() function from esens_subsetting.py.
        Returns the subset members based on the subset technique
        parameters of the Subset obj.
        """
        S = self.getSens()
        print(os.path.isfile(self._analysis))
        if (os.path.isfile(self._analysis) == False):
            if self._analysis_type == "RAP":
                print("{} does not exist. Running RAP interpolation.".format(self._analysis))
                self.interpRAP()
            elif self._analysis_type == "WRF":
                print("{} does not exist. Running WRF post-process.".format(self._analysis))
                if self._wrfpath_to_post_proc is not None:
                    half_post_processed = True
                else:
                    half_post_processed = False
                self.processWRFAnalysis(half_post_processed)
        if (os.path.isfile(S.getWRFSensFile()) == False):
            print("{} does not exist. Running sensitivity module with Sens obj.".format(S.getWRFSensFile()))
            S.runSENS()
        if (os.path.isfile(self._sensvalfile) == False):
            print("{} does not exist. Running store sens vals with Sens obj.".format(self._sensvalfile))
            S.storeSENSvals()
        self._subset, sensstrings = ensSubset(S.getWRFSensFile(), self._analysis,
                                              self._sensvalfile, S.getEnsnum(),
                                              self._subsize, method=self._method,
                                              sens_thresh=self._percent,
                                              sensvars=self.getSensVars())
        # Save subset members in textfile
        np.savetxt(S.getDir() + 'subset_members.txt', np.array(self._subset, dtype=int))
        return

    def calcProbs(self, members):
        """
        Runs the fortran executable calcprobSUBSET or sixhrprobcalc to calculate
        probabilities for any number of ensemble members. Takes
        an input file with ensemble number, ensemble members,
        response time, neighborhood, and prob output path.
        """
        S = self.getSens()

        ############### calcProbs-specific error-handling methods. ############
        def checkSuccess(sens_obj, args):
            """
            Checks to see if the word "SUCCESSFUL"
            appears in the probcalc log file, indicating
            the probability calculation was completed
            correctly.
            """
            with open('probs.out') as out:
                lowcase_out = [o.lower() for o in out]
                success = (np.array(['successful' in s for s in lowcase_out]).any())
            return success

        def reRun(probpath, sens_obj, args):
            """
            Re-runs calcProbs.
            """
            if os.path.exists(probpath):
                os.popen('rm {}'.format(probpath))
            os.popen("cp {} {}".format(sens_obj.getRefFileD2(), probpath))
            subprocess_cmd(args)
            return
        ############## End Error-Handling Methods ############################

        # Create or navigate into probs directory
        probdir = S.getDir() + "probs/"
        direxists = os.path.exists(probdir)
        if (direxists == False):
            os.mkdir(probdir)
        os.chdir(probdir)
        print("Calculating probs for: ", members)
        # If semi-idealized, remove truth member from membership
        if self._semi_idealized or self._idealized:
            self._fullens = self._fullens[self._fullens != self._truth_member]
        if len(members) == len(self._fullens):
            fname = "fullens_probs.in"
            probout = self._fensprob
        else:
            fname = "subset_probs.in"
            probout = self._subprob

        # Create list of args to create input file with
        args = ["{}/create_probcalc_inputfile.bash".format(package_dir)]
        args.append(fname)
        args.append(str(len(members)))
        args.append(' '.join([str(mem) for mem in members]))
        args.append(str(S.getRTime()))
        args.append(str(self._nbr))
        args.append(probout)
        print(args)
        subprocess.check_call(args)

        # Initialize probfile
        if os.path.exists(probout):
            os.popen('rm {}'.format(probout))
        os.popen("cp {} {}".format(S.getRefFileD2(), probout))

        # Run fortran executable
        if S.getSixHour():
            print("Calculating six hour probabilities...")
            probcalc_exec = "sixhrprobcalc"
        else:
            print("Calculating one hour probabilities...")
            probcalc_exec = "probcalcSUBSET"
        # Build command
        probcalcpath = os.path.join(package_dir, probcalc_exec)
        args = "{} <{} >probs.out".format(probcalcpath, fname)
        # Execute
        try:
            subprocess_cmd(args)
        except OSError:
            # Try again - sometimes this executable is finicky
            try:
                reRun(probout, S, args)
            except:
                print('Probcalc failed twice. Continue to error checks.')

        # Check to make sure it ran correctly
        if os.path.exists(probout) == False:
            print('Probcalc failed. Re-running.')
            reRun(probout, S, args)

        # This can get stuck in an infinite loop if something is
        # fundamentally wrong with the configuration... But it's
        # a risk we're gonna have to take
        success = checkSuccess(S, args)
        while success == False:
            print(args)
            reRun(probout, S, args)
            success = checkSuccess(S, args)

        return

    def researchPlotProbs(self, use_subset):
        """
        Calls plotProbs from research_plotting library. Passes
        mainly file paths and some metadata and function will
        handle the rest.
        """
        S = self.getSens()
        direc = S.getDir() + "probs/"
        fullenspath = direc + self._fensprob
        subsetpath = direc + self._subprob
        if use_subset:
            path = subsetpath
        else:
            path = fullenspath
        wrfrefpath = S.getRefFileD2()

        plotProbs(path, wrfrefpath,
                   S.getRbox(), S.getRTime(), self._nbr, outpath=direc,
                   subset=use_subset)
        return

    def researchPlotDiffs(self, verif_day=False):
        """
        Calculates and plots delta probabilities between the 1-hr full ensemble
        probs and its corresponding 1-hr subset probs. Calls plotDiff from
        reserch_plotting library. Passes only file paths and outside function
        will handle the rest. If verif_day is set to True, will overlay
        full day's storm reports onto difference plot.
        """
        S = self.getSens()
        direc = S.getDir() + 'probs/'
        fullenspath = direc + self._fensprob
        subsetpath = direc + self._subprob
        wrfrefpath = S.getRefFileD2()

        # Format date for storm reports
        date = S.getRunInit()
        SPCdate = str(date)[2:10].replace('-','')

        plotDiff(fullenspath, subsetpath, wrfrefpath,
                 S.getRbox(), S.getRTime(), SPCdate,
                 stormreports=verif_day)
        return

    def plotSixPanels(self, storm_reports=True):
        S = self.getSens()
        # Make sure subset has already been calculated
        yr, mo, day, hr = fromDatetime(S.getRunInit(), interp=False)
        rfunclabel = S.getRString().replace(' ','').lower()
        dirdate = str(yr) + str(mo) + str(day) + str(hr)
        plotSixPanels(dirdate, storm_reports,
                      self.getSubMembers(), sixhour=False,
                      time=S.getRTime(), subsettype=rfunclabel,
                      nbrhd=self._nbr)
        return

    # def storeUHStatsCSV(self, outpath, pperfpath, reliabilityobpath):
    #     """
    #     Stores UH verification stats to csv file.
    #
    #     Inputs
    #     ------
    #     outpath ----------- absolute output filepath as a string
    #     pperfpath --------- absolute filepath for practically perfect
    #                         gridded data
    #     reliabilityobpath - absolute filepath for reliability observational
    #                         point data on ensemble grid
    #
    #     Outputs
    #     -------
    #     returns NULL and stores UH verification stats to outpath as csv
    #     """
    #     S = self.getSens()
    #     fensprobpath = S.getDir() + "probs/" + self._fensprob
    #     subprobpath = S.getDir() + "probs/" + self._subprob
    #     rtimedate = S.getRunInit() + timedelta(hours=S.getRTime())
    #     if os.path.exists(pperfpath):
    #         fens_fss = FSSnetcdf(fensprobpath, pperfpath, rtimedate, var='updraft_helicity',
    #                   thresh=self._thresh, rboxpath=S.getDir()+'esens.in')
    #         sub_fss = FSSnetcdf(subprobpath, pperfpath, rtimedate, var='updraft_helicity',
    #                   thresh=self._thresh, rboxpath=S.getDir()+'esens.in')
    #         fens_reliability = ReliabilityRbox(fensprobpath, S.getRunInit(),
    #                                       S.getRTime(),
    #                                       obpath=reliabilityobpath,
    #                                       var='updraft_helicity',
    #                                       thresh=self._thresh,
    #                                       rboxpath=S.getDir()+'esens.in',
    #                                       sixhr=False, nbrhd=self._nbr)
    #         sub_reliability = ReliabilityRbox(subprobpath,
    #                                       S.getRunInit(), S.getRTime(),
    #                                       obpath=reliabilityobpath,
    #                                       var='updraft_helicity',
    #                                       thresh=self._thresh,
    #                                       rboxpath=S.getDir()+'esens.in',
    #                                       sixhr=False, nbrhd=self._nbr)
    #         # prob_bins will stay the same, so OK to clobber
    #         prob_bins, f_fcstfreq_rbox, f_ob_hr_rbox = fens_reliability
    #         prob_bins, s_fcstfreq_rbox, s_ob_hr_rbox = sub_reliability
    #         fens_rel_dict = {'rel_bins' : prob_bins,
    #                         'fcst_freq_rbox' : f_fcstfreq_rbox,
    #                         'ob_hr_rbox' : f_ob_hr_rbox}
    #         sub_rel_dict = {'rel_bins' : prob_bins,
    #                         'fcst_freq_rbox' : s_fcstfreq_rbox,
    #                         'ob_hr_rbox' : s_ob_hr_rbox}
    #         # fens_rel_dict = np.asarray([prob_bins, f_fcstfreq_rbox,
    #         #                             f_ob_hr_rbox], dtype=float)
    #         # sub_rel_dict = np.asarray([prob_bins, s_fcstfreq_rbox,
    #         #                             s_ob_hr_rbox], dtype=float)
    #         f_fss_tot, f_fss_rbox, sig = fens_fss
    #         s_fss_tot, s_fss_rbox, sig = sub_fss
    #     else:
    #         raise FileNotFoundError('Please run storePracPerf() or set correct obpath.')
    #     if os.path.exists(outpath):
    #         mode = 'a'
    #     else:
    #         mode = 'w'
    #     # Open file with specified mode
    #     with open(outpath, mode=mode) as outfile:
    #         if mode == 'w':
    #             # If new file, add header row
    #             cols = ['Run_Init', 'Sens_Time', 'Subset_Size',
    #                     'Analysis', 'Sens_Vars', 'Subset_Method',
    #                     'Sens_Threshold', 'Response_Func',
    #                     'Response_Time', 'Response_Box',
    #                     'Full_Ens_FSS_Total', 'Full_Ens_FSS_Rbox',
    #                     'Subset_FSS_Total', 'Subset_FSS_Rbox',
    #                     'Prac_Perf_Sigma', 'Neighborhood',
    #                     'Response_Thresh', 'Sub_Members',
    #                     'Full_Ens_Reliability_Rbox',
    #                     'Subset_Reliability_Rbox']
    #             outfile_writer = csv.writer(outfile, delimiter=',')
    #             outfile_writer.writerow(cols)
    #         # Add row valid for current subset
    #         entry = [str(S.getRunInit()), S.getSensTime(), self._subsize,
    #                  self._analysis_type, self._sensvars,
    #                  list(self._methodchoices.keys())[self._method-1],
    #                  self._percent, S.getRString(), S.getRTime(), S.getRbox(),
    #                  f_fss_tot, f_fss_rbox, s_fss_tot, s_fss_rbox,
    #                  sig[0], self._nbr, self._thresh, self.getSubMembers(),
    #                  fens_rel_dict, sub_rel_dict]
    #         outfile_writer = csv.writer(outfile, delimiter=',')
    #         outfile_writer.writerow(entry)
    #     return

    def storeUHStatsNetCDF(self, outpath, pperfpath, reliabilityobpath,
                            bssobpath):
        """
        Stores UH verification stats to netCDF4 file using xarray.

        Inputs
        ------
        outpath ----------- absolute output filepath as a string
        pperfpath --------- absolute filepath for practically perfect
                            gridded data
        reliabilityobpath - absolute filepath for reliability observational
                            point data on ensemble grid
        bssobpath --------- absolute filepath for Brier skill score
                            observation grid

        Outputs
        -------
        returns NULL and stores UH verification stats to outpath as netCDF4
        """
        S = self.getSens()
        fensprobpath = S.getDir() + "probs/" + self._fensprob
        subprobpath = S.getDir() + "probs/" + self._subprob
        outrel = S.getDir() + "reliability_out.nc"
        rtimedate = S.getRunInit() + timedelta(hours=S.getRTime())

        def checkSuccess():
            """
            Checks to see if the word "SUCCESSFUL"
            appears in the reliabilitycalc log file, indicating
            the reliability calculation was completed
            correctly.
            """
            with open(S.getDir()+'reliability.out') as out:
                lowcase_out = [o.lower() for o in out]
                success = (np.array(['successful' in s for s in lowcase_out]).any())
            return success

        if os.path.exists(pperfpath):
            fens_fss = FSSnetcdf(fensprobpath, pperfpath, rtimedate, var='updraft_helicity',
                      thresh=self._thresh, rboxpath=S.getDir()+'esens.in')
            sub_fss = FSSnetcdf(subprobpath, pperfpath, rtimedate, var='updraft_helicity',
                      thresh=self._thresh, rboxpath=S.getDir()+'esens.in')
            f_fss_tot, f_fss_rbox, sig = fens_fss
            s_fss_tot, s_fss_rbox, sig = sub_fss
            # If fortran-produced reliability file doesn't exist
            #  then create it.
            if os.path.exists(reliabilityobpath) == False:
                print(reliabilityobpath + " does not exist...")
                print("Running storeNearestNeighborFortran()...")
                storeNearestNeighborFortran(S.getRunInit(), S.getRTime(),
                                                reliabilityobpath,
                                                sixhour=S.getSixHour(),
                                                wrfrefpath=S.getRefFileD2())
            # Process the reliability observation file by calculating
            #  reliabilty in fortran - those results are stored to
            #  a netCDF file called 'reliability_out.nc'
            print("Processing reliability of {} with Fortran...".format(subprobpath))
            success = False
            i = 0
            while (success == False):
                storeReliabilityRboxFortran(S.getDir(), S.getRTime(),
                            subprobpath, reliabilityobpath, outrel,
                            rboxpath=S.getDir()+'esens.in',
                            sixhour=S.getSixHour(),
                            variable="updraft_helicity",
                            rthresh=self._thresh, nbrhd=self._nbr,
                            wrfrefpath=S.getRefFileD2())
                time.sleep(.300)
                print(i)
                success = checkSuccess()
                i += 1
            # Process result with xarray
            sub_reliability = xr.open_dataset(outrel)
            prob_bins = sub_reliability["prob_bins"]
            s_fcstfreq_rbox = sub_reliability["fcst_frequency"]
            s_ob_hr_rbox = sub_reliability["ob_hit_rate"]
            sub_reliability.close()
            ####################################################################
            # Calculate BSS using probabilistic forecast and binary observations
            # Same as FSS except binary obs
            ####################################################################
            fens_prob_ds = xr.open_dataset(fensprobpath)
            sub_prob_ds = xr.open_dataset(subprobpath)
            # Choose correct indices based on variable and threshold
            probinds = {'reflectivity': {40: 0, 50: 5},
                        'updraft_helicity': {25: 1, 40: 2, 100: 3},
                        'wind_speed': {40: 4}}
            index = probinds['updraft_helicity'][self._thresh]
            # Pull probabilities
            fens_probs = np.asarray(fens_prob_ds["P_HYD"])[0, index]
            sub_probs = np.asarray(sub_prob_ds["P_HYD"])[0, index]
            fens_prob_ds.close(); sub_prob_ds.close()
            # Open binary observation dataset
            bss_ob_ds = xr.open_dataset(bssobpath)
            bss_obs = np.asarray(bss_ob_ds["nearest_neighbor"][0])
            bss_ob_ds.close()
            # Mask obs/fcst
            ref_ds = xr.open_dataset(S.getRefFileD2())
            lons = np.asarray(ref_ds["XLONG"].values[0], dtype=float)
            lats = np.asarray(ref_ds["XLAT"].values[0], dtype=float)
            ref_ds.close()
            llon, ulon, llat, ulat = np.array(S.getRbox(), dtype=float)
            lonmask = (lons > llon) & (lons < ulon)
            latmask = (lats > llat) & (lats < ulat)
            mask = lonmask & latmask
            del(lons); del(lats)
            masked_fens_probs = np.ma.masked_array(fens_probs/100., mask=~mask)
            masked_sub_probs = np.ma.masked_array(sub_probs/100., mask=~mask)
            masked_obs = np.ma.masked_array(bss_obs, mask=~(mask & (bss_obs <= 9e36)))
            # print("NANMAX BSS OBS", np.nanmax(masked_obs))
            # print("Rbox probs", masked_fens_probs)
            # Calculate BSS!
            f_bss_rbox = FSS(masked_fens_probs, masked_obs)
            print(f"Full Ens BSS: {f_bss_rbox}")
            s_bss_rbox = FSS(masked_sub_probs, masked_obs)
            print(f"Subset BSS: {s_bss_rbox}")
            ####################################################################
        else:
            raise FileNotFoundError('Please run storePracPerf() or \
            set correct obpath.')
        # Create entry as xarray dataset
        sens_vars = np.empty((1, 30), dtype="U30")
        sens_vars[0, :len(self._sensvars)] = np.asarray(self._sensvars[:])
        sub_mems = np.zeros((1,len(self._fullens)), dtype=int) * np.NaN
        sub_mems[0, :len(self.getSubMembers())] = np.asarray(self.getSubMembers())
        ds = xr.Dataset({'Sens_Time': (['subset'], np.atleast_1d(S.getSensTime())),
                 'Subset_Size': (['subset'], np.atleast_1d(self._subsize)),
                 'Analysis': (['subset'], np.atleast_1d(self._analysis_type)),
                 'Sens_Vars':  (['subset', 'sensvar'],
                                np.atleast_2d(sens_vars)),
                 'Subset_Method': (['subset'],
                    np.atleast_1d(list(self._methodchoices.keys())[self._method-1])),
                 'Sens_Threshold': (['subset'], np.atleast_1d(self._percent)),
                 'Response_Func': (['subset'], np.atleast_1d(S.getRString())),
                 'Response_Time': (['subset'], np.atleast_1d(S.getRTime())),
                 'Response_Box': (['subset', 'rbox'],
                                    np.atleast_2d(S.getRbox())),
                 'Full_Ens_FSS_Total': (['subset'], np.atleast_1d(f_fss_tot)),
                 'Full_Ens_FSS_Rbox': (['subset'], np.atleast_1d(f_fss_rbox)),
                 'Subset_FSS_Total': (['subset'], np.atleast_1d(s_fss_tot)),
                 'Subset_FSS_Rbox': (['subset'], np.atleast_1d(s_fss_rbox)),
                 'Prac_Perf_Sigma': (['subset'], np.atleast_1d(sig[0])),
                 'Full_Ens_BSS_Rbox': (['subset'], np.atleast_1d(f_bss_rbox)),
                 'Subset_BSS_Rbox': (['subset'], np.atleast_1d(s_bss_rbox)),
                 'Neighborhood': (['subset'], np.atleast_1d(self._nbr)),
                 'Response_Thresh': (['subset'], np.atleast_1d(self._thresh)),
                 'Subset_Fcst_Freq_Rbox': (['subset', 'prob_bins'],
                                            np.atleast_2d(s_fcstfreq_rbox)),
                 'Subset_Ob_Hit_Rate_Rbox': (['subset', 'prob_bins'],
                                             np.atleast_2d(s_ob_hr_rbox)),
                 'Subset_Members': (['subset', 'submems'],
                                    sub_mems)},
                 coords={'run': S.getRunInit(),
                         # 'subset': np.arange(1),
                         'bins': prob_bins,
                         'rbox': ['llon', 'ulon', 'llat', 'ulat']})
        if os.path.exists(outpath):
            og_ds = xr.open_dataset(outpath)
            print("Opened original dataset...")
            # Check first to see if we need to add full ens reliability
            print("Current Rthresh",
                self._thresh, "Current Full Ens Thresh Vals",
                og_ds.Full_Ens_Response_Thresh.values)
            if (self._thresh not in og_ds.Full_Ens_Response_Thresh.values):
                print("Processing full ens reliability for new threshold...")
                success = False
                while (success == False):
                    storeReliabilityRboxFortran(S.getDir(), S.getRTime(),
                                fensprobpath, reliabilityobpath, outrel,
                                rboxpath=S.getDir()+'esens.in',
                                sixhour=S.getSixHour(),
                                variable="updraft_helicity",
                                rthresh=self._thresh, nbrhd=self._nbr,
                                wrfrefpath=S.getRefFileD2())
                    time.sleep(.300)
                    success = checkSuccess()
                fens_reliability = xr.open_dataset(outrel)
                prob_bins = fens_reliability["prob_bins"]
                f_fcstfreq_rbox = fens_reliability["fcst_frequency"]
                f_ob_hr_rbox = fens_reliability["ob_hit_rate"]
                print("New threshold full ob hit rate", f_ob_hr_rbox)
                new_ds = xr.Dataset({'Full_Ens_Response_Thresh': (['rthresh'],
                                            np.atleast_1d(self._thresh)),
                                    'Full_Ens_Fcst_Freq_Rbox': (['rthresh', 'prob_bins'],
                                            np.atleast_2d(f_fcstfreq_rbox)),
                                    'Full_Ens_Ob_Hit_Rate_Rbox': (['rthresh', 'prob_bins'],
                                            np.atleast_2d(f_ob_hr_rbox))},
                                    coords={'run': S.getRunInit(),
                                            'bins': prob_bins,
                                            'rbox': ['llon', 'ulon', 'llat', 'ulat']})
                print(new_ds)
                fens_reliability.close()
                og_fens_ds = xr.concat([og_ds, new_ds], dim='rthresh',
                                data_vars=['Full_Ens_Response_Thresh',
                                'Full_Ens_Fcst_Freq_Rbox',
                                'Full_Ens_Ob_Hit_Rate_Rbox'])
                new_ds.close()
                print("Concatenated rel to appending dataset...\n", og_fens_ds)
                # Finally concatenate with existing netCDF dataset
                ds = xr.concat([og_fens_ds, ds], dim='subset', data_vars='minimal')
                og_fens_ds.close()
            else:
                ds = xr.concat([og_ds, ds], dim='subset', data_vars='minimal')
            og_ds.close()
            print("Appended new dataset to original dataset...\n\n", ds)
            # os.popen("rm {}".format(outpath))
        else:
            print("Adding full ens rel to ", ds)
            success = False
            while (success == False):
                storeReliabilityRboxFortran(S.getDir(), S.getRTime(),
                            fensprobpath, reliabilityobpath, outrel,
                            rboxpath=S.getDir()+'esens.in',
                            sixhour=S.getSixHour(),
                            variable="updraft_helicity",
                            rthresh=self._thresh, nbrhd=self._nbr,
                            wrfrefpath=S.getRefFileD2())
                time.sleep(.300)
                success = checkSuccess()
            fens_reliability = xr.open_dataset(outrel)
            prob_bins = fens_reliability["prob_bins"]
            f_fcstfreq_rbox = fens_reliability["fcst_frequency"]
            f_ob_hr_rbox = fens_reliability["ob_hit_rate"]
            fens_reliability.close()
            ds["Full_Ens_Response_Thresh"] = (('rthresh'),
                                        np.atleast_1d(self._thresh))
            ds["Full_Ens_Fcst_Freq_Rbox"] = (('rthresh', 'prob_bins'),
                                        np.atleast_2d(f_fcstfreq_rbox))
            ds["Full_Ens_Ob_Hit_Rate_Rbox"] = (('rthresh', 'prob_bins'),
                                        np.atleast_2d(f_ob_hr_rbox))
        print("FINAL DATASET TO WRITE:", ds)
        ds.to_netcdf("test.nc", unlimited_dims=['subset', 'rthresh'],
                    format="NETCDF4", mode='w')
        os.rename("test.nc", outpath)
        ds.close()

        return

    def storeReflStatsNetCDF(self, outpath, gridradfiles,
                            reliabilityobpath, og_gridradfiles=None):
        """
        Function for storing subset simulated reflectivity verification
        statistics using GridRad as the observation dataset. Currently
        stores stats to NetCDF using xarray datasets. Verifies response
        box reflectivity maxima and coverage.

        Inputs
        ------
        outpath ----------- absolute filename for statistics output
        gridradfiles ------ list of interpolated GridRad filepaths
        reliabilityobpath - file to store reliability observations
                            (binary hits/misses given a reflectivity
                            threshold)
        og_gridradfiles --- (optional) only need if evaluating refl
                            maxima - list of absolute paths to
                            original GridRad files

        Outputs
        -------
        returns NULL but will store an xarray dataset with verification stats
        to netCDF using the specified outpath parameter.
        """
        S = self.getSens()
        fensprobpath = S.getDir() + "probs/" + self._fensprob
        subprobpath = S.getDir() + "probs/" + self._subprob
        rvalspath = S.getDir() + "Rvals.nc"
        outrel = S.getDir() + "reliability_out_refl.nc"
        rtimedate = S.getRunInit() + timedelta(hours=S.getRTime())

        def checkSuccess():
            """
            Checks to see if the word "SUCCESSFUL"
            appears in the reliabilitycalc log file, indicating
            the reliability calculation was completed
            correctly.
            """
            with open(S.getDir()+'reliability.out') as out:
                lowcase_out = [o.lower() for o in out]
                success = (np.array(['successful' in s for s in lowcase_out]).any())
            return success

        if os.path.exists(fensprobpath):
            #fens_fss = FSSnetcdf(fensprobpath, pperfpath, rtimedate, var='reflectivity',
            #          thresh=self._thresh, rboxpath=S.getDir()+'esens.in')
            #sub_fss = FSSnetcdf(subprobpath, pperfpath, rtimedate, var='reflectivity',
            #          thresh=self._thresh, rboxpath=S.getDir()+'esens.in')
            #f_fss_tot, f_fss_rbox, sig = fens_fss
            #s_fss_tot, s_fss_rbox, sig = sub_fss
            if "Coverage" in S.getRString():
                resp = "DBZ_COV"
                dbz_ob, npts_rbox = calc_refl_cov_rbox(interpgridradfiles=gridradfiles,
                                            rboxpath=S.getDir()+'esens.in',
                                            zlev=0, refl_thresh=self._thresh)
            elif "Max" in S.getRString():
                resp = "DBZ_MAX"
                dbz_ob = calc_refl_max_rbox(gridradfiles=og_gridradfiles,
                                            rboxpath=S.getDir()+'esens.in',
                                            zlev=0)
            # fens_avg_rval_rbox = calc_subset_avg_response_rbox(member_list=self._fullens,
            #                             rvalues_ncfile=rvalspath,
            #                             rfuncstr=resp)
            # sub_avg_rval_rbox = calc_subset_avg_response_rbox(member_list=self._subset,
            #                             rvalues_ncfile=rvalspath,
            #                             rfuncstr=resp)
            # print("Full Ensemble Average {} Rbox: {}".format(resp, fens_avg_rval_rbox))
            # print("Subset Average {} Rbox: {}".format(resp, sub_avg_rval_rbox))
            rstring_to_rindex = {"6-hr Max Refl": "DBZ_MAX",
                                 "6-hr Refl Coverage": "DBZ_COV"}
            rvals_dat = xr.open_dataset(rvalspath)
            sub_mae = np.zeros_like(self.getSubMembers())
            fens_mae = np.zeros_like(self._fullens)
            # Calculate MAE for each subset member
            for ind, mem in enumerate(self.getSubMembers()):
                print("Subset Member:", mem)
                mem_rval = rvals_dat[rstring_to_rindex[S.getRString()]][mem-1]
                sub_mae[ind] = abs(mem_rval - dbz_ob)
            # Calculate MAE for each full ensemble member
            for ind, fens_mem in enumerate(self._fullens):
                print("Full Ens Mem:", fens_mem)
                mem_rval = rvals_dat[rstring_to_rindex[S.getRString()]][fens_mem-1]
                fens_mae[ind] = abs(mem_rval - dbz_ob)
            print("{} Observation Rbox: {}".format(resp, dbz_ob))
            fens_abs_err = np.mean(fens_mae)
            sub_abs_err = np.mean(sub_mae)
            print("Subset {} error: {}".format(S.getRString(), sub_abs_err))
            print("Full Ensemble {} error: {}".format(S.getRString(), fens_abs_err))

            # If fortran-produced reliability file doesn't exist
            #  then create it.
            if os.path.exists(reliabilityobpath) == False:
                print(reliabilityobpath + " does not exist...")
                print("Running storeNearestNeighborFortran()...")
                storeNearestNeighborFortran(S.getRunInit(), S.getRTime(),
                                                outpath=reliabilityobpath,
                                                variable='reflectivity',
                                                interp_gridradfiles=gridradfiles,
                                                reflthreshold=self._thresh,
                                                sixhour=S.getSixHour(),
                                                wrfrefpath=S.getRefFileD2())
            # Process the reliability observation file by calculating
            #  reliabilty in fortran - those results are stored to
            #  a netCDF file called 'reliability_out.nc'
            print("Processing reliability of {} with Fortran...".format(subprobpath))
            success = False
            i = 0
            while (success == False):
                storeReliabilityRboxFortran(S.getDir(), S.getRTime(),
                            subprobpath, reliabilityobpath, outrel,
                            rboxpath=S.getDir()+'esens.in',
                            sixhour=S.getSixHour(),
                            variable="reflectivity",
                            rthresh=self._thresh, nbrhd=self._nbr,
                            wrfrefpath=S.getRefFileD2())
                print(i)
                time.sleep(.300)
                success = checkSuccess()
                i += 1
            # Process result with xarray
            sub_reliability = xr.open_dataset(outrel)
            prob_bins = sub_reliability["prob_bins"]
            s_fcstfreq_rbox = sub_reliability["fcst_frequency"]
            s_ob_hr_rbox = sub_reliability["ob_hit_rate"]
            sub_reliability.close()
        else:
            raise FileNotFoundError('Please run storePracPerf() or \
            set correct obpath.')
        # Create entry as xarray dataset
        sens_vars = np.empty((1, 30), dtype="U30")
        sens_vars[0, :len(self._sensvars)] = np.asarray(self._sensvars[:])
        sub_mems = np.zeros((1,len(self._fullens)), dtype=int) * np.NaN
        sub_mems[0, :len(self.getSubMembers())] = np.asarray(self.getSubMembers())
        ds = xr.Dataset({'Sens_Time': (['subset'], np.atleast_1d(S.getSensTime())),
                 'Subset_Size': (['subset'], np.atleast_1d(self._subsize)),
                 'Analysis': (['subset'], np.atleast_1d(self._analysis_type)),
                 'Sens_Vars':  (['subset', 'sensvar'],
                                np.atleast_2d(sens_vars)),
                 'Subset_Method': (['subset'],
                    np.atleast_1d(list(self._methodchoices.keys())[self._method-1])),
                 'Sens_Threshold': (['subset'], np.atleast_1d(self._percent)),
                 'Response_Func': (['subset'], np.atleast_1d(S.getRString())),
                 'Response_Time': (['subset'], np.atleast_1d(S.getRTime())),
                 'Response_Box': (['subset', 'rbox'],
                                    np.atleast_2d(S.getRbox())),
                 'Full_Ens_MAE_Rbox' : (['subset'], np.atleast_1d(fens_abs_err)),
                 'Subset_MAE_Rbox': (['subset'], np.atleast_1d(sub_abs_err)),
                 'Neighborhood': (['subset'], np.atleast_1d(self._nbr)),
                 'Response_Thresh': (['subset'], np.atleast_1d(self._thresh)),
                 'Subset_Fcst_Freq_Rbox': (['subset', 'prob_bins'],
                                            np.atleast_2d(s_fcstfreq_rbox)),
                 'Subset_Ob_Hit_Rate_Rbox': (['subset', 'prob_bins'],
                                             np.atleast_2d(s_ob_hr_rbox)),
                 'Subset_Members': (['subset', 'submems'],
                                    sub_mems)},
                 coords={'run': S.getRunInit(),
                         # 'subset': np.arange(1),
                         'bins': prob_bins,
                         'rbox': ['llon', 'ulon', 'llat', 'ulat']})
        if os.path.exists(outpath):
            og_ds = xr.open_dataset(outpath)
            print("Opened original dataset...")
            # Check first to see if we need to add full ens reliability
            print("Current Rthresh",
                self._thresh, "Current Full Ens Thresh Vals",
                og_ds.Full_Ens_Response_Thresh.values)
            if (self._thresh not in og_ds.Full_Ens_Response_Thresh.values):
                print("Processing full ens reliability for new threshold...")
                success = False
                while (success == False):
                    storeReliabilityRboxFortran(S.getDir(), S.getRTime(),
                                fensprobpath, reliabilityobpath, outrel,
                                rboxpath=S.getDir()+'esens.in',
                                sixhour=S.getSixHour(),
                                variable="reflectivity",
                                rthresh=self._thresh, nbrhd=self._nbr,
                                wrfrefpath=S.getRefFileD2())
                    time.sleep(.300)
                    success = checkSuccess()
                fens_reliability = xr.open_dataset(outrel)
                prob_bins = fens_reliability["prob_bins"]
                f_fcstfreq_rbox = fens_reliability["fcst_frequency"]
                f_ob_hr_rbox = fens_reliability["ob_hit_rate"]
                print("New threshold full ob hit rate", f_ob_hr_rbox)
                new_ds = xr.Dataset({'Full_Ens_Response_Thresh': (['rthresh'],
                                            np.atleast_1d(self._thresh)),
                                    'Full_Ens_Fcst_Freq_Rbox': (['rthresh', 'prob_bins'],
                                            np.atleast_2d(f_fcstfreq_rbox)),
                                    'Full_Ens_Ob_Hit_Rate_Rbox': (['rthresh', 'prob_bins'],
                                            np.atleast_2d(f_ob_hr_rbox))},
                                    coords={'run': S.getRunInit(),
                                            'bins': prob_bins,
                                            'rbox': ['llon', 'ulon', 'llat', 'ulat']})
                print(new_ds)
                fens_reliability.close()
                og_fens_ds = xr.concat([og_ds, new_ds], dim='rthresh',
                                data_vars=['Full_Ens_Response_Thresh',
                                'Full_Ens_Fcst_Freq_Rbox',
                                'Full_Ens_Ob_Hit_Rate_Rbox'])
                new_ds.close()
                print("Concatenated rel to appending dataset...\n", og_fens_ds)
                # Finally concatenate with existing netCDF dataset
                ds = xr.concat([og_fens_ds, ds], dim='subset', data_vars='minimal')
                og_fens_ds.close()
            else:
                ds = xr.concat([og_ds, ds], dim='subset', data_vars='minimal')
            og_ds.close()
            print("Appended new dataset to original dataset...\n\n", ds)
            # os.popen("rm {}".format(outpath))
        else:
            print("Adding full ens rel to ", ds)
            success = False
            while (success == False):
                storeReliabilityRboxFortran(S.getDir(), S.getRTime(),
                            fensprobpath, reliabilityobpath, outrel,
                            rboxpath=S.getDir()+'esens.in',
                            sixhour=S.getSixHour(),
                            variable="reflectivity",
                            rthresh=self._thresh, nbrhd=self._nbr,
                            wrfrefpath=S.getRefFileD2())
                time.sleep(.300)
                success = checkSuccess()
            fens_reliability = xr.open_dataset(outrel)
            prob_bins = fens_reliability["prob_bins"]
            f_fcstfreq_rbox = fens_reliability["fcst_frequency"]
            f_ob_hr_rbox = fens_reliability["ob_hit_rate"]
            fens_reliability.close()
            ds["Full_Ens_Response_Thresh"] = (('rthresh'),
                                        np.atleast_1d(self._thresh))
            ds["Full_Ens_Fcst_Freq_Rbox"] = (('rthresh', 'prob_bins'),
                                        np.atleast_2d(f_fcstfreq_rbox))
            ds["Full_Ens_Ob_Hit_Rate_Rbox"] = (('rthresh', 'prob_bins'),
                                        np.atleast_2d(f_ob_hr_rbox))
        print("FINAL DATASET TO WRITE:", ds)
        ds.to_netcdf("test.nc", unlimited_dims=['subset', 'rthresh'],
                    format="NETCDF4", mode='w')
        os.rename("test.nc", outpath)
        ds.close()

        return

    def idealizedStoreRfuncStatsNetCDF(self, outpath, rvalspath="Rvals.nc"):
        """
        Stores response function verification stats for an idealized experiment
        in which a randomly selected ensemble member serves as truth and is
        used to compute response errors directly.

        Inputs
        ------
        outpath ----------- absolute output filepath as a string
        rvalspath --------- name of netCDF file produced from
                            sixhresens that stores member response values

        Outputs
        -------
        returns NULL and stores UH verification stats to outpath as netCDF4
        """
        if self._idealized:
            S = self.getSens()
            rstring_to_rindex = {"6-hr Max UH": "UH_MAX",
                                 "6-hr UH Coverage": "UH_COV",
                                 "6-hr Max Refl": "DBZ_MAX",
                                 "6-hr Refl Coverage": "DBZ_COV"}
            truthmem_ind = self._truth_member - 1
            # Open member response values
            rvals_dat = Dataset(S.getDir() + rvalspath)
            # Create entry as xarray dataset
            sens_vars = np.empty((1, 30), dtype="U30")
            sens_vars[0, :len(self._sensvars)] = np.asarray(self._sensvars[:])
            sub_mems = np.zeros((1,len(self._fullens)), dtype=int) * np.NaN
            sub_mems[0, :len(self.getSubMembers())-1] = np.asarray(self.getSubMembers())[self.getSubMembers() != self._truth_member]
            assert(self.getSubMembers()[0] == self._truth_member)
            truth_rval = rvals_dat[rstring_to_rindex[S.getRString()]][truthmem_ind]
            # Make sure to ignore truth member in error calculations
            # Truth member inclusion in subset gives the subset unfair advantage
            subRMS = np.zeros((self._subsize-1))
            fensRMS = np.zeros_like(self._fullens[self._fullens != self._truth_member])
            print(subRMS.shape, fensRMS.shape)
            print(truth_rval)
            print("CALCULATING SUBSET MEMBER RMSE'S...")
            print(self.getSubMembers())
            print(rstring_to_rindex[S.getRString()])
            responses = rvals_dat.variables[rstring_to_rindex[S.getRString()]][:]
            # Ignore truth member in spread calculations too
            sub_rvals_spread = np.var(responses[self.getSubMembers()[self.getSubMembers() != self._truth_member]-1])
            fens_rvals_spread = np.var(responses[~np.in1d(range(len(responses)),truthmem_ind)])
            # print(responses[truthmem_ind], responses[~np.in1d(range(len(responses)),truthmem_ind)])
            for ind, submem in enumerate(self.getSubMembers()
                                [self.getSubMembers() != self._truth_member]):
                rval = responses[submem-1]
                print("Member:", submem, "Response Val:", rval)
                subRMS[ind] = rmse(predictions=rval, targets=truth_rval)
            print("Subset Mean RMSE:", np.mean(subRMS))
            print("CALCULATING FULL ENS MEMBER RMSE'S...")
            for ind, fensmem in enumerate(self._fullens[self._fullens != self._truth_member]):
                rval = responses[fensmem-1]
                print("Member:", fensmem, "Response Val:", rval)
                fensRMS[ind] = rmse(predictions=rval, targets=truth_rval)
            print("Full Ens Mean RMSE:", np.mean(fensRMS))

            ds = xr.Dataset({'Sens_Time': (['subset'], np.atleast_1d(S.getSensTime())),
                     'Subset_Size': (['subset'], np.atleast_1d(self._subsize-1)),
                     'Sens_Vars':  (['subset', 'sensvar'],
                                    np.atleast_2d(sens_vars)),
                     'Subset_Method': (['subset'],
                        np.atleast_1d(list(self._methodchoices.keys())[self._method-1])),
                     'Sens_Threshold': (['subset'], np.atleast_1d(self._percent)),
                     'Response_Func': (['subset'], np.atleast_1d(S.getRString())),
                     'Response_Time': (['subset'], np.atleast_1d(S.getRTime())),
                     'Response_Box': (['subset', 'rbox'],
                                        np.atleast_2d(S.getRbox())),
                     'Full_Ens_RMSE_Rbox': (['subset'], np.atleast_1d(np.mean(fensRMS))),
                     'Full_Ens_Response_Spread': (['subset'], np.atleast_1d(fens_rvals_spread)),
                     'Subset_RMSE_Rbox': (['subset'], np.atleast_1d(np.mean(subRMS))),
                     'Subset_Response_Spread': (['subset'], np.atleast_1d(sub_rvals_spread)),
                     'RMSE_Diff': (['subset'], np.atleast_1d(np.mean(subRMS)-np.mean(fensRMS))),
                     'Neighborhood': (['subset'], np.atleast_1d(self._nbr)),
                     'Response_Thresh': (['subset'], np.atleast_1d(self._thresh)),
                     'Subset_Members': (['subset', 'submems'],
                                        sub_mems)},
                     coords={'run': S.getRunInit(),
                             # 'subset': np.arange(1),
                             'rbox': ['llon', 'ulon', 'llat', 'ulat']})
            if os.path.exists(outpath):
                og_ds = xr.open_dataset(outpath)
                print("Opened original dataset...")
                ds = xr.concat([og_ds, ds], dim='subset', data_vars='minimal')
                og_ds.close()
                print("Appended new dataset to original dataset...\n\n", ds)
            print("FINAL DATASET TO WRITE:", ds)
            ds.to_netcdf("test.nc", unlimited_dims=['subset', 'rthresh'],
                        format="NETCDF4", mode='w')
            os.rename("test.nc", outpath)
            ds.close()
        else:
            print("ERROR! Attempt to run idealized experiment with a real subset obj")
            raise

        return

    def oberrorAddedIdealizedStoreRfuncStatsNetCDF(self, outpath,
                                                rvalspath="Rvals.nc"):
        """
        Stores response function verification stats for an idealized experiment
        in which a randomly selected ensemble member serves as truth and is
        used to compute response errors directly. This version is meant to run
        with ob error added to the early analysis.

        Inputs
        ------
        outpath ----------- absolute output filepath as a string
        rvalspath --------- name of netCDF file produced from
                            sixhresens that stores member response values

        Outputs
        -------
        returns NULL and stores UH verification stats to outpath as netCDF4
        """
        if self._idealized:
            S = self.getSens()
            rstring_to_rindex = {"6-hr Max UH": "UH_MAX",
                                 "6-hr UH Coverage": "UH_COV",
                                 "6-hr Max Refl": "DBZ_MAX",
                                 "6-hr Refl Coverage": "DBZ_COV"}
            truthmem_ind = self._truth_member - 1
            # Open member response values
            rvals_dat = Dataset(S.getDir() + rvalspath)
            # Create entry as xarray dataset
            sens_vars = np.empty((1, 30), dtype="U30")
            sens_vars[0, :len(self._sensvars)] = np.asarray(self._sensvars[:])
            sub_mems = np.zeros((1,len(self._fullens)), dtype=int) * np.NaN
            sub_mems[0, :len(self.getSubMembers())-1] = np.asarray(self.getSubMembers())[self.getSubMembers() != self._truth_member]
            # assert(self.getSubMembers()[0] == self._truth_member)
            truth_rval = rvals_dat[rstring_to_rindex[S.getRString()]][truthmem_ind]
            # Make sure to ignore truth member in error calculations
            # Truth member inclusion in subset gives the subset unfair advantage
            subRMS = np.zeros((self._subsize-1))
            fensRMS = np.zeros_like(self._fullens[self._fullens != self._truth_member])
            print(subRMS.shape, fensRMS.shape)
            print(truth_rval)
            print("CALCULATING SUBSET MEMBER RMSE'S...")
            print(self.getSubMembers())
            print(rstring_to_rindex[S.getRString()])
            responses = rvals_dat.variables[rstring_to_rindex[S.getRString()]][:]
            # Ignore truth member in spread calculations too
            sub_rvals_spread = np.var(responses[self.getSubMembers()[self.getSubMembers() != self._truth_member]-1])
            fens_rvals_spread = np.var(responses[~np.in1d(range(len(responses)),truthmem_ind)])
            # print(responses[truthmem_ind], responses[~np.in1d(range(len(responses)),truthmem_ind)])
            for ind, submem in enumerate(self.getSubMembers()
                                [self.getSubMembers() != self._truth_member]):
                rval = responses[submem-1]
                print("Member:", submem, "Response Val:", rval)
                subRMS[ind] = rmse(predictions=rval, targets=truth_rval)
            print("Subset Mean RMSE:", np.mean(subRMS))
            print("CALCULATING FULL ENS MEMBER RMSE'S...")
            for ind, fensmem in enumerate(self._fullens[self._fullens != self._truth_member]):
                rval = responses[fensmem-1]
                print("Member:", fensmem, "Response Val:", rval)
                fensRMS[ind] = rmse(predictions=rval, targets=truth_rval)
            print("Full Ens Mean RMSE:", np.mean(fensRMS))

            ds = xr.Dataset({'Sens_Time': (['subset'], np.atleast_1d(S.getSensTime())),
                     'Subset_Size': (['subset'], np.atleast_1d(self._subsize)),
                     'Sens_Vars':  (['subset', 'sensvar'],
                                    np.atleast_2d(sens_vars)),
                     'Subset_Method': (['subset'],
                        np.atleast_1d(list(self._methodchoices.keys())[self._method-1])),
                     'Sens_Threshold': (['subset'], np.atleast_1d(self._percent)),
                     'Response_Func': (['subset'], np.atleast_1d(S.getRString())),
                     'Response_Time': (['subset'], np.atleast_1d(S.getRTime())),
                     'Response_Box': (['subset', 'rbox'],
                                        np.atleast_2d(S.getRbox())),
                     'Full_Ens_RMSE_Rbox': (['subset'], np.atleast_1d(np.mean(fensRMS))),
                     'Full_Ens_Response_Spread': (['subset'], np.atleast_1d(fens_rvals_spread)),
                     'Subset_RMSE_Rbox': (['subset'], np.atleast_1d(np.mean(subRMS))),
                     'Subset_Response_Spread': (['subset'], np.atleast_1d(sub_rvals_spread)),
                     'RMSE_Diff': (['subset'], np.atleast_1d(np.mean(subRMS)-np.mean(fensRMS))),
                     'Neighborhood': (['subset'], np.atleast_1d(self._nbr)),
                     'Response_Thresh': (['subset'], np.atleast_1d(self._thresh)),
                     'Subset_Members': (['subset', 'submems'],
                                        sub_mems)},
                     coords={'run': S.getRunInit(),
                             # 'subset': np.arange(1),
                             'rbox': ['llon', 'ulon', 'llat', 'ulat']})
            if os.path.exists(outpath):
                og_ds = xr.open_dataset(outpath)
                print("Opened original dataset...")
                ds = xr.concat([og_ds, ds], dim='subset', data_vars='minimal')
                og_ds.close()
                print("Appended new dataset to original dataset...\n\n", ds)
            print("FINAL DATASET TO WRITE:", ds)
            ds.to_netcdf("test.nc", unlimited_dims=['subset', 'rthresh'],
                        format="NETCDF4", mode='w')
            os.rename("test.nc", outpath)
            ds.close()
        else:
            print("ERROR! Attempt to run idealized experiment with a real subset obj")
            raise

        return

    def semiIdealizedStoreUHStatsNetCDF(self, outpath, pperfpath,
                                    reliabilityobpath, bssobpath, spc_grid,
                                    spc_sigma=None):
        """
        Stores UH verification stats for a semi-idealized experiment in
        which a randomly-chosen ensemble member serves as "truth" and
        generates a surrogate severe field. Stats stored to to netCDF4
        file using xarray.

        Inputs
        ------
        outpath ----------- absolute output filepath as a string
        pperfpath --------- absolute filepath for practically perfect
                            gridded data
        reliabilityobpath - absolute filepath for reliability observational
                            point data on ensemble grid
        spc_grid ---------- boolean specifying whether to use the SPC 211
                            grid and interpolate to WRF, or generate SSRs
                            on WRF grid throughout
        spc_sigma --------- if using SPC grid, may want to define
                            standard deviation of Gaussian kernel
                            (smoothing coefficient) explicitly.
                            Otherwise, will use subset neighborhood
                            and SPC grid-spacing

        Outputs
        -------
        returns NULL and stores UH verification stats to outpath as netCDF4
        """
        ####################### Checking for success ###########################
        def checkSuccess():
            """
            Checks to see if the word "SUCCESSFUL"
            appears in the reliabilitycalc log file, indicating
            the reliability calculation was completed
            correctly.
            """
            with open(S.getDir()+'reliability.out') as out:
                lowcase_out = [o.lower() for o in out]
                success = (np.array(['successful' in s for s in lowcase_out]).any())
            return success
        ########################################################################

        if self._semi_idealized:
            S = self.getSens()
            rstring_to_rindex = {"6-hr Max UH": "UH_MAX",
                                 "6-hr UH Coverage": "UH_COV",
                                 "6-hr Max Refl": "DBZ_MAX",
                                 "6-hr Refl Coverage": "DBZ_COV"}
            truthmem_ind = self._truth_member - 1

            # Create entry as xarray dataset
            sens_vars = np.empty((1, 30), dtype="U30")
            sens_vars[0, :len(self._sensvars)] = np.asarray(self._sensvars[:])
            sub_mems = np.zeros((1,len(self._fullens)), dtype=int) * np.NaN
            sub_mems[0, :len(self.getSubMembers())-1] = np.asarray(
                        self.getSubMembers())[self.getSubMembers() != self._truth_member]
            assert(self.getSubMembers()[0] == self._truth_member)

            # Generate surrogate severe fields
            wrfref = xr.open_dataset(S.getRefFileD2())
            wrflats, wrflons =  wrfref["XLAT"][0].values, wrfref["XLONG"][0].values
            zlevs = len(wrfref["P_HYD"][0,:,0,0])
            if spc_grid:
                # Get lats and lons for practically perfect grid
                ppfile = '{}/../calc/pperf_grid_template.npz'.format(package_dir)
                f = np.load(ppfile)
                lons = f["lon"]
                lats = f["lat"]
                f.close()
                dx = 81. # Approximate horizontal grid spacing of SPC 211 grid
                if spc_sigma is not None:
                    sigma = spc_sigma
                else:
                    sigma = self._nbr / dx
            else:
                lats, lons = wrfref["XLAT"][0].values, wrfref["XLONG"][0].values
                dx = self.getHorizGridSpacingD2()
                sigma = self._nbr / dx
            wrfref.close()
            print(np.shape(lats))
            fens_SSRs = np.zeros((len(self._fullens)-1, len(lats),
                                len(lats[0,:])))
            fens_SSPFs = np.zeros_like(fens_SSRs)
            sub_SSRs = np.zeros((len(self.getSubMembers())-1,
                                len(lats), len(lats[0,:])))
            sub_SSPFs = np.zeros_like(sub_SSRs)
            # Determine time range for SSRs
            if S.getSixHour():
                trange = np.arange(S.getRTime()-5, S.getRTime()+1)
            else:
                trange = np.arange(S.getRTime(), S.getRTime()+1)

            # Define prob and reliability paths
            subprobpath = S.getDir() + "probs/" + self._subprob
            fensprobpath = S.getDir() + "probs/" + self._fensprob
            outrel = S.getDir() + "reliability_out.nc"
            rtimedate = S.getRunInit() + timedelta(hours=S.getRTime())

            # Calculate FSS from idealized pperf
            if os.path.exists(pperfpath):
                fens_fss = FSSnetcdf(fensprobpath, pperfpath, rtimedate,
                            var='updraft_helicity', thresh=self._thresh,
                            rboxpath=S.getDir()+'esens.in')
                sub_fss = FSSnetcdf(subprobpath, pperfpath, rtimedate,
                            var='updraft_helicity', thresh=self._thresh,
                            rboxpath=S.getDir()+'esens.in')
                f_fss_tot, f_fss_rbox, sig = fens_fss
                s_fss_tot, s_fss_rbox, sig = sub_fss
                # Process the reliability observation file by calculating
                #  reliabilty in fortran - those results are stored to
                #  a netCDF file called 'reliability_out.nc'
                print("Processing reliability of {} with Fortran...".format(subprobpath))
                success = False
                i = 0
                while (success == False):
                    storeReliabilityRboxFortran(S.getDir(), S.getRTime(),
                                subprobpath, reliabilityobpath, outrel,
                                rboxpath=S.getDir()+'esens.in',
                                sixhour=S.getSixHour(),
                                variable="updraft_helicity",
                                rthresh=self._thresh, nbrhd=self._nbr,
                                wrfrefpath=S.getRefFileD2())
                    time.sleep(.300)
                    print("{}th attempt of calculating reliability...".format(i))
                    success = checkSuccess()
                    i += 1
                # Process result with xarray
                sub_reliability = xr.open_dataset(outrel)
                prob_bins = sub_reliability["prob_bins"]
                s_fcstfreq_rbox = sub_reliability["fcst_frequency"]
                s_ob_hr_rbox = sub_reliability["ob_hit_rate"]
                sub_reliability.close()
                ####################################################################
                # Calculate BSS using probabilistic forecast and binary observations
                # Same as FSS except binary obs
                ####################################################################
                fens_prob_ds = xr.open_dataset(fensprobpath)
                sub_prob_ds = xr.open_dataset(subprobpath)
                # Choose correct indices based on variable and threshold
                probinds = {'reflectivity': {40: 0, 50: 5},
                            'updraft_helicity': {25: 1, 40: 2, 100: 3},
                            'wind_speed': {40: 4}}
                index = probinds['updraft_helicity'][self._thresh]
                # Pull probabilities
                fens_probs = np.asarray(fens_prob_ds["P_HYD"])[0, index]
                sub_probs = np.asarray(sub_prob_ds["P_HYD"])[0, index]
                fens_prob_ds.close(); sub_prob_ds.close()
                # Open binary observation dataset
                bss_ob_ds = xr.open_dataset(bssobpath)
                bss_obs = np.asarray(bss_ob_ds["nearest_neighbor"][0])
                bss_ob_ds.close()
                # Mask obs/fcst
                ref_ds = xr.open_dataset(S.getRefFileD2())
                lons = np.asarray(ref_ds["XLONG"].values[0], dtype=float)
                lats = np.asarray(ref_ds["XLAT"].values[0], dtype=float)
                ref_ds.close()
                llon, ulon, llat, ulat = np.array(S.getRbox(), dtype=float)
                lonmask = (lons > llon) & (lons < ulon)
                latmask = (lats > llat) & (lats < ulat)
                mask = lonmask & latmask
                del(lons); del(lats)
                masked_fens_probs = np.ma.masked_array(fens_probs/100., mask=~mask)
                masked_sub_probs = np.ma.masked_array(sub_probs/100., mask=~mask)
                masked_obs = np.ma.masked_array(bss_obs, mask=~(mask & (bss_obs <= 9e36)))
                # print("NANMAX BSS OBS", np.nanmax(masked_obs))
                # print("Rbox probs", masked_fens_probs)
                # Calculate BSS!
                f_bss_rbox = FSS(masked_fens_probs, masked_obs)
                print(f"Full Ens BSS: {f_bss_rbox}")
                s_bss_rbox = FSS(masked_sub_probs, masked_obs)
                print(f"Subset BSS: {s_bss_rbox}")
                ####################################################################
            else:
                raise FileNotFoundError('Please run storePracPerf() or \
                set correct obpath.')
            # Create entry as xarray dataset
            sens_vars = np.empty((1, 30), dtype="U30")
            sens_vars[0, :len(self._sensvars)] = np.asarray(self._sensvars[:])
            sub_mems = np.zeros((1,len(self._fullens)), dtype=int) * np.NaN
            sub_mems[0, :len(self.getSubMembers())] = np.asarray(self.getSubMembers())
            ds = xr.Dataset({'Sens_Time': (['subset'], np.atleast_1d(S.getSensTime())),
                     'Subset_Size': (['subset'], np.atleast_1d(self._subsize)),
                     'Sens_Vars':  (['subset', 'sensvar'],
                                    np.atleast_2d(sens_vars)),
                     'Subset_Method': (['subset'],
                        np.atleast_1d(list(self._methodchoices.keys())[self._method-1])),
                     'Sens_Threshold': (['subset'], np.atleast_1d(self._percent)),
                     'Response_Func': (['subset'], np.atleast_1d(S.getRString())),
                     'Response_Time': (['subset'], np.atleast_1d(S.getRTime())),
                     'Response_Box': (['subset', 'rbox'],
                                        np.atleast_2d(S.getRbox())),
                     'Full_Ens_FSS_Total': (['subset'], np.atleast_1d(f_fss_tot)),
                     'Full_Ens_FSS_Rbox': (['subset'], np.atleast_1d(f_fss_rbox)),
                     'Subset_FSS_Total': (['subset'], np.atleast_1d(s_fss_tot)),
                     'Subset_FSS_Rbox': (['subset'], np.atleast_1d(s_fss_rbox)),
                     'Full_Ens_BSS_Rbox': (['subset'], np.atleast_1d(f_bss_rbox)),
                     'Subset_BSS_Rbox': (['subset'], np.atleast_1d(s_bss_rbox)),
                     'Prac_Perf_Sigma': (['subset'], np.atleast_1d(sig[0])),
                     'Neighborhood': (['subset'], np.atleast_1d(self._nbr)),
                     'Response_Thresh': (['subset'], np.atleast_1d(self._thresh)),
                     'Subset_Fcst_Freq_Rbox': (['subset', 'prob_bins'],
                                                np.atleast_2d(s_fcstfreq_rbox)),
                     'Subset_Ob_Hit_Rate_Rbox': (['subset', 'prob_bins'],
                                                 np.atleast_2d(s_ob_hr_rbox)),
                     'Subset_Members': (['subset', 'submems'],
                                        sub_mems)},
                     coords={'run': S.getRunInit(),
                             # 'subset': np.arange(1),
                             'bins': prob_bins,
                             'rbox': ['llon', 'ulon', 'llat', 'ulat']})
            if os.path.exists(outpath):
                og_ds = xr.open_dataset(outpath)
                print("Opened original dataset...")
                # Check first to see if we need to add full ens reliability
                print("Current Rthresh",
                    self._thresh, "Current Full Ens Thresh Vals",
                    og_ds.Full_Ens_Response_Thresh.values)
                if (self._thresh not in og_ds.Full_Ens_Response_Thresh.values):
                    print("Processing full ens reliability for new threshold...")
                    success = False
                    while (success == False):
                        storeReliabilityRboxFortran(S.getDir(), S.getRTime(),
                                    fensprobpath, reliabilityobpath, outrel,
                                    rboxpath=S.getDir()+'esens.in',
                                    sixhour=S.getSixHour(),
                                    variable="updraft_helicity",
                                    rthresh=self._thresh, nbrhd=self._nbr,
                                    wrfrefpath=S.getRefFileD2())
                        time.sleep(.300)
                        success = checkSuccess()
                    fens_reliability = xr.open_dataset(outrel)
                    prob_bins = fens_reliability["prob_bins"]
                    f_fcstfreq_rbox = fens_reliability["fcst_frequency"]
                    f_ob_hr_rbox = fens_reliability["ob_hit_rate"]
                    print("New threshold full ob hit rate", f_ob_hr_rbox)
                    new_ds = xr.Dataset({'Full_Ens_Response_Thresh': (['rthresh'],
                                                np.atleast_1d(self._thresh)),
                                        'Full_Ens_Fcst_Freq_Rbox': (['rthresh', 'prob_bins'],
                                                np.atleast_2d(f_fcstfreq_rbox)),
                                        'Full_Ens_Ob_Hit_Rate_Rbox': (['rthresh', 'prob_bins'],
                                                np.atleast_2d(f_ob_hr_rbox))},
                                        coords={'run': S.getRunInit(),
                                                'bins': prob_bins,
                                                'rbox': ['llon', 'ulon',
                                                        'llat', 'ulat']})
                    print(new_ds)
                    fens_reliability.close()
                    og_fens_ds = xr.concat([og_ds, new_ds], dim='rthresh',
                                    data_vars=['Full_Ens_Response_Thresh',
                                    'Full_Ens_Fcst_Freq_Rbox',
                                    'Full_Ens_Ob_Hit_Rate_Rbox'])
                    new_ds.close()
                    print("Concatenated rel to appending dataset...\n", og_fens_ds)
                    # Finally concatenate with existing netCDF dataset
                    ds = xr.concat([og_fens_ds, ds], dim='subset', data_vars='minimal')
                    og_fens_ds.close()
                else:
                    ds = xr.concat([og_ds, ds], dim='subset', data_vars='minimal')
                og_ds.close()
                print("Appended new dataset to original dataset...\n\n", ds)
                # os.popen("rm {}".format(outpath))
            else:
                print("Adding full ens rel to ", ds)
                success = False
                while (success == False):
                    storeReliabilityRboxFortran(S.getDir(), S.getRTime(),
                                fensprobpath, reliabilityobpath, outrel,
                                rboxpath=S.getDir()+'esens.in',
                                sixhour=S.getSixHour(),
                                variable="updraft_helicity",
                                rthresh=self._thresh, nbrhd=self._nbr,
                                wrfrefpath=S.getRefFileD2())
                    time.sleep(.300)
                    success = checkSuccess()
                fens_reliability = xr.open_dataset(outrel)
                prob_bins = fens_reliability["prob_bins"]
                f_fcstfreq_rbox = fens_reliability["fcst_frequency"]
                f_ob_hr_rbox = fens_reliability["ob_hit_rate"]
                fens_reliability.close()
                ds["Full_Ens_Response_Thresh"] = (('rthresh'),
                                            np.atleast_1d(self._thresh))
                ds["Full_Ens_Fcst_Freq_Rbox"] = (('rthresh', 'prob_bins'),
                                            np.atleast_2d(f_fcstfreq_rbox))
                ds["Full_Ens_Ob_Hit_Rate_Rbox"] = (('rthresh', 'prob_bins'),
                                            np.atleast_2d(f_ob_hr_rbox))
            print("FINAL DATASET TO WRITE:", ds)
            ds.to_netcdf("test.nc", unlimited_dims=['subset', 'rthresh'],
                        format="NETCDF4", mode='w')
            os.rename("test.nc", outpath)
            ds.close()

            return

    def store_uhcov_mae_w_lsrs(self, obpath, outpath):
        """
        Calculates and stores mean absolute error of UH coverage using
        number of gridded (onto native WRF grid) storm reports within
        response box as the observation.

        Inputs
        ------
        obpath ---- absolute file path containing nearest_neighbor variable
                    with gridded binary observations
        outpath --- absolute file path to store mean absolute error result

        Outputs
        -------
        Returns NULL, but stores MAE result to a csv file
        """
        # Get sens object and double-check response function is correct
        S = self.getSens()

        # Define full ensemble and subset members
        # Ensure subset has already been generated
        assert(len(self.getSubMembers()) > 0)
        sub_members = self.getSubMembers() # subset mems
        fens_members = self.getFullEns() # full ensemble mems

        wrfref_d2 = xr.open_dataset(S.getRefFileD2())
        lons = wrfref_d2["XLONG"][0]
        lats = wrfref_d2["XLAT"][0]
        print("Lat Shape:", np.shape(lats))
        wrfref_d2.close()

        # Get response box bounds for masking LSR grid
        llon, ulon, llat, ulat = S.getRbox()
        lonmask = (lons > llon) & (lons < ulon)
        latmask = (lats > llat) & (lats < ulat)
        mask = lonmask & latmask

        # Grab SPC reports for valid time frame
        lsr_grid = xr.open_dataset(obpath)["nearest_neighbor"].values[0]
        lsr_rbox = lsr_grid[mask]
        print("LSR Grid Rbox Shape:", np.shape(lsr_rbox))
        lsr_count_rbox = np.sum(lsr_rbox)
        print("LSR Count within Rbox:", lsr_count_rbox)

        # Pull UH Coverage values
        rvals = xr.open_dataset(S.getDir() + "Rvals.nc")
        forecasts = np.asarray(rvals["UH_COV"])
        full_ens_fcsts = forecasts[fens_members-1]
        print("Full ensemble coverage values:", full_ens_fcsts)

        # Calculate mean absolute error values
        fens_mae = calc.mae(predictions=full_ens_fcsts,
                            targets=np.ones_like(full_ens_fcsts)*lsr_count_rbox)
        print("Full ensemble MAE:", fens_mae)

        sub_fcsts = forecasts[sub_members - 1]
        sub_mae = calc.mae(predictions=sub_fcsts,
                                targets=np.ones_like(sub_fcsts)*lsr_count_rbox)
        print("Subset MAE:", sub_mae)
        rvals.close()

        # Create entry as xarray dataset
        sens_vars = np.empty((1, 30), dtype="U30")
        sens_vars[0, :len(self._sensvars)] = np.asarray(self._sensvars[:])
        sub_mems = np.zeros((1,len(self._fullens)), dtype=int) * np.NaN
        sub_mems[0, :len(self.getSubMembers())] = np.asarray(self.getSubMembers())
        sub_cov = np.zeros((1,len(self._fullens)), dtype=int) * np.NaN
        sub_cov[0, :len(self.getSubMembers())] = sub_fcsts
        ds = xr.Dataset({'Sens_Time': (['subset'], np.atleast_1d(S.getSensTime())),
                 'Subset_Size': (['subset'], np.atleast_1d(self._subsize)),
                 'Analysis': (['subset'], np.atleast_1d(self._analysis_type)),
                 'Sens_Vars':  (['subset', 'sensvar'],
                                np.atleast_2d(sens_vars)),
                 'Subset_Method': (['subset'],
                    np.atleast_1d(list(self._methodchoices.keys())[self._method-1])),
                 'Sens_Threshold': (['subset'], np.atleast_1d(self._percent)),
                 'Response_Func': (['subset'], np.atleast_1d(S.getRString())),
                 'Response_Time': (['subset'], np.atleast_1d(S.getRTime())),
                 'Response_Box': (['subset', 'rbox'],
                                    np.atleast_2d(S.getRbox())),
                 'LSR_Count_Rbox': (['subset'], np.atleast_1d(lsr_count_rbox)),
                 'Subset_UH_Coverage_Rbox': (['subset', 'submems'],
                                            np.atleast_2d(sub_cov)),
                 'Full_Ens_MAE_LSR_Rbox': (['subset'], np.atleast_1d(fens_mae)),
                 'Subset_MAE_LSR_Rbox': (['subset'], np.atleast_1d(sub_mae)),
                 'Neighborhood': (['subset'], np.atleast_1d(self._nbr)),
                 'Response_Thresh': (['subset'], np.atleast_1d(self._thresh)),
                 'Subset_Members': (['subset', 'submems'],
                                    sub_mems)},
                 coords={'run': S.getRunInit(),
                         # 'subset': np.arange(1),
                         'rbox': ['llon', 'ulon', 'llat', 'ulat']})
        if os.path.exists(outpath):
            og_ds = xr.open_dataset(outpath)
            print("Opened original dataset...")
            # Check first to see if we need to add full ens reliability
            print("Current Rthresh",
                self._thresh, "Current Full Ens Thresh Vals",
                og_ds.Full_Ens_Response_Thresh.values)
            if (self._thresh not in og_ds.Full_Ens_Response_Thresh.values):
                new_ds = xr.Dataset({'Full_Ens_Response_Thresh': (['rthresh'],
                                            np.atleast_1d(self._thresh))},
                                    coords={'run': S.getRunInit(),
                                            'rbox': ['llon', 'ulon', 'llat', 'ulat']})
                og_fens_ds = xr.concat([og_ds, new_ds], dim='rthresh',
                                data_vars=['Full_Ens_Response_Thresh'])
                new_ds.close()
                print("Concatenated rel to appending dataset...\n", og_fens_ds)
                # Finally concatenate with existing netCDF dataset
                ds = xr.concat([og_fens_ds, ds], dim='subset', data_vars='minimal')
                og_fens_ds.close()
            else:
                ds = xr.concat([og_ds, ds], dim='subset', data_vars='minimal')
            og_ds.close()
            print("Appended new dataset to original dataset...\n\n", ds)
            # os.popen("rm {}".format(outpath))
        else:
            ds["Full_Ens_Response_Thresh"] = (('rthresh'),
                                        np.atleast_1d(self._thresh))
        print("FINAL DATASET TO WRITE:", ds)
        ds.to_netcdf("test.nc", unlimited_dims=['subset', 'rthresh'],
                    format="NETCDF4", mode='w')
        os.rename("test.nc", outpath)
        ds.close()

        return
