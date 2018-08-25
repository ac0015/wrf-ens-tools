#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  8 10:11:46 2018

@author: aucolema
"""
import numpy as np
import os
from datetime import timedelta
from sens import Sens
from esens_subsetting import ensSubset
from interp_analysis import fromDatetime, interpRAPtoWRF, subprocess_cmd
import research_plotting 
from netCDF4 import Dataset
from calc import FSS, Reliability
from post_process import storePracPerf, gen_dict, process_wrf
    
#################
# Subset class
#################
class Subset:
    '''
    The Subset class encapsulates all the information pertaining
    to a forecast ensemble and its subset which is gathered from
    relevant sensitivity variables and response functions and their
    respective sensitivity fields.
    '''    
    def __init__(self, sens=Sens(), subset_size=21, subset_method='percent',
                 percent=70., sensvalfile="SENSvals.nc", nbrhd=32.1869, thresh=25.,
                 sensvars=['500_hPa_GPH', '700_hPa_T', '850_hPa_T', 'SLP'],
                 analysis_type="RAP"):
        '''
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
                        for probabilities in verification statistics.
        sensvars ------ list of strings containing the desired sensitivity variables
                        on which to base the subset. Must comply with naming conventions
                        of keys in the RAP interpolation from interp_analysis.py.
        '''
        self._sens = sens
        self._fullens = np.arange(1,sens.getEnsnum()+1,1)
        self._subsize = subset_size
        self._methodchoices = {'point' : 1, 'weight' : 2, 'percent' : 3}
        self._percent = percent
        self._sensvars = sensvars
        self._thresh = thresh
        self._nbr = nbrhd
        
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
        self._analysis_type = analysis_type
        if self._analysis_type.upper() == "RAP":
            self._analysis = "{}RAP_interp_to_WRF_{}{}{}{}.nc".format(sens.getDir(),
                          yr, mo, day, hr)
        elif self._analysis_type.upper() == "WRF":
            self._analysis = "{}WRF_analysis_{}{}{}{}.nc".format(sens.getDir(),
                          yr, mo, day, hr)
        else:
            raise ValueError('{} not a vaid option. RAP and WRF are supported analyisis types.'.format(self._analysis_type))
        
        self._sensvalfile = sens.getDir() + sensvalfile
        self._fensprob = "FULLENSwrfout_nbr{}_f{}.prob".format(int(self._nbr), self._sens.getRTime())
        self._subprob = "SUBSETwrfout_nbr{}_f{}.prob".format(int(self._nbr), self._sens.getRTime())
        
        # Set to True if calcProbs was last run for a subset.
        self._calcprobs_subset = False
        
    def __str__(self):
        return "Subset object with full ensemble of {} members, subset size of {}, and using the {} subsetting method \
        with a threshold of {}, neighborhood of {}, response threshold of {} and these sensitivity variables: {}. \
        Based on Sens object: \n {}".format(self._sens.getEnsnum(), self._subsize,
        self._method, str(self._percent), str(self._nbr), str(self._thresh), ','.join(self.getSensVars()), str(self.getSens()))
        
    def setSubsetMethod(self, subset_method):
        '''
        Set the ensemble subsetting method
        '''
        self._method = self._methodchoices[subset_method]
        return
    
    def setSubsetSize(self, subset_size):
        '''
        Set the ensemble subset size.
        '''
        self._subsize = subset_size
        return
    
    def setAnalysis(self, analysispath, analysistype):
        '''
        Set the absolute path of the interpolated analysis file
        to be used for verification and set the analysis type.
        Options for analysis type are RAP or WRF. Make sure the
        analysispath is of the same type described in analysistype.
        '''
        self._analysis = analysispath
        self._analysis_type = analysistype
          
    def getAnalysis(self):
        '''
        Returns path to analysis file at senstime.
        '''
        return self._analysis
    
    def getSubMembers(self):
        '''
        Returns a list of the members in the subset.
        '''
        return self._subset
    
    def getSens(self):
        '''
        Returns Sens object valid for subset.
        '''
        return self._sens
    
    def getSensVars(self):
        '''
        Returns sensitivity variables being
        used for subsetting.
        '''
        return self._sensvars
    
    def interpRAP(self):
        '''
        If using RAP analysis, must call this function to interpolate
        the analysis to our WRF grid before doing any subsetting. Returns
        NULL but will produce outfile with filepath as described by 
        the analysis attribute of the Subset instance.
        '''
        os.chdir(self.getSens().getDir())
        yr, mo, day, hr = fromDatetime(self._sensdate, interp=True)
        interpRAPtoWRF(yr, mo, day, hr, self.getSens().getRefFileD1()) 
        return 
 
    def processWRFAnalysis(self):
        '''
        If using WRF analysis, must call this function to post-process
        all necessary sensitivity variables for subsetting. Returns NULL
        but will produce netCDF outfile with filepath described by the
        nalysis attribute of the Subset instance.
        '''
        basedir = self.getSens().getDir()
        # Analysis will only be one file without subdirectory, but still
        #  need to give process_wrf a dictionary, so let gen_dict()
        #  take care of that.
        inpath = gen_dict(basedir, 1, subdir="", ntimes=1)
        # Use default vars for process_wrf and default naming conventions.
        process_wrf(inpath, outpath=self._analysis, reduced=True)
        
    def calcSubset(self):
        '''
        Calls the ensSubset() function from esens_subsetting.py. 
        '''
        S = self.getSens()
        if (os.path.isfile(self._analysis) == False):
            if self._analysis_type == "RAP":
                print("{} does not exist. Running RAP interpolation.".format(self._analysis))
                self.interpRAP()
            elif self._analysis_type == "WRF":
                print("{} does not exist. Running WRF post-process.".format(self._analysis))
                self.processWRFAnalysis()
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
    
    def calcProbs(self, members, gridspacing=None):
        '''
        Runs the fortran executable calcprobSUBSET to calculate
        probabilities for any number of ensemble members. Takes
        an input file with ensemble number
        '''
        S = self.getSens()
        
        ### calcProbs-specific error-handling methods. #######
        def checkOverHundred(probpath, sens_obj, args):
            '''
            Prints max probability and returns True if max prob
            is higher than 100.
            '''
            try:
                probs = Dataset(probpath)
                probvar = probs.variables['P_HYD'][0]
                fail = np.max(probvar[:4]) > 100.
                print("Max probability: ", np.max(probvar[:4]))
                probs.close()
            except:
                reRun(probpath, sens_obj, args)
                return checkOverHundred(probpath, sens_obj, args)
            return fail
        
        def checkZero(probpath, sens_obj, args):
            '''
            Prints max probability and returns True is max prob
            is zero.
            '''
            try:
                probs = Dataset(probout)
                probvar = probs.variables['P_HYD'][0]
                possible_fail = np.max(probvar[:4]) == 0.
                #print("Max probability: ", np.max(probvar[:4]))
                probs.close()
            except:
                reRun(probpath, sens_obj, args)
                return checkZero(probpath, sens_obj, args)
            return possible_fail

        def reRun(probpath, sens_obj, args):
            '''
            Re-runs calcProbs.
            '''
            if os.path.exists(probpath):
                os.popen('rm {}'.format(probpath))
            os.popen("cp {} {}".format(sens_obj.getRefFileD2(), probpath))
            subprocess_cmd(args)
            return

        # Create or navigate into probs directory
        probdir = S.getDir() + "probs/"
        direxists = os.path.exists(probdir)
        if (direxists == False):
            os.mkdir(probdir)
        os.chdir(probdir)
        print("Calculating probs for: ", members)
        if len(members) == len(self._fullens): 
            fname = "fullens_probs.in"
            probout = self._fensprob
            self._calcprobs_subset = False
        else: 
            fname = "subset_probs.in"
            probout = self._subprob
            self._calcprobs_subset = True
            
        # Format input file
        os.popen('echo {} > {}'.format(len(members), fname))
        os.popen('echo {} >> {}'.format(' '.join(np.array(members, dtype=str)), fname))
        os.popen('echo {} >> {}'.format(str(S.getRTime()), fname))
        os.popen('echo {} >> {}'.format(str(self._nbr), fname))
        os.popen('echo {} >> {}'.format(probout, fname))
        print("Storing to {}".format(probout))
        
        # Initialize probfile
        if os.path.exists(probout):
            os.popen('rm {}'.format(probout))
        os.popen("cp {} {}".format(S.getRefFileD2(), probout))
        
        # Run fortran executable
        # TO-DO: figure out why this throws a 'not netCDF error'
        # Doesn't throw this error when probcalcis executed directly
        args = "/lustre/work/aucolema/enkfDART/src/probcalc <{} >probs.out".format(fname)
        try:
            subprocess_cmd(args)
        except OSError:
            # Try again - sometimes this executable is finicky
            try:
                reRun(probout, S, args)
            except:
                print('probcalc failed twice. Continue to error checks.')
                
        # Check to make sure it ran correctly
        if os.path.exists(probout) == False:
            print('probcalc failed. Re-running.')
            reRun(probout, S, args)
        fail_thresh = checkOverHundred(probout, S, args)
        
        # If max prob is zero, check probs.out to make sure probcalc didn't
        #  somehow fail without error. (This exectuable is finicky as helll).
        possible_fail = checkZero(probout, S, args)
        
        # If max prob is over hundred, probcalcfailed and needs to re-run.      
        if fail_thresh:
            print('Max prob indicates prob calc did not run correctly')
            print('Attempting once more...')
            reRun(probout, S, args)
            if os.path.exists(probout) == False:
                reRun(probout, S, args)
            fail_thresh = checkOverHundred(probout, S, args)
            possible_fail = checkZero(probout, S, args)
            # Repeat checks one more time.
            if fail_thresh:
                print('Failed again. Re-running entire function.')
                return self.calcProbs(members)
            elif possible_fail:
                with open('probs.out') as out:
                    lowcase_out = [o.lower() for o in out]
                    print(lowcase_out)
                    # If the word 'error' in probs.out re-run, or if max is still zero
                    #  and we're calculating for the full ensemble, likely needs re-run.
                    fail = (np.array(['error' in s for s in lowcase_out]).any()) or (probout == self._fensprob)
                    if fail:
                        print("The rare and elusive netCDF error may have been detected. Re-run entire function")
                        return self.calcProbs(members)
        # If max prob is zero, probcalcSUBEST likely failed. Investigate.
        elif possible_fail:
            with open('probs.out') as out:
                lowcase_out = [o.lower() for o in out]
                print(lowcase_out)
                # If the word 'error' in probs.out re-run, or if max is still zero
                #  and we're calculating for the full ensemble, likely needs re-run.
                fail = (np.array(['error' in s for s in lowcase_out]).any()) or (probout == self._fensprob) or (len(np.array([s for s in lowcase_out])) < 5)
                #print(fail)
                if fail:
                    print("The rare and elusive netCDF error may have been detected. Re-run probcalcSUBSET.")
                    reRun(probout, S, args)
                    fail_thresh = checkOverHundred(probout, S, args)
                    possible_fail = checkZero(probout, S, args)
                    # Repeat checks one more time.
                    if fail_thresh:
                        print('Failed again. Re-running entire function.')
                        return self.calcProbs(members)
                    elif possible_fail:
                        with open('probs.out') as out:
                            lowcase_out = [o.lower() for o in out]
                            print(lowcase_out)
                            # If the word 'error' in probs.out re-run, or if max is still zero
                            #  and we're calculating for the full ensemble, likely needs re-run.
                            fail = (np.array(['error' in s for s in lowcase_out]).any()) or (probout == self._fensprob)
                            if fail:
                                print('Failed again. Re-running entire function.')
                                return self.calcProbs(members)
        return
    
    def researchPlotProbs(self, use_subset):
        '''
        Calls plotProbs from research_plotting library. Passes
        mainly file paths and some metadata and function will
        handle the rest.
        '''
        S = self.getSens()
        direc = S.getDir() + "probs/"
        fullenspath = direc + self._fensprob
        subsetpath = direc + self._subprob
        if use_subset:
            path = subsetpath
        else:
            path = fullenspath
        wrfrefpath = S.getRefFileD2()
 
        research_plotting.plotProbs(path, wrfrefpath, 
                                    S.getRbox(), S.getRTime(), outpath=direc)
        return
            
    def researchPlotDiffs(self, use_subset, verif_day=False):
        '''
        Calculates and plots delta probabilities between the 1-hr full ensemble
        probs and its corresponding 1-hr subset probs. Calls plotDiff from 
        reserch_plotting library. Passes only file paths and outside function 
        will handle the rest. If verif_day is set to True, will overlay 
        full day's storm reports onto difference plot.
        '''
        S = self.getSens()
        direc = S.getDir() + 'probs/'
        fullenspath = direc + self._fensprob
        subsetpath = direc + self._subprob
        wrfrefpath = S.getRefFileD2()
        if use_subset:
            path = subsetpath
        else:
            path = fullenspath
        
        # Format date for storm reports
        date = S.getRunInit()
        SPCdate = str(date)[2:10].replace('-','')
        
        research_plotting.plotDiff(path, wrfrefpath,
                                   S.getRbox(), S.getRTime(), SPCdate,
                                   stormreports=verif_day)
        return
    
    def plotSixPanels(self, storm_reports=True):
        S = self.getSens()
        # Make sure subset has already been calculated
        if self._subset == None:
            print("Identifying subset members for current Subset object...")
            self.calcSubset()
        yr, mo, day, hr = fromDatetime(S.getRunInit(), interp=False)
        rfunclabel = S.getRString().replace(' ','').lower()
        dirdate = str(yr) + str(mo) + str(day) + str(hr)
        research_plotting.plotSixPanels(dirdate, storm_reports, 
                                        self.getSubMembers(), sixhour=False,
                                        time=S.getRTime(), subsettype=rfunclabel)
        return
    
    def storePracPerf(self, pperfpath, sigma=2):
        '''
        Runs storePracPerf() from post_process module
        using time data from the subset obj and stores
        to netcdf path specified with pperfpath. Optionally
        alter sigma to be used in gaussian filter. 
        Defaults to 2 which is operational norm.        
        '''
        S = self.getSens()
        # Calculate pperf for hour before, of, and after response time
        fhrs = [S.getRTime()-1, S.getRTime(), S.getRTime()+1]
        storePracPerf(S.getRunInit(), fhrs, pperfpath, sigma=sigma)
        return
    
    def storeUHStats(self, outpath, pperfpath):
        '''
        Calculates and stores subset info along with
        six verification metrics, which are calculated
        for both the full ensemble and the subset within and
        outside of the response box. Stores in a netCDF
        file with path specified by outpath and uses 
        practically perfect stored in pperfpath for 
        metric calculations.
        '''
        S = self.getSens()
        fensprobpath = S.getDir() + "probs/" + self._fensprob
        subprobpath = S.getDir() + "probs/" + self._subprob
        rtimedate = S.getRunInit() + timedelta(hours=S.getRTime())
        if os.path.exists(outpath):
            statsout = Dataset(outpath, 'a')
        else:
            statsout = Dataset(outpath, 'w')
            statsout.createDimension('Times',  None)
            statsout.createDimension('vars', None)
            statsout.createDimension('rbox', 4)
            # Reliability diagram will need three variables
            statsout.createDimension('rel', 3)
            statsout.createDimension('bins', 10)
            statsout.createVariable('Run_Init', str, ('Times'))
            statsout.createVariable('Sens_Time', int, ('Times'))
            statsout.createVariable('Analysis', str, ('Times'))
            statsout.createVariable('Sens_Vars', str, ('Times', 'vars'))
            statsout.createVariable('Subset_Size', int, ('Times'))
            statsout.createVariable('Subset_Method', str, ('Times'))
            statsout.createVariable('Sens_Threshold', int, ('Times'))
            statsout.createVariable('Response_Func', str, ('Times'))
            statsout.createVariable('Response_Time', int, ('Times'))
            statsout.createVariable('Response_Box', float, ('Times', 'rbox'))
            statsout.createVariable('Full_Ens_FSS_Total', float, ('Times'))
            statsout.createVariable('Full_Ens_FSS_Rbox', float, ('Times'))
            statsout.createVariable('Subset_FSS_Total', float, ('Times'))
            statsout.createVariable('Subset_FSS_Rbox', float, ('Times'))
            statsout.createVariable('Prac_Perf_Sigma', int, ('Times'))
            statsout.createVariable('Neighborhood', float, ('Times'))
            statsout.createVariable('Response_Thresh', float, ('Times'))
            statsout.createVariable('Full_Ens_Reliability_Total', float, ('Times', 'rel', 'bins'))
            statsout.createVariable('Full_Ens_Reliability_Rbox', float, ('Times', 'rel', 'bins'))
            statsout.createVariable('Subset_Reliability_Total', float, ('Times', 'rel', 'bins'))
            statsout.createVariable('Subset_Reliability_Rbox', float, ('Times', 'rel', 'bins'))
        if os.path.exists(pperfpath):
            fens_fss = FSS(fensprobpath, pperfpath, rtimedate, var='updraft_helicity',
                      thresh=self._thresh, rboxpath=S.getDir()+'esens.in')
            sub_fss = FSS(subprobpath, pperfpath, rtimedate, var='updraft_helicity',
                      thresh=self._thresh, rboxpath=S.getDir()+'esens.in')
            fens_reliability = Reliability(fensprobpath, S.getRunInit(), S.getRTime(),
                                           obpath=None, var='updraft_helicity',
                                           thresh=self._thresh, rboxpath=S.getDir()+'esens.in',
                                           sixhr=False)
            sub_reliability = Reliability(subprobpath, S.getRunInit(), S.getRTime(),
                                           obpath=None, var='updraft_helicity',
                                           thresh=self._thresh, rboxpath=S.getDir()+'esens.in',
                                           sixhr=False)
            # prob_bins and ob hit rate variables will stay the same, so OK to clobber
            prob_bins, f_fcstfreq_tot, ob_hr_tot, f_fcstfreq_rbox, ob_hr_rbox = fens_reliability
            prob_bins, s_fcstfreq_tot, ob_hr_tot, s_fcstfreq_rbox, ob_hr_rbox = sub_reliability
            f_fss_tot, f_fss_rbox, sig = fens_fss
            s_fss_tot, s_fss_rbox, sig = sub_fss
            init = statsout.variables['Run_Init']
            senstime = statsout.variables['Sens_Time']
            subsize = statsout.variables['Subset_Size']
            analysis = statsout.variables['Analysis']
            sensvar = statsout.variables['Sens_Vars']
            submethod = statsout.variables['Subset_Method']
            thresh = statsout.variables['Sens_Threshold']
            rfunc = statsout.variables['Response_Func']
            resptime = statsout.variables['Response_Time']
            respbox = statsout.variables['Response_Box']
            fensfsstot = statsout.variables['Full_Ens_FSS_Total']
            fensfssrbox = statsout.variables['Full_Ens_FSS_Rbox']
            subfsstot = statsout.variables['Subset_FSS_Total']
            subfssrbox = statsout.variables['Subset_FSS_Rbox']
            pperfsig = statsout.variables['Prac_Perf_Sigma']
            nbr = statsout.variables['Neighborhood']
            respthresh = statsout.variables['Response_Thresh']
            fens_rel_tot = statsout.variables['Full_Ens_Reliability_Total']
            fens_rel_rbox = statsout.variables['Full_Ens_Reliability_Rbox']
            sub_rel_tot = statsout.variables['Subset_Reliability_Total']
            sub_rel_rbox = statsout.variables['Subset_Reliability_Rbox']
            n = len(init)
            #print(n, str(S.getRunInit()))
            init[n] = str(S.getRunInit())
            senstime[n] = S.getSensTime()
            analysis[n] = self._analysis_type
            sensvar[n,:] = np.atleast_2d(self._sensvars)[:]
            subsize[n] = self._subsize
            submethod[n] = list(self._methodchoices.keys())[self._method-1]
            thresh[n] = self._percent
            rfunc[n] = S.getRString()
            resptime[n] = S.getRTime()
            respbox[n,:] = S.getRbox()[:]
            respthresh[n] = self._thresh
            fensfsstot[n] = f_fss_tot
            fensfssrbox[n] = f_fss_rbox
            subfsstot[n] = s_fss_tot
            subfssrbox[n] = s_fss_rbox
            pperfsig[n] = sig
            nbr[n] = self._nbr
            fens_rel_tot[n,:,:] = np.atleast_2d(np.vstack((prob_bins, f_fcstfreq_tot, ob_hr_tot)))[:]
            fens_rel_rbox[n,:,:] = np.atleast_2d(np.vstack((prob_bins, f_fcstfreq_rbox, ob_hr_rbox)))[:]
            sub_rel_tot[n,:,:] = np.atleast_2d(np.vstack((prob_bins, s_fcstfreq_tot, ob_hr_tot)))[:]
            sub_rel_rbox[n,:,:] = np.atleast_2d(np.vstack((prob_bins, s_fcstfreq_rbox, ob_hr_rbox)))[:]
            statsout.close()
        else:
            raise FileNotFoundError('Please run storePracPerf() or set correct obpath.')
        return
    
        
                                              
        
        
        
        
        
        
        
        
