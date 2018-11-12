#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 10:07:04 2018

@author: aucolema
"""
import numpy as np
import os
from datetime import datetime
from datetime import timedelta
from interp_analysis import subprocess_cmd, fromDatetime

#################
# Sens class
#################
class Sens:
    '''
    The Sens class encapsulates the information pertaining to
    a forecast ensemble, response function, response function
    box, sensitivity time, and response time for running a
    sensitivity code package and producing a sensitivity
    object to be used by other applications.
    '''
    def __init__(self, infile=False, gui=False, ensbasepath=None,
                 llat=31.578, ulat=35.578, llon=-105.855, ulon=-98.855,
                 rfuncstr="Max 1h UH", senstime=6, rtime=24, run=None,
                 submit=datetime.utcnow()):
        '''
        Constructor for an instance of the Sens class. Defaults to
        a response function box centered around Lubbock, TX with
        the Max UH as the response and most recent ensemble to
        calculate sensitivity.
        '''
        # Set run initialization time
        if run == None:
            # If no run is provided, use previous real-time run
            if submit.hour < 12:
                hr, day = 0, -1
            else:
                hr, day = 12, 0
            dayinit = submit + timedelta(days=day)
            self._run = datetime(dayinit.year, dayinit.month,
                                 dayinit.day, hr)
        else:
            self._run = run

        # Base directory containing ensemble run
        yr, mo, day, hr = fromDatetime(self._run)
        if ensbasepath is None:
            # Assume Austin's using this on quanah
            fullpath = "/lustre/research/bancell/aucolema/HWT2016runs/"
            self._dir = fullpath + "{}{}{}{}/".format(yr, mo, day, hr)
        else:
            self._dir = ensbasepath

        # If esens.in or subsetGUI.txt available, use them
        if infile:
            if gui:
                esens_in = np.genfromtxt(self._dir + 'subsetGUI.txt', dtype=str)
                self._rtime = int(esens_in[0])
                self._llon = float(esens_in[1])
                self._ulon = float(esens_in[2])
                self._llat = float(esens_in[3])
                self._ulat = float(esens_in[4])
                self._ensnum = 42
                self._senstime = senstime
                self._rfuncstr = rfuncstr
            else:
                esens_in = np.genfromtxt(self._dir + 'esens.in', dtype=str)
                self._ensnum = int(esens_in[0])
                self._senstime = int(esens_in[1])
                self._rtime = int(esens_in[2])
                self.setRFunc(int(esens_in[3]))
                self._llon = float(esens_in[4])
                self._ulon = float(esens_in[5])
                self._llat = float(esens_in[6])
                self._ulat = float(esens_in[7])
        else:
            self._llat = llat
            self._ulat = ulat
            self._llon = llon
            self._ulon = ulon
            self._rfuncstr = rfuncstr
            self._rtime = int(rtime)
            # TO-DO: Decide how/if these should be chosen
            self._ensnum = 42
            self._senstime = senstime

        self._fhrs = 48
        self._date = submit

        # Default in/outfile paths and reference file paths
        self._sensin = self._dir + "esens.in"
        self._wrfrefd1 = self._dir + "wrfoutREF"
        self._wrfrefd2 = self._dir + "wrfoutREFd2"
        self._wrfsens = self._dir + "wrfout.sens"


    def __str__(self):
        rbox = "Lats: {} to {}, Lons: {} to {}".format(self._llat, self._ulat,
                      self._llon, self._ulon)
        return "Run initialized at: {} \n Response Box: {} \n  \
        Response Function: {} at f{}".format(self._run, rbox,
        self._rfuncstr, self._rtime)

    def setDir(self, dirpath):
        '''
        Set base directory from which sensitivity code
        will run.
        '''
        try:
            if os.path.exists(dirpath):
                self._dir = dirpath
            else:
                raise NameError(dirpath + " does not exist.")
        except:
            raise
        return

    def setInfile(self, absinfilepath):
        '''
        Set the absolute path of the input file supplied to the
        sensitivity code.
        '''
        self._sensin = absinfilepath
        return

    def setRefFileD1(self, absrefpathd1):
        '''
        Set the WRF reference file path for the outer domain.
        '''
        self._wrfrefd1 = absrefpathd1
        return

    def setRefFileD2(self, absrefpathd2):
        '''
        Set the WRF reference file path for the inner domain.
        Also resets grid-spacing parameter.
        '''
        self._wrfrefd2 = absrefpathd2
        return

    def setRFunc(self, index):
        '''
        Set response function with an index.
        '''
        rfuncstrs = {1 : "Avg Refl", 2 : "Max Refl",
                     3 : "Avg 1h UH", 4 : "Max 1h UH",
                     5 : "Accum PCP", 6 : "Avg Wind Spd"}
        self._rfuncstr = rfuncstrs[index]

    def getRIndex(self):
        '''
        Returns the integer corresponding to the response function
        string provided to the constructor for use with the sensitivity
        code. These strings are specific to those used in the subset GUI.
        '''
        rfuncinds = {"Avg Refl" : 1, "Max Refl" : 2,
                     "Avg 1h UH" : 3, "Max 1h UH" : 4,
                     "Accum PCP" : 5, "Avg Wind Spd" : 6}
        return rfuncinds[self._rfuncstr]


    def getRString(self):
        '''
        Returns response function string
        '''
        return self._rfuncstr

    def getRunInit(self):
        '''
        Returns the initialization time of ensemble run as a datetime
        object.
        '''
        return self._run

    def getSensTime(self):
        '''
        Returns the sensitivity time as an integer of n forecast hrs
        from the run's initialization.
        '''
        return self._senstime

    def getRTime(self):
        '''
        Returns the response function time as an integer of n forecast
        hours from the run's initialization.
        '''
        return self._rtime

    def getEnsnum(self):
        '''
        Returns ensemble size as an integer.
        '''
        return self._ensnum

    def getDir(self):
        '''
        Returns base directory for ensemble run as a string.
        '''
        return self._dir

    def getRbox(self):
        '''
        Returns the response function box bounds as a tuple
        containing lower latitude, upper latitude, lower longitude,
        and upper longitude respectively.
        '''
        return self._llon, self._ulon, self._llat, self._ulat

    def getInfile(self):
        '''
        Returns the input file path provided to the sensitivity code.
        '''
        return self._sensin

    def getRefFileD1(self):
        '''
        Returns the WRF outer domain reference file path as a string.
        '''
        return self._wrfrefd1

    def getRefFileD2(self):
        '''
        Returns the WRF inner domain reference file path as a string.
        '''
        return self._wrfrefd2

    def getWRFSensFile(self):
        '''
        Returns the WRF sens file used by the Fortran sensitivity
        code. This is hardcoded in the Fortran, so only change
        if the Fortran naming conventions are changed.
        '''
        return self._wrfsens

    def createInfile(self, rvals=True):
        '''
        Creates (or overwrites if it already exists) an input file for
        the sensitivity code to use. If rvals is true in infile,
        sensitivity code will write response function values for each
        member to netCDF in 'Rvals.nc'.
        '''
        fpath = self._sensin
        if os.path.isfile(fpath):
            os.remove(fpath)
        if rvals:
            rval = ".true."
        else:
            rval = ".false."
        args = [str(self._ensnum), str(self._senstime), str(self._rtime), str(self.getRIndex()),
                str(self._llon), str(self._ulon), str(self._llat), str(self._ulat), rval]
        np.savetxt(fpath, args, fmt="%s", delimiter='\n')
        return

    def renameWRFOUT(self, restored=False):
        '''
        Renames wrfout files for current ensemble run

        Inputs
        ------
        restored - optional boolean specifying if wrfout files
                    are restored to their full capacity

        '''
        for i in range(self._ensnum):
            subdir = self._dir + "mem" + str(i+1) + "/"
            for t in range(self._fhrs + 1):
                current = self._run + timedelta(hours = t)

                # Take care of leading zeros to build path string
                if len(str(current.month)) == 1: mo = "0" + str(current.month)
                else: mo = str(current.month)
                if len(str(current.day)) == 1: day = "0" + str(current.day)
                else: day = str(current.day)
                if len(str(current.hour)) == 1: hr = "0" + str(current.hour)
                else: hr = str(current.hour)

                # Build path string
                domain_strs = ['SENS', 'R']
                dom = 1

                # If restored, files have different naming conventions
                if restored:
                    wrfout_strs = ['{}wrfout_d0{}_red_{}-{}-{}_{}:00:00',
                                   '{}res_wrfout_d0{}_{}-{}-{}_{}:00:00']
                else:
                    wrfout_strs = ['{}wrfout_d0{}_red_{}-{}-{}_{}:00:00',
                                   '{}wrfout_d0{}_red_{}-{}-{}_{}:00:00']
                for d in range(len(domain_strs)):
                    tmp = wrfout_strs[d].format(subdir,
                           str(dom), str(current.year), mo, day, hr)
                    new = '{}{}{}_{}.out'.format(subdir, domain_strs[d], str(i+1), str(t))
                    tmpexists, newexists = os.path.isfile(tmp), os.path.isfile(new)
                    if (tmpexists == True) and (newexists == False):
                        os.rename(tmp, new)
                        print('Renaming {} to {}'.format(tmp, new))
                        #print('yeet')
                    dom += 1
        return

    def runMeanCalc(self):
        '''
        Runs the fortran program 'meancalcSENS' to calculate
        and store means of sensitivity variables for use with
        sensitivity code.
        '''
        os.chdir(self._dir)
        tmpmem = "mem1/SENS1_0.out"
        tmpRmem = "mem1/R1_0.out"
        SENSoutfile = "SENSmean.out"
        Routfile = "Rmean.out"

        try:
            # Need to copy a sensivity member into SENSmean and Rmean outfiles
            sens_mem = os.path.relpath(tmpmem)
            r_mem = os.path.relpath(tmpRmem)
        except:
            print(os.sys.exc_info()[0])
            raise
        else:
            os.popen('cp {} {}'.format(sens_mem, SENSoutfile))
            os.popen('cp {} {}'.format(r_mem, Routfile))

        args = "module load intel; /lustre/work/aucolema/enkfDART/src/meancalcSENS >meancalcSENS.out"
        subprocess_cmd(args)
        return

    def runSENS(self):
        '''
        Runs the fortran sensitivity code and stores sensitivity
        netCDF file as the outfile dictated by the sensitivity object.
        '''
        os.chdir(self._dir)

        try:
            if (os.path.isfile(self._sensin)):
                os.popen('cp {} {}'.format(self._wrfrefd1, self._wrfsens))
            else:
                raise IOError("{} does not exist. Running createInfile()".format(self._sensin))
        except IOError as msg:
            print(msg)
            self.createInfile()
        except:
            print(os.sys.exc_info()[0])
            raise

        args = "module load intel; /lustre/work/aucolema/enkfDART/src/esensSPC <{} >esens.out".format(self._sensin)
        subprocess_cmd(args)
        return

    def storeSENSvals(self):
        '''
        Runs fortran code to store individual member sensitivity variable values
        to 'SENSvals.nc'.
        '''
        os.chdir(self._dir)
        args = "module load intel; /lustre/work/aucolema/enkfDART/src/sensvector <{} >sensvector.out".format(self._sensin)
        subprocess_cmd(args)
        return

    def runAll(self):
        '''
        Run entire suite of code that performs the following in order using
        default options:
            1. Renames WRF outfiles for sens and subset library
                naming conventions.
            2. Create input text file for fortran code. Defaults
                to 'esens.in'
            3. Initialize 'SENSmean.out' and 'Rmean.out' WRF
                outfiles, and then run meancalcSENS fortran
                executable.
            4. Run esensSPC fortran executable which
                currently stores output to 'wrfout.sens'.
                Also creates 'Rvals.nc' if last line in
                input file is set to True (is by default).
            5. Run sensvector fortran executable, which
                stores all individual member outer domain
                variables for sens time in a single file.
        Returns NULL.
        '''
        os.chdir(self._dir)
        print("Renaming WRF outfiles...")
        self.renameWRFOUT()
        print("Done with renaming!")
        self.createInfile()
        print("Created input file!")
        print("Starting runMeanCalc()...")
        self.runMeanCalc()
        print("Done with mean calculations!")
        print("Starting runSENS()...")
        self.runSENS()
        print("Done with sensitivity code!")
        print("Running storeSENSvals()...")
        self.storeSENSvals()
        print("Done storing SENS member values!")
        print("Succesfully completed all default SENS functions!")
        return
