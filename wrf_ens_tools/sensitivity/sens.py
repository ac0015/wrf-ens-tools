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
from wrf_ens_tools.post import subprocess_cmd, fromDatetime

package_dir = os.path.dirname(os.path.abspath(__file__))

#################
# Sens class
#################
class Sens:
    """
    The Sens class encapsulates the information pertaining to
    a forecast ensemble, response function, response function
    box, sensitivity time, and response time for running a
    sensitivity code package and producing a sensitivity
    object to be used by other applications.
    """
    def __init__(self, infile=False, gui=False, ensbasepath=None,
                 llat=31.578, ulat=35.578, llon=-105.855, ulon=-98.855,
                 rfuncstr="Max UH", sixhr=False, senstime=6, rtime=24,
                 rthresh=25., run=None, ensnum=42,
                 submit=datetime.utcnow()):
        """
        Constructor for an instance of the Sens class. Defaults to
        a response function box centered around Lubbock, TX with
        Max UH as the response and the most recent ensemble run to
        calculate sensitivity.

        Inputs
        ------
        infile -------------------- optional boolean specifying existence of
                                    sensitivity input file containing ensemble
                                    size, sensitivity time, response time,
                                    response function index, response box,
                                    boolean to store response values, and
                                    the response thresholds to use for UH
                                    and simulated reflectivity coverage
                                    respectively.
        gui ----------------------- optional boolean specifying existence of
                                    SENSEI GUI output file, which is then
                                    used to create the sensitivity input file.
        ensbasepath --------------- optional string specifying absolute path to
                                    the base directory of the ensemble.
        llat ---------------------- lower latitude of response box as a float.
                                    If infile or gui is set to True, this is
                                    pulled from the input file instead.
        ulat ---------------------- upper latitude of response box as a float.
                                    If infile or gui is set to True, this is
                                    pulled from the input file instead.
        llon ---------------------- lower longitude of response box as a float.
                                    If infile or gui is set to True, this is
                                    pulled from the input file instead.
        ulon ---------------------- upper longitude of response box as a float.
                                    If infile or gui is set to True, this is
                                    pulled from the input file instead.
        rfuncstr ------------------ string specifying response function to use.
                                    This is used to identify response index to
                                    write to the sensitivity input file.
        sixhr --------------------- boolean identifying whether to calculate
                                    a six-hour (True) or one-hour (False)
                                    response function.
        senstime ------------------ integer describing the sensitivity time
                                    in number of forecast hours from
                                    initialization.
        rtime --------------------- integer describing the response time in
                                    number of forecast hours from initialization.
        rthresh ------------------- response function threshold (only necessary
                                    for coverage response functions - magnitude
                                    response functions will ignore this
                                    parameter). For example, to calculate
                                    UH (>25) Coverage, you would need an rthresh
                                    of 25 m2/s2.
        run ----------------------- datetime object describing initializtion
                                    time of the run.
        ensnum -------------------- full ensemble size as an integer. This is
                                    written to the sensitivity input file.
        submit -------------------- records the time submitted.
        """
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
            # TO-DO: Replace this lazy feature
            fullpath = "/lustre/research/bancell/aucolema/HWT2016runs/"
            self._dir = fullpath + "{}{}{}{}/".format(yr, mo, day, hr)
        else:
            self._dir = ensbasepath

        # Initialize a few attributes
        self._rfuncstr = rfuncstr
        self._rthresh = rthresh

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
                # Only have coverage responses built for 6-hr rfuncs so far
                if sixhr:
                    self._uhthresh = float(esens_in[9])
                    self._dbzthresh = float(esens_in[10])
                    self._rfuncstr = self.getRString()
                    if 'UH' in self._rfuncstr:
                        self._rthresh = self._uhthresh
                    elif 'Refl' in self._rfuncstr:
                        self._rthresh = self._dbzthresh

        # If not, use user-specified inputs
        else:
            self._llat = llat
            self._ulat = ulat
            self._llon = llon
            self._ulon = ulon
            self._rtime = int(rtime)
            self._rfuncstr = rfuncstr
            self._ensnum = ensnum
            self._senstime = senstime
            self._rthresh = rthresh

        if sixhr:
            self._sixhr = True
            self._rfuncstr = "6-hr " + self._rfuncstr
        else:
            self._sixhr = False
            self._rfuncstr = "1-hr " + self._rfuncstr

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
        tostr = "Run initialized at: {} \n Response Box: {}\nResponse Function: {} at f{}".format(self._run, rbox,
        self._rfuncstr, self._rtime)
        if "Coverage" in self._rfuncstr:
            tostr += "\nResponse Threshold {}".format(self._rthresh)
        return tostr

    def setDir(self, dirpath):
        """
        Set base directory from which sensitivity code
        will run.
        """
        try:
            if os.path.exists(dirpath):
                self._dir = dirpath
            else:
                raise NameError(dirpath + " does not exist.")
        except:
            raise
        return

    def setInfile(self, absinfilepath):
        """
        Set the absolute path of the input file supplied to the
        sensitivity code.
        """
        self._sensin = absinfilepath
        return

    def setRefFileD1(self, absrefpathd1):
        """
        Set the WRF reference file path for the outer domain.
        """
        self._wrfrefd1 = absrefpathd1
        return

    def setRefFileD2(self, absrefpathd2):
        """
        Set the WRF reference file path for the inner domain.
        Also resets grid-spacing parameter.
        """
        self._wrfrefd2 = absrefpathd2
        return

    def setRFunc(self, index):
        """
        Set response function with an index.
        """
        rfuncstrs = {1 : "Avg Refl", 2 : "Max Refl",
                     3 : "Avg UH", 4 : "Max UH",
                     5 : "Accum PCP", 6 : "Avg Wind Spd",
                     7 : "UH Coverage", 8 : "Refl Coverage"}
        self._rfuncstr = rfuncstrs[index]

    def setSixHour(self, sixhr):
        """
        Set six hour boolean indicating whether to calculate
        the sensitivity for a one hour or six hour response
        function. If you've already run any sensitivity code,
        you will need to rerun it after changing this boolean.
        """
        self._sixhr = sixhr

    def getRIndex(self):
        """
        Returns the integer corresponding to the response function
        string provided to the constructor for use with the sensitivity
        code. These strings are specific to those used in the subset GUI.
        """
        rfuncinds = {"Avg Refl" : 1, "1-hr Max Refl" : 2,
                     "6-hr Max Refl" : 2,
                     "1-hr Avg UH" : 3, "1-hr Max UH" : 4,
                     "6-hr Avg UH" : 3, "6-hr Max UH" : 4,
                     "Accum PCP" : 5, "Avg Wind Spd" : 6,
                     "1-hr UH Coverage" : 7, "6-hr UH Coverage" : 7,
                     "1-hr Refl Coverage" : 8, "6-hr Refl Coverage" : 8}
        return rfuncinds[self._rfuncstr]


    def getRString(self):
        """
        Returns response function string
        """
        return self._rfuncstr

    def getRunInit(self):
        """
        Returns the initialization time of ensemble run as a datetime
        object.
        """
        return self._run

    def getSensTime(self):
        """
        Returns the sensitivity time as an integer of n forecast hrs
        from the run's initialization.
        """
        return self._senstime

    def getRTime(self):
        """
        Returns the response function time as an integer of n forecast
        hours from the run's initialization.
        """
        return self._rtime

    def getEnsnum(self):
        """
        Returns ensemble size as an integer.
        """
        return self._ensnum

    def getDir(self):
        """
        Returns base directory for ensemble run as a string.
        """
        return self._dir

    def getRbox(self):
        """
        Returns the response function box bounds as a tuple
        containing lower latitude, upper latitude, lower longitude,
        and upper longitude respectively.
        """
        return self._llon, self._ulon, self._llat, self._ulat

    def getInfile(self):
        """
        Returns the input file path provided to the sensitivity code.
        """
        return self._sensin

    def getRefFileD1(self):
        """
        Returns the WRF outer domain reference file path as a string.
        """
        return self._wrfrefd1

    def getRefFileD2(self):
        """
        Returns the WRF inner domain reference file path as a string.
        """
        return self._wrfrefd2

    def getWRFSensFile(self):
        """
        Returns the WRF sens file used by the Fortran sensitivity
        code. This is hardcoded in the Fortran, so only change
        if the Fortran naming conventions are changed.
        """
        return self._wrfsens

    def getSixHour(self):
        """
        Returns the boolean indicating whether response function
        is over a one hour or six hour time frame. If false,
        assume six hours.
        """
        return self._sixhr

    def getResponseThreshold(self):
        """
        Returns response function threshold. Only matters if sensitivity
        is being calculated for a coverage response function.
        """
        return self._rthresh

    def createInfile(self, rvals=True):
        """
        Creates (or overwrites if it already exists) an input file for
        the sensitivity code to use. If rvals is true in infile,
        sensitivity code will write response function values for each
        member to netCDF in 'Rvals.nc'.
        """
        fpath = self._sensin
        if os.path.isfile(fpath):
            os.remove(fpath)
        if rvals:
            rval = ".true."
        else:
            rval = ".false."
        args = [str(self._ensnum), str(self._senstime),
                str(self._rtime), str(self.getRIndex()),
                str(self._llon), str(self._ulon), str(self._llat),
                str(self._ulat), rval, str(self._rthresh),
                str(self._rthresh)]
        np.savetxt(fpath, args, fmt="%s", delimiter='\n')
        return

    def renameWRFOUT(self, restored=False):
        """
        Renames wrfout files for current ensemble run

        Inputs
        ------
        restored - optional boolean specifying if wrfout files
                    have been restored to their full capacity
                    from a previously reduced output file

        """
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
        """
        Runs the fortran program 'meancalcSENS' to calculate
        and store means of sensitivity variables for use with
        sensitivity code.
        """
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

        meancalcpath = os.path.join(package_dir, 'meancalcSENS')
        args = "{} >meancalcSENS.out".format(meancalcpath)
        subprocess_cmd(args)
        return

    def runSENS(self):
        """
        Runs the fortran sensitivity code and stores sensitivity
        netCDF file as the outfile dictated by the sensitivity object.
        """
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

        if self._sixhr:
            sens_exec = "sixhresens"
        else:
            sens_exec = "esensSPC"

        esenspath = os.path.join(package_dir, sens_exec)
        args = "{} <{} >esens.out".format(esenspath, self._sensin)
        subprocess_cmd(args)
        return

    def storeSENSvals(self):
        """
        Runs fortran code to store individual member sensitivity variable values
        to 'SENSvals.nc'.
        """
        os.chdir(self._dir)
        sensvecpath = os.path.join(package_dir, 'sensvector')
        args = "{} <{} >sensvector.out".format(sensvecpath,
                                                                self._sensin)
        subprocess_cmd(args)
        return

    def runAll(self):
        """
        Run entire suite of code that performs the following in order using
        default options:
            1. Renames WRF outfiles for sens and subset library
                naming conventions.
            2. Creates input text file for fortran code. Defaults
                to 'esens.in'
            3. Initializes 'SENSmean.out' and 'Rmean.out' WRF
                outfiles, and then runs then meancalcSENS fortran
                executable.
            4. Runs esensSPC (or sixhresens if six-hour response)
                fortran executable which currently stores
                output to 'wrfout.sens'.
                Also creates 'Rvals.nc' if last line in
                input file is set to True (which it is by default).
            5. Runs sensvector fortran executable, which
                stores all individual member outer domain
                variables for sens time in a single file.
        Returns NULL.
        """
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
