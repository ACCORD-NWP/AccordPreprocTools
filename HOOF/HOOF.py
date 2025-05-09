#!/usr/bin/python3
# *****************************************************************************************************************************************
# HOOF.py - version 1.10:
#    tool for homogenizing, dealiasing and superobing OPERA radar files 
#    P. Smerkol - May 2025
# *****************************************************************************************************************************************
import sys
import os
import re
import datetime
import warnings
import traceback
import time
import numpy as np
import scipy.optimize as optimize
import h5py

# *****************************************************************************************************************************************
# class HOOFHdfRadarAtribute: 
#    object representing a radar attribute from a HDF5 file 
# *****************************************************************************************************************************************
class HOOFHdfRadarAttribute(object):

   # ------------ constructor
   def __init__(self, site, groups, name, value):

      self.site = site             # name of the radar site that the attribute belongs to
      self.groups = list(groups)   # a list of group names to which the attribute belongs to
      self.name = name             # name of the attribute
      self.value = value           # value of the attribute

# *****************************************************************************************************************************************
# class HOOFSettings:
#    Class that holds settings read from the namelist file
#    Has static attributes so can be called from everywhere without initialization 
# *****************************************************************************************************************************************
class HOOFSettings(object):

   # ------- initialization of members (out of the constructor so we get static attributes)
   inFolder = ""                       # relative path to the folder with input hdf5 files
   outFolder = ""                      # relative path to the folder where the output files will be written
   settingsFileName = ""               # filename of the config file (must be in the same folder as HOOF)
   fileExtensions = []                 # tuple of file extensions representing valid HDF5 radar files
   superobQualTask = ""                # shortened name of task for the quality group used in superobing
   
   warningTag = ""                     # text printed next to warnings, to make them searchable
   errorTag = ""                       # text printed next to errors, to make them searchable
   printConsoleWarnings = True         # flag for writing warnings to the console output
   printLogWarnings = True             # flag for writing warnings to the log
   printConsoleErrors = True           # flag for writing errors to the console output
   printConsoleTiming = True           # flag for writing timing to the console output
   
   dbzNames = []                       # list of radar moment names that contain DBZ measurements
   thNames = []                        # list of radar moment names that contain TH measurements
   vradNames = []                      # list of radar moment names that contain VRAD measurements
   qualNames = []                      # names of quality groups attached to DBZ that we would like to keep (ROPO, BLOCK, SAT or TOTAL)
   reqAtts = []                        # list of HOOFHdfRadarAttribute objects that contain required attributes
   specAtts = []                       # list of lists of HOOFHdfRadarAttribute objects that contain radar specific attributes

   dealiasing  = True                  # flag to perform dealiasing
   windModelZSectorSize = -1           # size of z bin sector in wind model calculation
   windModelMinGoodPoints = -1         # minimum number of good points in a sector to calculate the wind model
   dealiasingMaxWind = -1              # maximum wind speed in m/s allowed after dealiasing procedure
   
   superobing = True                   # flag to perform superobing
   rangeBinFactor = -1                 # range bin multiplication factor for superobing
   rayAngleFactor = -1                 # ray angle multiplication factor for superobing
   maxArcSize = -1                     # maximum allowed arc size in metres for any superobing cell
   minQuality = -1                     # minimum quality to accept
   dbzThClearsky = -1                  # Dbz measurement threshold for clear sky in dB
   dbzPercentage = -1                  # Percentage of good points needed for each superob point for DBZ measurements
   vradPercentage = -1                 # Percentage of good points needed for each superob point for VRAD measurements
   vradMaxStd = -1                     # maximum stadard deviation of points allowed for each superob point for VRAD measurements

   reqRootGroups = []                  # list of unique required root groups
   reqDatasetGroups = []               # list of unique required dataset groups
   reqDataGroups = []                  # list of unique required data groups
   reqQualGroups = []                  # list of unique required quality groups

   earthRadius = 6371200               # Earth radius in metres
   Ke = 4/3                            # factor in Equivalent Earth calculation   

   # --------- converts a string to float
   @classmethod
   def __convert(cls, str):
      if str == "None":
         return None

      try:
         return float(str)
      except ValueError:
         print(cls.errorTag, " : conversion error in parsing namelist in parameter ", str)
         sys.exit()

   # --------- parses a boolean variable
   @classmethod
   def __parseBool(cls, line):
      if line == "TRUE" or line == "T":
         return True
      else:
         return False

   # --------- parses a multi parameter line
   @classmethod
   def __parseMultiLine(cls, line):
      return list(filter(None, line.replace("{", "").replace("}", "").split()))

   # --------- parses a line with a radar attribute
   @classmethod
   def __parseRadarAtt(cls, line):
      attList = list(filter(None, re.split("/| ", line)))
      name = attList[-3]
      value = cls.__convert(attList[-1])
      groups = attList[:-3]
      return name, value, groups

   # --------- parses the namelist and sets the appropriate variables
   @classmethod
   def Parse(cls, configFileName, inFolder, outFolder):

      # save the IO folders and settings file read from command line arguments
      cls.inFolder = inFolder
      cls.outFolder = outFolder
      cls.settingsFileName = configFileName

      # open the config file, read all lines and save the uncommented ones
      allLines = open(configFileName, 'r').read().splitlines()
      lines = [line.strip() for line in allLines if line.strip()[0] != '#']

      # find keyword indexes
      keywordIndexes = []
      for i in range(len(lines)):
         if '[' in lines[i] and ']' in lines[i]:
            keywordIndexes.append(i)
      keywordIndexes.append(len(lines))

      # loop on keyword indexes
      for i in range(len(keywordIndexes)-1):
         
         # get the current keyword and indexes
         index = keywordIndexes[i]
         nextindex = keywordIndexes[i+1]
         keyword = lines[index].strip(" []")
         
         # parse settings according to keyword
         if keyword == "File extensions to read":
            cls.fileExtensions = tuple(cls.__parseMultiLine(lines[index+1]))
         if keyword == "Superobing output quality task name":
            cls.superobQualTask = lines[index+1]
         if keyword == "Log keywords":
            for j in range(index+1, nextindex):
               line = lines[j].split()
               if line[0] == "WarningTag":
                  cls.warningTag = line[2]
               if line[0] == "ErrorTag":
                  cls.errorTag = line[2]
         if keyword == "Print warnings to console":
            cls.printConsoleWarnings = cls.__parseBool(lines[index+1])
         if keyword == "Print warnings to log":
            cls.printLogWarnings = cls.__parseBool(lines[index+1])
         if keyword == "Print errors to console":
            cls.printConsoleErrors = cls.__parseBool(lines[index+1])
         if keyword == "Print timing to console":
            cls.printConsoleTiming = cls.__parseBool(lines[index+1])            
         if keyword == "Radar moment names to save":
            for j in range(index+1, nextindex):
               line = cls.__parseMultiLine(lines[j])
               if line[0] == "DBZ":
                  cls.dbzNames = line[2:]
               if line[0] == "TH":
                  cls.thNames = line[2:]
               if line[0] == "VRAD":
                  cls.vradNames = line[2:]
         if keyword == "Required DBZ moment quality groups":
            cls.qualNames = cls.__parseMultiLine(lines[index+1])
         if keyword == "Attributes and default values":
            for j in range(index+1, nextindex):
               name, value, groups = cls.__parseRadarAtt(lines[j])
               attribute = HOOFHdfRadarAttribute("required", groups, name, value)
               cls.reqAtts.append(attribute)
         if "Attributes and default values -" in keyword:
            site = keyword.split()[-1]
            specAttributes = []
            for j in range(index+1, nextindex):
               name, value, groups = cls.__parseRadarAtt(lines[j])
               attribute = HOOFHdfRadarAttribute(site, groups, name, value)
               specAttributes.append(attribute)
            cls.specAtts.append(specAttributes)
         if keyword == "Dealiasing":
            cls.dealiasing = cls.__parseBool(lines[index+1])
         if keyword == "Height sector size in m":
            cls.windModelZSectorSize = cls.__convert(lines[index+1])
         if keyword == "Minimum good points in height sector":
            cls.windModelMinGoodPoints = cls.__convert(lines[index+1])
         if keyword == "Maximum dealiased wind speed in m/s":
            cls.dealiasingMaxWind = cls.__convert(lines[index+1])
         if keyword == "Superobing":
            cls.superobing = cls.__parseBool(lines[index+1])
         if keyword == "Range bin factor":
            cls.rangeBinFactor = cls.__convert(lines[index+1])
         if keyword == "Ray angle factor":
            cls.rayAngleFactor = cls.__convert(lines[index+1])
         if keyword == "Max arc size in m":
            cls.maxArcSize = cls.__convert(lines[index+1])
         if keyword == "Quality group to use":
            cls.superobQualName = lines[index+1]
         if keyword == "DBZ min quality":
            cls.minQuality = cls.__convert(lines[index+1])
         if keyword == "DBZ clear sky threshold":
            cls.dbzThClearsky = cls.__convert(lines[index+1])
         if keyword == "DBZ min percentage of good points":
            cls.dbzPercentage = cls.__convert(lines[index+1])
         if keyword == "VRAD min percentage of good points":
            cls.vradPercentage = cls.__convert(lines[index+1])
         if keyword == "VRAD max standard deviation":
            cls.vradMaxStd = cls.__convert(lines[index+1])

      # check if superob quality name exists in selected homogenization group and set the original task
      if cls.superobing:
         if cls.superobQualName in cls.qualNames:
            if cls.superobQualName == "ROPO":
               cls.superobQualName = "ropo"
            elif cls.superobQualName == "BLOCK":
               cls.superobQualName = "beamblockage"
            elif cls.superobQualName == "SAT":
               cls.superobQualName = "satfilter"
            elif cls.superobQualName == "TOTAL":
               cls.superobQualName = "qi_total"
         else:
            print(cls.errorTag, ": superobing is true, but the superob quality group is not in the selected homogenization groups")
            sys.exit()

      # determine the unique required subgroups in homogenization for each level
      for att in cls.reqAtts:
         if len(att.groups) == 1:
            cls.reqRootGroups.append(att.groups[0])
         if len(att.groups) == 2:
            cls.reqDatasetGroups.append(att.groups[1])
         if len(att.groups) == 3 and att.groups[1] == "data":
            cls.reqDataGroups.append(att.groups[2])
         if len(att.groups) == 3 and att.groups[1] == "quality":
            cls.reqQualGroups.append(att.groups[2])
      cls.reqRootGroups = list(set(cls.reqRootGroups))
      cls.reqDatasetGroups = list(set(cls.reqDatasetGroups))
      cls.reqDataGroups = list(set(cls.reqDataGroups))
      cls.reqQualGroups = list(set(cls.reqQualGroups))
         
      return

# *****************************************************************************************************************************************
# class HOOFWorker: 
#    general class that handles errors and warnings, used as base class for all worker classes
# *****************************************************************************************************************************************
class HOOFWorker(object):

   # ----------- constructor
   def __init__(self):

      self.errors = []          # list of error texts to display
      self.warnings = []        # list of warning texts to display

   # ------------ writes the warnings and errors to the standard output and log file
   def writeLog(self, logFile):
   
      # print warnings if set
      for warning in self.warnings:
         if HOOFSettings.printLogWarnings:
            logFile.write(HOOFSettings.warningTag + ": " + '%s\n'%(warning))
         if HOOFSettings.printConsoleWarnings:
            print("    " + HOOFSettings.warningTag + ": " + '%s'%(warning))
     
      # print errors
      for error in self.errors:
         logFile.write(HOOFSettings.errorTag + ": " + '%s\n'%(error))
         if HOOFSettings.printConsoleErrors:
            print("    " + HOOFSettings.errorTag + ": " + '%s'%(error))  

      return

# *****************************************************************************************************************************************
# class HOOFMeasurement: 
#    object that holds all data regarding one measurement (DBZ or VRAD) from one data volume
# *****************************************************************************************************************************************
class HOOFMeasurement(object):

   # ------------ constructor
   def __init__(self):

      self.datasets = None   # names of dataset groups belonging to the measurement in hdf5 file 
      self.qualdatas = None  # names of data groups for the TOTAL quality belonging to measurement
      self.nel = None        # number of elevations (datasets) in the HDF file after homogenization
      self.elangles = None   # elevation angles in radians for all elevations
      self.naz = None        # number of azimuths (rays) for all elevations
      self.azimuths = None   # azimuths in radians for all rays and elevations
      self.nr = None         # number of range bins for all elevations
      self.ranges = None     # bin ranges in metres for all range bins and elevations
      self.rscales = None    # range bin scales for all elevations
      self.rstarts = None    # range bin starts for all elevations
      self.vnys = None       # Nyquist velocities for all elevations, azimuths and range bins
      self.meas = None       # measurements of DBZ or VRAD for all elevations, azimuths and range bins
      self.ths = None        # values of TH corresponding to DBZ for all elevations, azimuths and range bins
      self.quals = None      # TOTAL quality values for all elevations, azimuths and range bins 
      self.zs = None         # heights for all elevations, ras and range bins

# *****************************************************************************************************************************************
# class HOOFData: 
#    object that holds all data transferred between different worker objects
# *****************************************************************************************************************************************
class HOOFData(object):

   # ------------ constructor
   def __init__(self): 

      self.radarsite = None             # name of the radar site taken from the OPERA filename
      self.radarheight = None           # height of the radar site in metres

      self.dbz = HOOFMeasurement()      # all data from DBZ measurements
      self.vrad = HOOFMeasurement()     # all data from VRAD measurements      
      
      self.zsectorstarts = None         # start heights of determined height ranges for dealiasing
      self.zsectorends = None           # end heights of determined height ranges for dealiasing
      self.zgoodindexes = None          # array of (el, az, r) indexes that are good for dealiasing for all height ranges in dynamic z height allocation
      self.windmodels = None            # values of the dealiasing wind model for all VRAD elevations, azimuths and range bins
      self.ns = None                    # Nyquist multipliers in dealiasing for all VRAD elevations, azimuths and range bins
      self.dvrads = None                # dealiased VRAD measurements for all VRAD elevations, azimuths and range bins

      self.sdbz = HOOFMeasurement()     # all data from superobed DBZ measurements
      self.svrad = HOOFMeasurement()    # all data from superobed VRAD measurements

   # ------------- stores homogenized data from a h5py object to this object
   def storeHomogenizedData(self, hdf):

      # get names and length of all VRAD and DBZ dataset groups and
      # corresponding quality data group to use in superobing
      self.dbz.datasets = []
      self.dbz.qualdatas = []
      self.vrad.datasets = []
      for rootGroup in hdf:
         if "dataset" in rootGroup:
            qtyName = hdf[rootGroup]["data1"]["what"].attrs["quantity"].decode()
            if qtyName == "DBZH":
               self.dbz.datasets.append(rootGroup)
               if HOOFSettings.superobing:
                  for dataGroup in hdf[rootGroup]:
                     if "quality" in dataGroup:
                        if HOOFSettings.superobQualName in hdf[rootGroup][dataGroup]["how"].attrs["task"].decode():
                           self.dbz.qualdatas.append(dataGroup)
                           break               
               continue
            if qtyName == "VRAD":
               self.vrad.datasets.append(rootGroup)
               continue
      self.dbz.nel = len(self.dbz.datasets)
      self.vrad.nel = len(self.vrad.datasets)

      # get radar site height
      self.radarheight = hdf["where"].attrs["height"]

      # fill DBZ quantity related data
      if self.dbz.nel > 0:

         # get dimensions
         self.dbz.naz = np.zeros(self.dbz.nel)
         self.dbz.nr = np.zeros(self.dbz.nel)
         for i in range(self.dbz.nel):
            self.dbz.naz[i] = hdf[self.dbz.datasets[i]]["where"].attrs["nrays"]
            self.dbz.nr[i] = hdf[self.dbz.datasets[i]]["where"].attrs["nbins"]
         self.dbz.naz = self.dbz.naz.astype(int)
         self.dbz.nr = self.dbz.nr.astype(int)
         nel = self.dbz.nel
         naz = np.max(self.dbz.naz)
         nr = np.max(self.dbz.nr)

         # read data from homogenized hdf object
         self.dbz.elangles = np.zeros(nel)
         self.dbz.azimuths = np.full((nel, naz), np.nan)
         self.dbz.ranges = np.full((nel, nr), np.nan)
         self.dbz.rstarts = np.zeros(nel)
         self.dbz.rscales = np.zeros(nel)
         self.dbz.meas = np.full((nel, naz, nr), np.nan)
         self.dbz.ths = np.full((nel, naz, nr), np.nan)
         self.dbz.quals = np.full((nel, naz, nr), np.nan)
         for i in range(self.dbz.nel):
            a = self.dbz.naz[i]
            r = self.dbz.nr[i]
            dataset = self.dbz.datasets[i]
            self.dbz.elangles[i] = hdf[dataset]["where"].attrs["elangle"] * np.pi / 180.
            self.dbz.azimuths[i,0:a] = np.linspace(0,360,a,endpoint=False) * np.pi / 180.
            rstart = hdf[dataset]["where"].attrs["rstart"]
            rscale = hdf[dataset]["where"].attrs["rscale"]
            self.dbz.rstarts[i] = rstart 
            self.dbz.rscales[i] = rscale
            self.dbz.ranges[i,0:r] = np.linspace(rstart, rstart + rscale*r, r, endpoint=False)
            dbzgain = hdf[dataset]["data1"]["what"].attrs["gain"]
            dbzoffset = hdf[dataset]["data1"]["what"].attrs["offset"]
            dbznodata = hdf[dataset]["data1"]["what"].attrs["nodata"]
            dbzraw = hdf[dataset]["data1"]["data"]
            dbzraw = np.where(dbzraw[:] == dbznodata, np.nan, dbzraw)
            self.dbz.meas[i,0:a,0:r] = dbzgain * dbzraw + dbzoffset
            thgain = hdf[dataset]["data2"]["what"].attrs["gain"]
            thoffset = hdf[dataset]["data2"]["what"].attrs["offset"]
            thnodata = hdf[dataset]["data2"]["what"].attrs["nodata"]
            thraw = hdf[dataset]["data2"]["data"]
            thraw = np.where(thraw[:] == thnodata, np.nan, thraw)
            self.dbz.ths[i,0:a,0:r] = thgain * thraw + thoffset
            if HOOFSettings.superobing:
               quality = self.dbz.qualdatas[i]
               qgain = hdf[dataset][quality]["what"].attrs["gain"]
               qoffset = hdf[dataset][quality]["what"].attrs["offset"]
               qraw = hdf[dataset][quality]["data"]
               qraw = np.where(dbzraw[:] == dbznodata, np.nan, qraw)
               self.dbz.quals[i,0:a,0:r] = qgain * qraw + qoffset

      # fill VRAD quantity related data
      if self.vrad.nel > 0:

         # get dimensions
         self.vrad.naz = np.zeros(self.vrad.nel)
         self.vrad.nr = np.zeros(self.vrad.nel)
         for i in range(self.vrad.nel):
            self.vrad.naz[i] = hdf[self.vrad.datasets[i]]["where"].attrs["nrays"]
            self.vrad.nr[i] = hdf[self.vrad.datasets[i]]["where"].attrs["nbins"]
         self.vrad.naz = self.vrad.naz.astype(int)
         self.vrad.nr = self.vrad.nr.astype(int)
         nel = self.vrad.nel
         naz = np.max(self.vrad.naz)
         nr = np.max(self.vrad.nr)

         # read data from homogenized hdf object
         self.vrad.elangles = np.zeros(nel)
         self.vrad.azimuths = np.full((nel, naz), np.nan)
         self.vrad.ranges = np.full((nel, nr), np.nan)
         self.vrad.rstarts = np.zeros(nel)
         self.vrad.rscales = np.zeros(nel)
         self.vrad.vnys = np.zeros(nel)
         self.vrad.meas = np.full((nel, naz, nr), np.nan)
         for i in range(self.vrad.nel):
            a = self.vrad.naz[i]
            r = self.vrad.nr[i]
            dataset = self.vrad.datasets[i]

            self.vrad.elangles[i] = hdf[dataset]["where"].attrs["elangle"] * np.pi / 180.
            self.vrad.azimuths[i,0:a] = np.linspace(0,360,a,endpoint=False) * np.pi / 180.
            rstart = hdf[dataset]["where"].attrs["rstart"]
            rscale = hdf[dataset]["where"].attrs["rscale"] 
            self.vrad.rstarts[i] = rstart
            self.vrad.rscales[i] = rscale
            self.vrad.ranges[i,0:r] = np.linspace(rstart, rstart + rscale*r, r, endpoint=False)
            self.vrad.vnys[i] = hdf[dataset]["how"].attrs["NI"]
            gain = hdf[dataset]["data1"]["what"].attrs["gain"]
            offset = hdf[dataset]["data1"]["what"].attrs["offset"]
            undetect = hdf[dataset]["data1"]["what"].attrs["undetect"]
            nodata = hdf[dataset]["data1"]["what"].attrs["nodata"]
            raw = hdf[dataset]["data1"]["data"]
            raw = np.where((raw[:] == undetect) | (raw[:] == nodata), np.nan, raw)
            self.vrad.meas[i,0:a,0:r] = gain * raw + offset

         # handle vrad measurements if we have VRAD data relative to Nyquist velocities
         if not np.isnan(self.vrad.meas).all():
            if np.nanmax(self.vrad.meas) <= 1.0:
               self.vrad.meas = self.vrad.meas * self.vrad.vnys.reshape(nel,1,1)

         # calculate heights for all VRAD measurements
         Ke = HOOFSettings.Ke
         R = HOOFSettings.earthRadius
         id = np.ones((nel,naz,nr))
         el = self.vrad.elangles.reshape(nel,1,1) * id
         r = self.vrad.ranges.reshape(nel,1,nr) * id
         self.vrad.zs = np.sqrt(np.square(r) + R*R*Ke*Ke + 2*r*R*Ke*np.sin(el)) - Ke*R + self.radarheight

# *****************************************************************************************************************************************
# class HOOFHomogenizationQuantity: 
#    object representing a measured quantity used in homogenization (written in a dataset/data or a dataset/quality group) that will be written to the output
#    stores the map between the input groups of the quantity and output groups of the quantity
# *****************************************************************************************************************************************
class HOOFHomogenizationQuantity(object):

   # ------------ constructor
   def __init__(self, name, date, angle, task, origGroups):

      self.name = name                     # name of the quantity (DBZH, TH, VRAD, QUALITYn)
      self.date = date                     # datetime object of the start date of the dataset to which the quantity belonged
      self.angle = angle                   # elevation angle of the dataset to which the quantity belonged
      self.task = task                     # shortened task name deduced from the quality how group
      self.origGroups = list(origGroups)   # list of group names that the quantity belonged to
      self.newGroups = []                  # list of group names that the quantity will belong to

# *****************************************************************************************************************************************
# class HOOFHomogenizer: 
#    object that homogenizes one HDF5 file with OPERA radar data
# *****************************************************************************************************************************************
class HOOFHomogenizer(HOOFWorker):

   # ------------ constructor
   def __init__(self, hdf, data):

      HOOFWorker.__init__(self)          # initialize the master class

      self.hdf = hdf                               # the h5py file object to homogenize
      self.data = data                             # object that holds all relevant data from worker objects
      self.qtys = []                               # list of all homogenization quantities found in the original file

   # -------------- transforms a group list to a form used in the specific and required attributes
   # -------------- (without numbers attached to dataset, data and quality groups)
   def __transformToReq(self, groups):

      reqGroups = []
      for group in groups:
         reqGroup = group
         if "dataset" in group:
            reqGroup = "dataset"
         if "data" in group and "dataset" not in group:
            reqGroup = "data"
         if "quality" in group:
            reqGroup = "quality"
         reqGroups.append(reqGroup)

      return reqGroups 

   # ------------- checks if a group exists in the HDF file
   def __checkGroup(self, groups):

      groupPath = "/".join(groups)
      if groupPath not in self.hdf:
         return False

      return True

   # ------------- checks if an attribute exists in the file, specific or required attributes
   def __checkAtt(self, groups, name):

      groupPath = "/".join(groups)
      reqGroups = self.__transformToReq(groups)

      # check if attribute exists in file
      if groupPath in self.hdf:
         if name in self.hdf[groupPath].attrs:
            return True

      # check if attribute exists in specific attributes
      specAtts = [attList for attList in HOOFSettings.specAtts if len(attList) > 0 if attList[0].site == self.data.radarsite]
      if len(specAtts) == 1:
         specAtt = [att for att in specAtts[0] if att.groups == reqGroups if att.name == name]
         if len(specAtt) == 1:
            if specAtt[0].value is not None:
               return True

      # check if attribute exists in required attributes
      comAtt = [att for att in HOOFSettings.reqAtts if att.groups == reqGroups if att.name == name]
      if len(comAtt) == 1:
         if comAtt[0].value is not None:
            return True

      # if the attribute does not exist anywhere
      return False

   # ------------- check all group attributes for existence
   def __checkGroupAtts(self, groups):

      reqGroups = self.__transformToReq(groups)
      names = [att.name for att in HOOFSettings.reqAtts if att.groups == reqGroups]
      failedAttNames = []
      for name in names:
         if not self.__checkAtt(groups, name):
               failedAttNames.append(name)

      return failedAttNames

   # ------------- get a specific attribute from a from a specific group represented by a group list
   # ------------- if it is not found, return default value taken from specific attributes
   # ------------- if it is also not found there, return the default value taken from required attributes
   # ------------- if it is not found at all, return None
   def __getAtt(self, groups, name):

      value = None

      # if the group and attribute are in the file, take the value from the file and return it
      groupPath = "/".join(groups)
      if groupPath in self.hdf:
         if name in self.hdf[groupPath].attrs:
            value = self.hdf[groupPath].attrs[name]
            return value

      # transform groups so that they can be searched in the specific and common parameters
      reqGroups = self.__transformToReq(groups)

      # else, take the attribute from specific attributes
      specAtts = [attList for attList in HOOFSettings.specAtts if attList[0].site == self.data.radarsite]
      if len(specAtts) == 1:
          specAtt = [att for att in specAtts[0] if att.groups == reqGroups if att.name == name]
          if len(specAtt) == 1:
             if specAtt[0].value is not None:
                value = specAtt[0].value
                return value

      # else, take the attribute from required attributes
      reqAtt = [att for att in HOOFSettings.reqAtts if att.groups == reqGroups if att.name == name]
      if len(reqAtt) == 1:
         if reqAtt[0].value is not None:
            value = reqAtt[0].value
            return value

      # else, return None
      return value

   # ------------- gets the start date and time as a datetime object from a dataset
   def __getStartDate(self, dataset):

      # get the startdate and starttime attributes
      dateGroups = [dataset, "what"]
      startDate = self.__getAtt(dateGroups, "startdate")
      startTime = self.__getAtt(dateGroups, "starttime")
      text = startDate + startTime

      # convert to a datetime object and return
      if startDate is not None and startTime is not None:
         year = int(text[0:4])
         month = int(text[4:6])
         day = int(text[6:8])
         hour = int(text[8:10])
         minute = int(text[10:12])
         if minute > 59:              # because we actually can get startdatetimes where this is true?!
            minute = 59
         second = int(text[12:14])
         if second > 59:              # because we actually can get startdatetimes where this is true?!
            second = 59

         # if we have a strange date
         if day < 1 or month < 1 or year < 1:
            return None

         dt = datetime.datetime(year, month, day, hour, minute, second)
         return dt

      return None   
      
   # ------------- gets the elevation angle rounded to 1 decimal place from a dataset
   def __getElAngle(self, dataset):

      # get the elangle attribute
      angleGroups = [dataset, "where"]
      angle = self.__getAtt(angleGroups, "elangle")

      # return elangle rounded to one decimal
      if angle is not None:
         angle = round(angle, 1)
         return angle

      return None

   # ------------- gets the shortened task name a quality group represented by a group list
   def __getTask(self, groups):

      taskGroups = list(groups)
      taskGroups.append("how")
      origTask = self.__getAtt(taskGroups, "task").decode()
      task = None
      if "ropo" in origTask:
         task = "ROPO"
      elif "beamblockage" in origTask:
         task = "BLOCK"
      elif "satfilter" in origTask:
         task = "SAT"
      elif "qi_total" in origTask:
         task = "TOTAL"

      return task

   # ------------- checks if the quantity in a group list is DBZ, TH or VRAD
   def __isQty(self, groups, qtyName):

      qtyNames = []
      if qtyName == "DBZH":
         qtyNames = HOOFSettings.dbzNames
      elif qtyName == "TH":
         qtyNames = HOOFSettings.thNames
      elif qtyName == "VRAD":
         qtyNames = HOOFSettings.vradNames

      qtyGroups = list(groups)
      qtyGroups.append("what")
      qty = self.__getAtt(qtyGroups, "quantity").decode()
      if qty in qtyNames:
         return True

      return False

   # ------------ validate a group and all its attributes for existence and write errors
   def __validateGroupAndAtts(self, groups):

      # validate group 
      if not self.__checkGroup(groups):
         self.errors.append("homogenization - required group /" + "/".join(groups) + " not found in file")
         return False

      # validate groups attributes
      failedAttNames = self.__checkGroupAtts(groups)
      if len(failedAttNames) > 0:
         self.errors.append("homogenization - attributes " + ",".join(failedAttNames) + " in group /" + "/".join(groups) + " not found in file, specific or required attributes")
         return False

      return True

   # ------------ Find all the DBZ, TH, VRAD and QUALITY quantities in the file satisfying the criteria.
   def sort(self):

      # loop over all "dataset" groups in the file
      dbzs = []
      ths = []
      vrads = []
      quals = []
      datasets = [group for group in self.hdf if "dataset" in group]
      for dataset in datasets:

         # skip the dataset if no date or elevation angle is found
         angle = self.__getElAngle(dataset)
         date = self.__getStartDate(dataset)
         if angle is None or date is None:
            self.warnings.append("homogenization - no date or elevation angle in dataset " + dataset + " - skipping")
            continue

         # get all "data" and "quality" groups in the current "dataset" group
         currDatasetGroup = self.hdf[dataset]
         dataGroups = [group for group in currDatasetGroup if "data" in group]
         qualGroups = [group for group in currDatasetGroup if "quality" in group]

         # save data groups into dbz, th and vrad quantities
         for dataGroup in dataGroups:
            origGroups = [dataset, dataGroup]
            if self.__isQty(origGroups, "DBZH"):
               nExistingDbz = len([dbz for dbz in dbzs if dbz.angle == angle and dbz.date == date])
               if nExistingDbz == 0:
                  dbzs.append(HOOFHomogenizationQuantity("DBZH", date, angle, None, origGroups))
            if self.__isQty(origGroups, "TH"):
               nExistingTh = len([th for th in ths if th.angle == angle and th.date == date])
               if nExistingTh == 0:
                  ths.append(HOOFHomogenizationQuantity("TH", date, angle, None, origGroups))
            if self.__isQty(origGroups, "VRAD"):
               nExistingVrad = len([vrad for vrad in vrads if vrad.angle == angle and vrad.date == date])
               if nExistingVrad == 0:
                  vrads.append(HOOFHomogenizationQuantity("VRAD", date, angle, None, origGroups))

         # save the required quality groups in quality quantities
         qualNumber = 0
         for qualGroup in qualGroups:
            origGroups = [dataset, qualGroup]
            task = self.__getTask(origGroups)
            nExistingQual = len([qual for qual in quals if qual.task == task and qual.angle == angle and qual.date == date])
            if nExistingQual == 0 and task in HOOFSettings.qualNames:
               qualNumber = qualNumber + 1
               quals.append(HOOFHomogenizationQuantity("QUALITY"+str(qualNumber), date, angle, task, origGroups))

      # sort dbz and vrad quantities by date
      dbzs.sort(key=lambda x: x.date)
      vrads.sort(key=lambda x: x.date)

      # determine the new groups for dbz and vrad quantities
      datasetCounter = 0
      for dbz in dbzs:
         datasetCounter = datasetCounter + 1
         dbz.newGroups = ["dataset"+str(datasetCounter), "data1"]
      for vrad in vrads:
         datasetCounter = datasetCounter + 1
         vrad.newGroups = ["dataset"+str(datasetCounter), "data1"]

      # put TH groups inside the corresponding DBZ group datasets (same elevation angle and start date)
      # if a TH group with no corresponding DBZ group is found, skip it
      ths2 = []
      for th in ths:
         dbzList = [dbz for dbz in dbzs if dbz.angle == th.angle and dbz.date == th.date]
         if len(dbzList) == 1:
            th.newGroups = [dbzList[0].newGroups[0], "data2"]
            ths2.append(th)
         else:
            self.warnings.append("homogenization - TH quantity in /" + "/".join(th.origGroups) + " has no matching DBZ group, omitting")
      
      # check that each DBZ quantity has a corresponding TH quantity, and skip it if not
      dbzs2 = []
      for dbz in dbzs:
         nCorrespTh = len([th for th in ths2 if th.newGroups[0] == dbz.newGroups[0]])
         if nCorrespTh == 0:
            self.warnings.append("homogenization - DBZ quantity " + "/".join(dbz.origGroups) + " has no corresponding TH group, omitting")
         else:
            dbzs2.append(dbz)

      # put QUALITY groups inside the corresponding DBZ group datasets (same elevation angle and start date)
      # if a QUALITY group with no corresponding DBZ group is found, skip it
      quals2 = []
      for qual in quals:
         dbzList = [dbz for dbz in dbzs if dbz.angle == qual.angle and dbz.date == qual.date]
         if len(dbzList) == 1:
            qual.newGroups = [dbzList[0].newGroups[0], "quality" + qual.name[7:]]
            quals2.append(qual)
         else:
            self.warnings.append("homogenization - QUALITY quantity in /" + "/".join(qual.origGroups) + " has no matching DBZ group, omitting")

      # check that each DBZ quantity has the required quality groups - otherwise omit the dataset that the DBZ group belongs to
      dbzs3 = []
      ths3 = []
      quals3 = []
      for dbz in dbzs2:
         qualNames = [qual.task for qual in quals2 if qual.angle == dbz.angle if qual.date == dbz.date]
         hasRequiredQualGroups = True
         for name in HOOFSettings.qualNames:
            if name not in qualNames:
               hasRequiredQualGroups = False
         if hasRequiredQualGroups:
            correspQual = [qual for qual in quals2 if qual.newGroups[0] == dbz.newGroups[0]]
            correspTh = [th for th in ths2 if th.newGroups[0] == dbz.newGroups[0]]
            dbzs3.append(dbz)
            ths3.extend(correspTh)
            quals3.extend(correspQual)
         else:
            self.warnings.append("homogenization - DBZ quantity in /" + "/".join(dbz.origGroups) + " does not have the required quality groups, omitting whole dataset")

      # recount the datasets and correct ther newGroups, so they always start with dataset 1
      finalDbzs = []
      finalThs = []
      finalQuals = []
      finalVrads = []
      datasetCounter = 0
      for dbz in dbzs3:
         datasetCounter = datasetCounter + 1
         newDatasetName = "dataset"+ str(datasetCounter)
         currDbz = dbz
         currThs = [th for th in ths3 if th.newGroups[0] == dbz.newGroups[0]]
         currQuals = [qual for qual in quals3 if qual.newGroups[0] == dbz.newGroups[0]]
         currDbz.newGroups = [newDatasetName, "data1"]
         for currTh in currThs:
            currTh.newGroups = [newDatasetName, "data2"]
         for currQual in currQuals:
            currQual.newGroups = [newDatasetName, "quality" + currQual.name[7:]]
         finalDbzs.append(currDbz)
         finalThs.extend(currThs)
         finalQuals.extend(currQuals)
      for vrad in vrads:
         datasetCounter = datasetCounter + 1
         newDatasetName = "dataset"+ str(datasetCounter)
         currVrad = vrad
         currVrad.newGroups = [newDatasetName, "data1"]
         finalVrads.append(currVrad)

      # save all quantity lists into the quantities list
      self.qtys.extend(finalDbzs)
      self.qtys.extend(finalThs)
      self.qtys.extend(finalQuals)
      self.qtys.extend(vrads)
      return

   # ------------ validates the file by checking that the required groups and attributes exist.
   def validate(self):

      # check the Conventions attribute
      conventions = None
      if "Conventions" in self.hdf.attrs:
         conventions = self.hdf.attrs["Conventions"].decode()
      if conventions != "ODIM_H5/V2_2":
         self.warnings.append("homogenization - Conventions attribute not set to ODIM_H5/V2_2, resetting.")

      # check the required root subgroups and attributes
      for group in HOOFSettings.reqRootGroups:
         currGroups = [group]
         self.__validateGroupAndAtts(currGroups)
      
      # check if there are any unique datasets to write (it is possible that all were omitted)
      datasets = [qty.origGroups[0] for qty in self.qtys]
      if len(datasets) == 0:
         self.errors.append("homogenization - no quantities to write to output file.")
         return

      # loop on all unique original datasets from quantities 
      datasets = set(datasets)
      for dataset in datasets:

         # check the required dataset subgroups and attributes
         for group in HOOFSettings.reqDatasetGroups:
            currGroups = [dataset, group]
            self.__validateGroupAndAtts(currGroups)

         # loop on all unique data groups from current unique dataset
         dataGroups = [qty.origGroups[1] for qty in self.qtys if dataset == qty.origGroups[0] if "data" in qty.origGroups[1]]
         dataGroups = set(dataGroups)
         for dataGroup in dataGroups:

            # check the required data subgroups and attributes
            for group in HOOFSettings.reqDataGroups:
               currGroups = [dataset, dataGroup, group]
               self.__validateGroupAndAtts(currGroups)

         # loop on all unique quality groups from current unique dataset
         qualGroups = [qty.origGroups[1] for qty in self.qtys if dataset == qty.origGroups[0] if "quality" in qty.origGroups[1]]
         qualGroups = set(qualGroups)
         for qualGroup in qualGroups:

            # check the required quality subgroups and attributes
            for group in HOOFSettings.reqQualGroups:
               currGroups = [dataset, qualGroup, group]
               self.__validateGroupAndAtts(currGroups)

      return

   # ------------ writes the homogenized data to a h5py object
   def writeTo(self, outHdf):

      # write the conventions attribute
      outHdf.attrs["Conventions"] = np.string_("ODIM_H5/V2_2")

      # create the required root subgroups and attributes
      for group in HOOFSettings.reqRootGroups:
         outHdf.create_group(group)
         atts = [att for att in HOOFSettings.reqAtts if att.groups == [group]]
         for att in atts:
            outHdf[group].attrs[att.name] = self.__getAtt([group], att.name)

      # loop on all quantities to write
      for qty in self.qtys:

         # create the dataset for the current quantity if it doesnt exist
         dataset = qty.newGroups[0]
         if dataset not in outHdf:
            outHdf.create_group(dataset)

         # create the required dataset subgroups if they dont exist
         for group in HOOFSettings.reqDatasetGroups:
            if group not in outHdf[dataset]:
               outHdf[dataset].create_group(group)

            # create the attributes of the dataset subgroups
            atts = [att for att in HOOFSettings.reqAtts if att.groups == ["dataset", group]]
            for att in atts:
               if att.name not in outHdf[dataset][group].attrs:
                  origGroups = [qty.origGroups[0], group]
                  outHdf[dataset][group].attrs[att.name] = self.__getAtt(origGroups, att.name)

         # create the data group in the dataset if it doesnt exist
         if "data" in qty.newGroups[1]:
            data = qty.newGroups[1]
            if data not in outHdf[dataset]:
               outHdf[dataset].create_group(data)

            # create the required data subgroups if they dont exist
            for group in HOOFSettings.reqDataGroups:
               if group not in outHdf[dataset][data]:
                  outHdf[dataset][data].create_group(group)

               # create the attributes of the data subgroups 
               atts = [att for att in HOOFSettings.reqAtts if att.groups == ["dataset", "data", group]]
               for att in atts:
                  if att.name not in outHdf[dataset][data][group].attrs:
                     origGroups = [qty.origGroups[0], qty.origGroups[1], group]
                     # for quantity attributes, we must set them to qty.name, because some TH quantities have the same original groups as DBZ quantities
                     if att.name == "quantity":
                        outHdf[dataset][data][group].attrs[att.name] = np.string_(qty.name)
                     # for other attributes, just copy from file
                     else:
                        outHdf[dataset][data][group].attrs[att.name] = self.__getAtt(origGroups, att.name)

            # create and copy the data array of the data group to the new file
            outHdf[dataset][data].create_dataset("data", data=self.hdf[qty.origGroups[0]][qty.origGroups[1]]["data"], compression="gzip", compression_opts=5)

         # create the quality group in the dataset if it doesnt exist
         if "quality" in qty.newGroups[1]:
            qual = qty.newGroups[1]
            if qual not in outHdf[dataset]:
               outHdf[dataset].create_group(qual)

            # create the required quality subgroups if they dont exist
            for group in HOOFSettings.reqQualGroups:
               if group not in outHdf[dataset][qual]:
                  outHdf[dataset][qual].create_group(group)

               # create the attributes of the quality subgroups
               atts = [att for att in HOOFSettings.reqAtts if att.groups == ["dataset", "quality", group]]
               for att in atts:
                  if att.name not in outHdf[dataset][qual][group].attrs:
                     origGroups = [qty.origGroups[0], qty.origGroups[1], group]
                     outHdf[dataset][qual][group].attrs[att.name] = self.__getAtt(origGroups, att.name)

            # create and copy the data array of the quality group to the new file
            outHdf[dataset][qual].create_dataset("data", data=self.hdf[qty.origGroups[0]][qty.origGroups[1]]["data"], compression="gzip", compression_opts=5)      

# *****************************************************************************************************************************************
# class HOOFDealiaser: 
#    object that calculates the dealiased VRAD fields
# *****************************************************************************************************************************************
class HOOFDealiaser(HOOFWorker):

   # ------------ constructor
   def __init__(self, hdf, data):

      HOOFWorker.__init__(self)      # initialize the master class

      self.hdf = hdf                           # h5py object to dealias
      self.data = data                         # object that holds all relevant data from worker objets
      self.cosEls3D = None                     # auxiliary quantity for faster code
      self.cosAzs3D = None                     # auxiliary quantity for faster code
      self.sinAzs3D = None                     # auxiliary quantity for faster code
      self.vNys3D = None                       # auxiliary quantity for faster code
      self.As = None                           # dealiasing a coefficients
      self.Bs = None                           # dealiasing b coefficients
      self.Ds = None                           # dealiasing D derivatives

   # ----------- the wind model function used in minimization
   def __windModelFunction(self, ab, u, v):
      
      a,b = ab
      return -a*u + b*v

   # ----------- checks the data if it is ok for dealiasing
   def checkData(self):

      if len(self.data.vrad.datasets) == 0:
         self.errors.append("dealiasing - no VRAD datasets in file")
         return

      if np.isnan(self.data.vrad.meas).all():
         self.errors.append("dealiasing - all VRAD data is NaN")
         return

   # ----------- calculates auxiliary and theoretical quantities for the wind model (a, b and D)
   def calculateAuxAndTheoreticalQtys(self):

      nel = self.data.vrad.nel
      naz = np.max(self.data.vrad.naz)
      nr = np.max(self.data.vrad.nr)
      id = np.ones((nel, naz, nr))
      self.cosEls3D = np.cos(self.data.vrad.elangles).reshape(nel,1,1) * id
      self.cosAzs3D = np.cos(self.data.vrad.azimuths).reshape(nel,naz,1) * id
      self.sinAzs3D = np.sin(self.data.vrad.azimuths).reshape(nel,naz,1) * id
      self.vNys3D = self.data.vrad.vnys.reshape(nel,1,1) * id
      self.As = self.cosEls3D * self.cosAzs3D * np.sin(np.pi * self.data.vrad.meas / self.vNys3D)
      self.Bs = self.cosEls3D * self.sinAzs3D * np.sin(np.pi * self.data.vrad.meas / self.vNys3D)
      f3s = self.vNys3D * np.cos(np.pi * self.data.vrad.meas / self.vNys3D) / np.pi
      df3s = np.roll(f3s,1,1) - np.roll(f3s,-1,1)
      azs = self.data.vrad.azimuths.reshape(nel,naz,1) * id
      dazs = np.roll(azs,1,1) - np.roll(azs,-1,1)
      dazs[:,0,:] = dazs[:,0,:] - 2*np.pi
      dazs[:,-1,:] = dazs[:,-1,:] - 2*np.pi
      self.Ds = df3s / dazs

   # ---------- determine height sectors
   def determineHeightSectors(self):
 
      # prepare some variables
      minGood = int(HOOFSettings.windModelMinGoodPoints)
      dZ = HOOFSettings.windModelZSectorSize
      zmax = np.nanmax(self.data.vrad.zs)
      zstart = self.data.radarheight
      zend = self.data.radarheight + dZ
      sectorstarts = []
      sectorends = []
      goodIndexes = []
      
      # determine where the nans are in measurements and D matrix
      nanIndexes = ~np.isnan(self.data.vrad.meas) & ~np.isnan(self.Ds)

      # determine points by height sectors  
      while True:

         # get good radar indexes for current height slice
         currIndexes = np.where(nanIndexes & (self.data.vrad.zs >= zstart) & (self.data.vrad.zs < zend))

         # if we have enough good points, declare a sector, save indexes and move to next slice
         if len(currIndexes[0]) > minGood:
            sectorstarts.append(zstart)
            sectorends.append(zend)
            goodIndexes.append(currIndexes)
         
         zstart = zend
         zend = zstart + dZ

         # if we come to the end, break
         if zend >= zmax:
            break

      # store sectors and indexes for later
      self.data.zsectorstarts = sectorstarts
      self.data.zsectorends = sectorends
      self.data.zgoodindexes = goodIndexes

   # ---------- calculate wind model
   def calculateWindModel(self):

      # prepare variables
      dvMax = float(HOOFSettings.dealiasingMaxWind)      
      self.data.windmodels = np.full(self.data.vrad.meas.shape, np.nan)

      # loop on height sectors
      newsectorstarts = []
      newsectorends = []
      goodIndexes = []      
      for z in range(len(self.data.zsectorstarts)-1):
         zstart = self.data.zsectorstarts[z]
         zend = self.data.zsectorends[z]
         currIndexes = self.data.zgoodindexes[z]
         As = self.As[currIndexes]
         Bs = self.Bs[currIndexes]
         Ds = self.Ds[currIndexes]
         cosEls = self.cosEls3D[currIndexes]
         cosAzs = self.cosAzs3D[currIndexes]
         sinAzs = self.sinAzs3D[currIndexes]

         with warnings.catch_warnings():
            warnings.filterwarnings('error')         

            # calculate u and v from minimization for each height sector
            try:
               uv, pcov = optimize.curve_fit( self.__windModelFunction, xdata=(As,Bs), ydata=Ds, p0=(0,0))
               u = uv[0]
               v = uv[1]
               vm = cosEls * (u * sinAzs + v * cosAzs)

               # exclude windmodels above max velocity and save sectors that have not been excluded               
               if np.min(vm) > -dvMax and np.max(vm) < dvMax:
                  self.data.windmodels[currIndexes] = vm
                  newsectorstarts.append(zstart)
                  newsectorends.append(zend)
                  goodIndexes.append(currIndexes)  

            except Warning as e:
               self.warnings.append("Fit not converged at height sector " + str(z) + "(" + str(zstart) + " " + str(zend) + ")")

      # store accepted sectors and indexes for later    
      self.data.zsectorstarts = newsectorstarts
      self.data.zsectorends = newsectorends
      self.data.zgoodindexes = goodIndexes

   # ------------ calculates dealiased VRAD
   def dealias(self):

      # prepare arrays
      self.data.ns = np.full(self.data.vrad.meas.shape, np.nan)      
      self.data.dvrads = np.full(self.data.vrad.meas.shape, np.nan)

      # calculate maximum possible Nyquist number
      dvMax = float(HOOFSettings.dealiasingMaxWind)
      vnyMin = np.min(self.data.vrad.vnys)
      nyMax = np.floor(dvMax/vnyMin)
      mn = np.full(self.data.vrad.meas.shape, 1000000.0)

      # dealias
      for m in np.arange(-nyMax, nyMax+1):
         mn_curr = np.abs(self.data.vrad.meas + 2.0*self.vNys3D*m - self.data.windmodels)
         self.data.ns = np.where(mn_curr < mn, m, self.data.ns)
         mn = np.where(mn_curr < mn, mn_curr, mn)
      self.data.dvrads = self.data.vrad.meas + 2.0 * self.data.ns * self.vNys3D

   # ---------- writes the dealiased data in the VRAD data group
   def writeTo(self, hdf):

      newVradDatasets = []
      for i in range(len(self.data.vrad.datasets)):
         dataset = self.data.vrad.datasets[i]
         del hdf[dataset]["data1"]["data"]
         values = self.data.dvrads[i,0:self.data.vrad.naz[i],0:self.data.vrad.nr[i]]
         gain = 1.
         offset = 0.0
         if not np.isnan(values).all():
            minvalue = np.nanmin(values)
            maxvalue = np.nanmax(values)
            gain = (maxvalue - minvalue) / 253.0
            if gain == 0.0:
               gain = 1.0
            offset = (254.0*minvalue - maxvalue) / 253.0
         hdf[dataset]["data1"]["what"].attrs["gain"] = gain
         hdf[dataset]["data1"]["what"].attrs["offset"] = offset
         hdf[dataset]["data1"]["what"].attrs["nodata"] = 255.0
         hdf[dataset]["data1"]["what"].attrs["undetect"] = 0.0
         data = (values - offset + 0.5*gain) / gain
         data = np.nan_to_num(data, nan=255)
         hdf[dataset]["data1"].create_dataset("data", data=np.uint8(data), compression="gzip", compression_opts=6)
         newVradDatasets.append(dataset)
      self.data.vrad.datasets = newVradDatasets

      return hdf           

# *****************************************************************************************************************************************
# class HOOFSuperOber: 
#    object that determines the super observations for DBZ and VRAD fields
# *****************************************************************************************************************************************
class HOOFSuperOber(HOOFWorker):

   # ------------ constructor
   def __init__(self, hdf, data):

      HOOFWorker.__init__(self)          # initialize the master class

      self.hdf = hdf                               # h5py object to superob
      self.data = data                             # object that holds all data
      self.dbzsNan = False                         # if all DBZ data is nan
      self.vradsNan = False                        # if all VRAD data is nan

   # ------------ checks data if it is ok for superobing
   def checkData(self):

      if self.data.dbz.nel > 0:
         if np.isnan(self.data.dbz.meas).all():
            self.dbzsNan = True

      if self.data.vrad.nel > 0:
         if np.isnan(self.data.vrad.meas).all():
            self.vradsNan = True
      
      if self.data.dbz.nel == 0 and self.data.vrad.nel == 0:
         self.errors.append("superobing - no data to superob")
         return

      if self.dbzsNan and self.vradsNan:
         self.errors.append("superobing - all data is NaN")
         return   

      if self.dbzsNan and not self.vradsNan:
         self.warnings.append("superobing - all DBZ data is NaN")

      if not self.dbzsNan and self.vradsNan:
         self.warnings.append("superobing - all VRAD data is NaN")

   # ------------ get the superobed data
   def getSuperobData(self):

      # short aliases
      binFactor = HOOFSettings.rangeBinFactor
      rayFactor = HOOFSettings.rayAngleFactor

      # DBZ measurements
      self.data.sdbz.nel = self.data.dbz.nel
      if self.data.dbz.nel > 0:

         # get dimensions
         self.data.sdbz.nr = (self.data.dbz.nr / binFactor).astype(int)
         self.data.sdbz.naz = (self.data.dbz.naz / rayFactor).astype(int)
         nel = self.data.sdbz.nel
         naz = np.max(self.data.sdbz.naz)
         nr = np.max(self.data.sdbz.nr)
 
         # store the needed data
         self.data.sdbz.elangles = self.data.dbz.elangles
         self.data.sdbz.azimuths = np.full((nel, naz), np.nan)
         self.data.sdbz.rscales = np.zeros(nel)
         self.data.sdbz.rstarts = np.zeros(nel)
         self.data.sdbz.ranges = np.full((nel, nr), np.nan)
         for i in range(nel):
            a = self.data.sdbz.naz[i]
            r = self.data.sdbz.nr[i]
            rstart = self.data.dbz.rstarts[i]
            rscale = self.data.dbz.rscales[i]
            self.data.sdbz.azimuths[i,0:a] = np.linspace(0,360,a,endpoint=False) * np.pi / 180.
            self.data.sdbz.rstarts[i] = rstart
            self.data.sdbz.rscales[i] = binFactor * rscale
            self.data.sdbz.ranges[i,0:r] = np.linspace(rstart, rstart + binFactor*rscale*r, r, endpoint=False)

      # VRAD measurements
      self.data.svrad.nel = self.data.vrad.nel      
      if self.data.vrad.nel > 0:

         # get dimensions
         self.data.svrad.nr = (self.data.vrad.nr / binFactor).astype(int)
         self.data.svrad.naz = (self.data.vrad.naz / rayFactor).astype(int)
         nel = self.data.svrad.nel
         naz = np.max(self.data.svrad.naz)
         nr = np.max(self.data.svrad.nr)
 
         # store the needed data
         self.data.svrad.elangles = self.data.vrad.elangles
         self.data.svrad.azimuths = np.full((nel, naz), np.nan)
         self.data.svrad.rscales = np.zeros(nel)
         self.data.svrad.rstarts = np.zeros(nel)
         self.data.svrad.ranges = np.full((nel, nr), np.nan)
         for i in range(nel):
            a = self.data.svrad.naz[i]
            r = self.data.svrad.nr[i]
            rstart = self.data.vrad.rstarts[i]
            rscale = self.data.vrad.rscales[i]
            self.data.svrad.azimuths[i,0:a] = np.linspace(0,360,a,endpoint=False) * np.pi / 180.
            self.data.svrad.rstarts[i] = rstart
            self.data.svrad.rscales[i] = binFactor * rscale
            self.data.svrad.ranges[i,0:r] = np.linspace(rstart, rstart + binFactor*rscale*r, r, endpoint=False)

   # ------------ do the superobing
   def superob(self):

      # short aliases
      arclim = HOOFSettings.maxArcSize
      clearsky = HOOFSettings.dbzThClearsky
      qualth = HOOFSettings.minQuality
      dbzpercent = HOOFSettings.dbzPercentage
      vradpercent = HOOFSettings.vradPercentage
      vradmaxstd = HOOFSettings.vradMaxStd
      dbzMinimum = np.nanmin(self.data.dbz.meas)
      Fb = int(HOOFSettings.rangeBinFactor)
      Fr = int(HOOFSettings.rayAngleFactor)
      Zmax = int((Fr-1)/2)

      # superob the DBZ measurements
      if self.data.dbz.nel > 0:

         # prepare the superobed arrays
         self.data.sdbz.meas = np.full((self.data.sdbz.nel,np.max(self.data.sdbz.naz),np.max(self.data.sdbz.nr)), np.nan)
         self.data.sdbz.ths = np.full((self.data.sdbz.nel,np.max(self.data.sdbz.naz),np.max(self.data.sdbz.nr)), np.nan)
         self.data.sdbz.quals = np.full((self.data.sdbz.nel,np.max(self.data.sdbz.naz),np.max(self.data.sdbz.nr)), 0.)

         # roll the matrices so we get the right ray positions
         meas = np.roll(self.data.dbz.meas, Zmax, 1)
         ths = np.roll(self.data.dbz.ths, Zmax, 1)
         quals = np.roll(self.data.dbz.quals, Zmax, 1)
         
         # loop on elevations
         for el in range(self.data.sdbz.nel):    

            # short aliases
            naz = self.data.dbz.naz[el]
            nr = self.data.dbz.nr[el]
            rscale = self.data.dbz.rscales[el]
            L = 360.*360.*arclim/(2.*np.pi*naz*Fb*rscale)
            snaz = self.data.sdbz.naz[el]
            snr = self.data.sdbz.nr[el]

            # calculate bin and ray borders for superob points
            sbins = np.array(range(0, nr + Fb, Fb)).astype(int)
            sbins = sbins[sbins <= nr]
            origstartrays = np.array(range(0, naz, Fr))
            origendrays = np.array(range(Fr, naz+Fr, Fr))          
            limindexes = [0]
            facsubs = []
            for z in range(Zmax+1):
               newfac = 2*(Zmax-z) + 1
               limindex = np.floor(L / newfac - 1) + 1
               if limindex >= len(sbins):
                  limindex = len(sbins)
               if limindex != limindexes[-1]:
                  limindexes.append(limindex)
                  facsubs.append(z)  
            if limindexes[-1] < len(sbins):
               limindexes[-1] = len(sbins)
            limindexes = np.array(limindexes)
            facsubs = np.array(facsubs)
            startrays = np.full((snr, snaz), np.nan).astype(int)
            endrays = np.full((snr, snaz), np.nan).astype(int)
            currfacsub = -1
            for i in range(snr):
               for j in range(len(facsubs)):
                  if i >= limindexes[j] and i < limindexes[j+1]:
                     currfacsub = facsubs[j]
               startrays[i,:] = origstartrays + currfacsub
               endrays[i, :] = origendrays - currfacsub

            # loop on superobing bins
            for i in range(snr):
               startbin = sbins[i]
               endbin = sbins[i+1]                  
               
               for j in range(snaz):               
                  startray = startrays[i,j]
                  endray = endrays[i,j]

                  # get the matrices for one superob point
                  currmeas = meas[el, startray:endray, startbin:endbin]
                  currths = ths[el, startray:endray, startbin:endbin]
                  currquals = quals[el, startray:endray, startbin:endbin]
                  nrays = endray-startray
                  nbins = endbin-startbin

                  # count the rainy and dry points and calculate average of rainy points
                  nrain = 0
                  ndry = 0
                  avg = 0.0
                  nrainth = 0
                  avgth = 0.0
                  for k in range(nrays):
                     for l in range(nbins):
                        q = currquals[k,l]
                        m = currmeas[k,l]
                        t = currths[k,l]
                        if q > qualth:
                           if m > clearsky:
                              nrain = nrain + 1
                              avg = avg + m
                              if t < 100000.0:
                                 nrainth = nrainth + 1
                                 avgth = avgth + t
                           else:
                              ndry = ndry + 1

                  # superob
                  if nrain > dbzpercent * nrays * nbins:
                     self.data.sdbz.meas[el, j, i] = avg/nrain
                     if nrainth > 0:
                        self.data.sdbz.ths[el, j, i] = avgth/nrainth
                     self.data.sdbz.quals[el,j, i] = 1.                     
                  else:
                     if ndry > 0:
                        self.data.sdbz.meas[el, j, i] = dbzMinimum
                        self.data.sdbz.quals[el,j, i] = 1.                     

      # superob the VRAD measurements
      if self.data.svrad.nel > 0:

         # prepare the superobed arrays
         self.data.svrad.meas = np.full((self.data.svrad.nel,np.max(self.data.svrad.naz),np.max(self.data.svrad.nr)), np.nan)
         self.data.svrad.quals = np.full((self.data.svrad.nel,np.max(self.data.svrad.naz),np.max(self.data.svrad.nr)), 0.)

         # roll the matrices so we get the right ray positions
         meas = None
         if not HOOFSettings.dealiasing:
            meas = np.roll(self.data.vrad.meas, Zmax, 1)
         else:
            meas = np.roll(self.data.dvrads, Zmax, 1)

         # loop on elevations
         for el in range(self.data.vrad.nel):

            # short aliases
            naz = self.data.vrad.naz[el]
            nr = self.data.vrad.nr[el]
            rscale = self.data.vrad.rscales[el]
            L = 360.*360.*arclim/(2.*np.pi*naz*Fb*rscale)
            snaz = self.data.svrad.naz[el]
            snr = self.data.svrad.nr[el]

            # calculate bin and ray borders for superob points
            sbins = np.array(range(0, nr + Fb, Fb)).astype(int)
            sbins = sbins[sbins <= nr]
            origstartrays = np.array(range(0, naz, Fr)).astype(int)
            origendrays = np.array(range(Fr, naz+Fr, Fr)).astype(int)            
            limindexes = [0]
            facsubs = []
            for z in range(Zmax+1):
               newfac = 2*(Zmax-z) + 1
               limindex = np.floor(L / newfac - 1) + 1
               if limindex >= len(sbins):
                  limindex = len(sbins)
               if limindex != limindexes[-1]:
                  limindexes.append(limindex)
                  facsubs.append(z)  
            if limindexes[-1] < len(sbins):
               limindexes[-1] = len(sbins)
            limindexes = np.array(limindexes)
            facsubs = np.array(facsubs)
            startrays = np.full((snr, snaz), np.nan).astype(int)
            endrays = np.full((snr, snaz), np.nan).astype(int)
            currfacsub = -1
            for i in range(snr):
               for j in range(len(facsubs)):
                  if i >= limindexes[j] and i < limindexes[j+1]:
                     currfacsub = facsubs[j]
               startrays[i,:] = origstartrays + currfacsub
               endrays[i, :] = origendrays - currfacsub

            # loop on superobing bins
            for i in range(snr):
               startbin = sbins[i]
               endbin = sbins[i+1] 

               for j in range(snaz):
                  startray = startrays[i,j]
                  endray = endrays[i,j]

                  # get the matrix for one superob point
                  currmeas = meas[el, startray:endray, startbin:endbin]
                  nrays = endray-startray
                  nbins = endbin-startbin

                  # count the good points and calculate average and standard deviation
                  ngood = 0
                  s = 0
                  s2 = 0
                  avg = 0
                  std = vradmaxstd + 1.0                  
                  for k in range(nrays):
                     for l in range(nbins):
                        m = currmeas[k,l]
                        if m < 1000000.0:   # this check is for nan values (only nan is false in this case)
                           ngood = ngood + 1
                           s = s + m
                           s2 = s2 + m*m
                  if ngood > 0:
                     avg = s / ngood
                     std2 = (s2 - s*s/ngood) / ngood
                     std = np.sqrt(std2) if std2 > 0.0 else 0.0

                  # sueprob
                  if ngood > vradpercent * nrays * nbins and std < vradmaxstd:
                     self.data.svrad.meas[el, j, i] = avg
                     self.data.svrad.quals[el, j, i] = 1.

   # --- write the superobed data to file
   def writeTo(self, hdf):

      for i in range(len(self.data.dbz.datasets)):
         dataset = self.data.dbz.datasets[i]
         qualityGroup = self.data.dbz.qualdatas[i]
         nr = self.data.sdbz.nr[i]
         naz = self.data.sdbz.naz[i]
         rscale = self.data.sdbz.rscales[i]
         hdf[dataset]["where"].attrs["nbins"] = nr
         hdf[dataset]["where"].attrs["nrays"] = naz
         hdf[dataset]["where"].attrs["rscale"] = rscale
         # DBZ group
         del hdf[dataset]["data1"]["data"]
         values = self.data.sdbz.meas[i,0:naz,0:nr]
         nodata = hdf[dataset]["data1"]["what"].attrs["nodata"]
         gain = 1.0
         offset = 0.0
         if not np.isnan(values).all():
            minvalue = np.nanmin(values)
            maxvalue = np.nanmax(values)
            gain = (maxvalue - minvalue) / 254.0
            if gain == 0.0:
               gain = 1.0            
            offset = (254.0*minvalue - maxvalue) / 253.0         
         data = (values - offset + 0.5*gain) / gain
         data = np.nan_to_num(data, nan=nodata)
         hdf[dataset]["data1"]["what"].attrs["gain"] = gain
         hdf[dataset]["data1"]["what"].attrs["offset"] = offset         
         hdf[dataset]["data1"].create_dataset("data", data=np.uint8(data), compression="gzip", compression_opts=6)
         # TH group
         del hdf[dataset]["data2"]["data"]
         values = self.data.sdbz.ths[i,0:naz,0:nr]
         nodata = hdf[dataset]["data2"]["what"].attrs["nodata"]
         gain = 1.0
         offset = 0.0
         if not np.isnan(values).all():
            minvalue = np.nanmin(values)
            maxvalue = np.nanmax(values)
            gain = (maxvalue - minvalue) / 254.0
            if gain == 0.0:
               gain = 1.0            
            offset = (254.0*minvalue - maxvalue) / 253.0           
         data = (values - offset + 0.5*gain) / gain
         data = np.nan_to_num(data, nan=nodata)
         hdf[dataset]["data2"]["what"].attrs["gain"] = gain
         hdf[dataset]["data2"]["what"].attrs["offset"] = offset         
         hdf[dataset]["data2"].create_dataset("data", data=np.uint8(data), compression="gzip", compression_opts=6)            
         # quality group
         for group in hdf[dataset]:
            if "quality" in group:
               if qualityGroup in group:
                  del hdf[dataset][group]['data']
               else:
                  del hdf[dataset][group]
         values = self.data.sdbz.quals[i,0:naz,0:nr]
         gain = 1. / 255.
         offset = 0.0
         data = (values - offset + 0.5*gain) / gain
         data = np.nan_to_num(data, nan=0.0)
         hdf[dataset][qualityGroup]["what"].attrs["gain"] = gain
         hdf[dataset][qualityGroup]["what"].attrs["offset"] = offset            
         hdf[dataset][qualityGroup].create_dataset("data", data=np.uint8(data), compression="gzip", compression_opts=6)
         hdf[dataset][qualityGroup]["how"].attrs["task"] = np.string_(HOOFSettings.superobQualTask)
         if qualityGroup != "quality1":
            hdf[dataset].move(qualityGroup, "quality1")              

      for i in range(len(self.data.vrad.datasets)):
         dataset = self.data.vrad.datasets[i]
         nr = self.data.svrad.nr[i]
         naz = self.data.svrad.naz[i]
         rscale = self.data.svrad.rscales[i]
         hdf[dataset]["where"].attrs["nbins"] = nr
         hdf[dataset]["where"].attrs["nrays"] = naz
         hdf[dataset]["where"].attrs["rscale"] = rscale
         # VRAD group
         del hdf[dataset]["data1"]["data"]
         values = self.data.svrad.meas[i,0:naz,0:nr]
         gain = 1.0
         offset = 0.0
         if not np.isnan(values).all():
            minvalue = np.nanmin(values)
            maxvalue = np.nanmax(values)
            gain = (maxvalue - minvalue) / 254.0
            if gain == 0.0:
               gain = 1.0
            offset = (254.0*minvalue - maxvalue) / 253.0
         hdf[dataset]["data1"]["what"].attrs["gain"] = gain
         hdf[dataset]["data1"]["what"].attrs["offset"] = offset
         hdf[dataset]["data1"]["what"].attrs["nodata"] = 255.0
         hdf[dataset]["data1"]["what"].attrs["undetect"] = 0.0
         data = (values - offset + 0.5*gain) / gain
         data = np.nan_to_num(data, nan=255)
         hdf[dataset]["data1"].create_dataset("data", data=np.uint8(data), compression="gzip", compression_opts=6)
         # quality group
         values = self.data.svrad.quals[i,0:naz,0:nr]
         gain = 1. / 255.
         offset = 0.0
         hdf[dataset].create_group("quality1")
         hdf[dataset]["quality1"].create_group("how")
         hdf[dataset]["quality1"].create_group("what")
         hdf[dataset]["quality1"]["what"].attrs["gain"] = gain
         hdf[dataset]["quality1"]["what"].attrs["offset"] = offset
         hdf[dataset]["quality1"]["how"].attrs["task"] = np.string_(HOOFSettings.superobQualTask)
         data = (values - offset) / gain
         data = np.nan_to_num(data, nan=0)
         hdf[dataset]["quality1"].create_dataset("data", data=np.uint8(data), compression="gzip", compression_opts=6)           

# --- handles errors in analysis from worker classes
def handleAnalysisErrors(worker, inHdf, outHdf, logFile, outFilePath):

   if len(worker.errors) > 0:
      worker.writeLog(logFile)
      inHdf.close()
      outHdf.close()
      logFile.close()
      os.remove(outFilePath)
      return False

   return True 

# *****************************************************************************************************************************************
#  main function
# *****************************************************************************************************************************************
if __name__=="__main__":
   
   # get the command line arguments and parse settings or print the instructions
   if len(sys.argv) == 4:
      configFileName = sys.argv[1]
      inFolder = sys.argv[2]
      outFolder = sys.argv[3]
      HOOFSettings.Parse(configFileName, inFolder, outFolder)
   else:
      print("Wrong number of command line arguments, the syntax is:")
      print( "./HOOF2.py <namelist file> <input folder> <output folder>")
      print("Input filenames must be in Opera format:")
      print("    first 15 characters - radar identifier (T_PAZZ42_C_LJLM)")
      print("    next 14 characters - date and time in YYYYMMDDhhmmss (20181027120000)")
      print("    next 5 characters - Opera NOD identifier (sipas)")
      print("    separated by _ and appended with an extension from the settings")      
      sys.exit()

   startTime = time.time()

   # get the list of input filenames and loop on them
   fileCounter = 0
   goodFileCounter = 0
   inputFileNames = [name for name in os.listdir(inFolder) if name.endswith(HOOFSettings.fileExtensions)]
   for inputFileName in inputFileNames:

      # print some running info
      fileCounter = fileCounter + 1
      print("----------- processing file ", inputFileName, "[", str(fileCounter), "/", str(len(inputFileNames)), "] -----------") 

      beginTime = time.time()
      t = [None] * 16
      
      try: 
         # read input file, prepare output file and log file and get radar site from current filename
         t[0] = time.perf_counter()
         print("    Reading input file ...")
         inFilePath = HOOFSettings.inFolder + "/" + inputFileName
         outFilePath = HOOFSettings.outFolder + "/" + inputFileName
         logFilePath = os.path.splitext(outFilePath)[0] + '.log'
         inHdf = h5py.File(inFilePath, 'r', driver='core', backing_store = False)
         outHdf = h5py.File(outFilePath, "w")
         logFile = open(logFilePath, "w")
         data = HOOFData()
         data.radarsite = inputFileName[31:36]

         # homogenize data
         t[1] = time.perf_counter()
         homogenizer = HOOFHomogenizer(inHdf, data)
         print("    Sorting datasets for homogenization ...")
         homogenizer.sort()
         t[2] = time.perf_counter()
         print("    Validating homogenization data ... ")
         homogenizer.validate()
         if not handleAnalysisErrors(homogenizer, inHdf, outHdf, logFile, outFilePath):
            continue
         t[3] = time.perf_counter()
         print("    Writing homogenized data ...")
         homogenizer.writeTo(outHdf)
         homogenizer.writeLog(logFile)
         t[4] = time.perf_counter()

         # store homogenized data for further use
         if HOOFSettings.dealiasing or HOOFSettings.superobing:
            print("    Storing homogenized data for further use ...")
            data.storeHomogenizedData(outHdf)
            t[5] = time.perf_counter()

         # dealias radial velocities
         if HOOFSettings.dealiasing:
            dealiaser = HOOFDealiaser(outHdf, data)
            print("    Checking data for dealiasing ...")
            dealiaser.checkData()
            if not handleAnalysisErrors(dealiaser, inHdf, outHdf, logFile, outFilePath):
               continue
            t[6] = time.perf_counter()
            print("    Calculating theoretical quantities for dealiasing ...")
            dealiaser.calculateAuxAndTheoreticalQtys()
            t[7] = time.perf_counter()
            print("    Determining height sectors for dealiasing ...")
            dealiaser.determineHeightSectors()
            t[8] = time.perf_counter()
            print("    Calculating dealiasing wind model ....")
            dealiaser.calculateWindModel()
            if not handleAnalysisErrors(dealiaser, inHdf, outHdf, logFile, outFilePath):
               continue
            t[9] = time.perf_counter()
            print("    Dealiasing ...")
            dealiaser.dealias()
            if not handleAnalysisErrors(dealiaser, inHdf, outHdf, logFile, outFilePath):
               continue
            t[10] = time.perf_counter()
            print("    Writing dealiased data ...")
            dealiaser.writeTo(outHdf)
            dealiaser.writeLog(logFile)
            t[11] = time.perf_counter()

         # superobing
         if HOOFSettings.superobing:
            superober = HOOFSuperOber(outHdf, data)
            print("    Checking data for superobing ...")
            superober.checkData()
            if not handleAnalysisErrors(superober, inHdf, outHdf, logFile, outFilePath):
               continue
            t[12] = time.perf_counter()
            print("    Preparing data for superobing ...")
            superober.getSuperobData()
            t[13] = time.perf_counter()
            print("    Superobing ...")
            superober.superob()
            t[14] = time.perf_counter()
            print("    Writing superobed data ...")
            superober.writeTo(outHdf)
            superober.writeLog(logFile)
            t[15] = time.perf_counter()

      except Exception as e:
         logFile.write(HOOFSettings.errorTag + ": unknown error in analysis, stack is:")
         logFile.write(traceback.format_exc())
         print(HOOFSettings.errorTag + ": unknown error in analysis, stack is:")
         print(traceback.format_exc())         
         inHdf.close()
         outHdf.close()
         logFile.close()
         os.remove(outFilePath)
         continue

      finishTime = time.time()

      # print timing to console
      if HOOFSettings.printConsoleTiming:
         tsum = finishTime - beginTime
         print("Time analysis:")
         print("   file read:        ", np.round(t[1]-t[0], 3),   "   ", np.round((t[1]-t[0]) / tsum, 3))
         print("   dataset sort:     ", np.round(t[2]-t[1], 3),   "   ", np.round((t[2]-t[1]) / tsum, 3))
         print("   data validate:    ", np.round(t[3]-t[2], 3),   "   ", np.round((t[3]-t[2]) / tsum, 3))
         print("   data write:       ", np.round(t[4]-t[3], 3),   "   ", np.round((t[4]-t[3]) / tsum, 3))
         if HOOFSettings.dealiasing or HOOFSettings.superobing:         
            print("   data store:       ", np.round(t[5]-t[4], 3),   "   ", np.round((t[5]-t[4]) / tsum, 3))
         if HOOFSettings.dealiasing:
            print("   dealias check:    ", np.round(t[6]-t[5], 3),   "   ", np.round((t[6]-t[5]) / tsum, 3))
            print("   theory calculate: ", np.round(t[7]-t[6], 3),   "   ", np.round((t[7]-t[6]) / tsum, 3))
            print("   height sectors:   ", np.round(t[8]-t[7], 3),   "   ", np.round((t[8]-t[7]) / tsum, 3))
            print("   wind model:       ", np.round(t[9]-t[8], 3),   "   ", np.round((t[9]-t[8]) / tsum, 3))
            print("   dealias:          ", np.round(t[10]-t[9], 3),  "   ", np.round((t[10]-t[9]) / tsum, 3))
            print("   dealias write:    ", np.round(t[11]-t[10], 3), "   ", np.round((t[11]-t[10]) / tsum, 3))
         if HOOFSettings.superobing:
            print("   superob check:    ", np.round(t[12]-t[11], 3), "   ", np.round((t[12]-t[11]) / tsum, 3))
            print("   superob prepare:  ", np.round(t[13]-t[12], 3), "   ", np.round((t[13]-t[12]) / tsum, 3))
            print("   superob:          ", np.round(t[14]-t[13], 3), "   ", np.round((t[14]-t[13]) / tsum, 3))
            print("   file write:       ", np.round(t[15]-t[14], 3), "   ", np.round((t[15]-t[14]) / tsum, 3))

      # if everything is ok, write output and close the files and remove the log file if empty
      endTime = time.time()
      goodFileCounter = goodFileCounter + 1
      print("    Analysis done in", np.round(endTime-beginTime, 1), "seconds.")
      inHdf.close()
      outHdf.close()
      logFile.close()
      if os.path.getsize(logFilePath) == 0:
         os.remove(logFilePath)
 
   endTime = time.time()
   print("HOOF sucessfully analysed", goodFileCounter, "out of", fileCounter, "files in ", np.round(endTime-startTime, 1), "seconds.")
