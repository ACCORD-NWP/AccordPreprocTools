#!/usr/bin/env python
# prepopera.py does a sanity check of the OPERA radar data files. Wrong units are corrected 
# missing values are filled in and so on.
# Each elevation is checked som that it does not overlap with any other elevation.
# Super observations are constructed for both wind and reflectivity.
#
# prepopera.py can handle files from OIFS and ODE.
#
# usage: python prepopera.py -d date (yyyymmddhh) -i indir -o outdir -t variable
#
# The variable can be either refl, wind or comb. The latter is a combinaton of both to be used 
# in the data assimilation.
#
import os
import sys,getopt,string
import numpy as np
from array import *
import h5py
import glob
import time
import random
from multiprocessing import Process, Queue, current_process, freeze_support, cpu_count


# General Settings:
# Specify which countries to process, e.g. [se, dk, fi]
#countries           = ["no"] 
#countries           = ["dk","se","fi","no","ee","fr","uk","ie","de","nl","be","is","pl","es"] 
countries           = ["dk","se","fi","no","ee","fr","uk","ie","de","es","nl","be","pl","pt","is"] 
nocpu=cpu_count()
#NUMBER_OF_PROCESSES = nocpu-2  # Specify the maximum allowed parallel processes used
NUMBER_OF_PROCESSES = 1        # Specify the maximum allowed parallel processes used
datasource          = "nimbus"   # Datasource ode/oifs
newrscale_dbz       = 6000     # Bin size in metres
arclim_dbz          = 6000     # Maximum size in ray direction
newrscale_dow       = 3000     # Bin size in metres
arclim_dow          = 3000     # Maximum size in ray direction
newrscale_comb      = 6000     # Bin size in metres
arclim_comb         = 6000     # Maximum size in ray direction
newrayscale_dbz     = 3        # Azimuth angle in degrees
newrayscale_dow     = 2        # Azimuth angle in degrees
newrayscale_comb    = 3        # Azimuth angle in degrees
restorethresh       = 0.55     # Level of quality to accept
clearsky_dbz        = 0        # Below this level, in dBz, is clear sky
sopercentlimit      = 0.3      # This many pecent of the SO need to contain good observations
windlimit           = 60       # Radial winds more than +-windlimit m/s is disregarded
windstd             = 5        # Maximum variation within the wind SO (m/s)
encoding=('utf-8')             # Encoding for p3 compatibility for conditional checks  
NI_min=30.0                    # Minimum Nyquist velocity that is allowed!


#------------------------
#
# Function run by worker processes
#

def worker(input):
    for func, args in iter(input.get, 'STOP'):
        calculate(func, args)

#
# Function used to calculate result
#

def calculate(func, args):
    func(*args)
    return '%s says that %s%s ' % \
        (current_process().name, func.__name__, args)

#------------------------
def find_attr(fid,attr_name):
# find number of scans!
     number_of_scans=scan_nr(fid)    
     attr_group_list=["how","what","where"]
     result=0
     scan_list=[]
     for scan_number in range(number_of_scans+1):
         if (scan_number==0):
             scan_list.append("") 
         else:
             scan_name="/dataset%d/" % (scan_number)
             scan_list.append(scan_name)
     for path_item in scan_list:
         for group_item in attr_group_list:
             pathname=path_item +  group_item 
             if pathname in fid.keys():
                 attr_local = fid[pathname]
                 for item in attr_local.attrs.keys(): 
                     if (item==attr_name):
                         result = attr_local.attrs[attr_name]
     return result
#------------------------
def find_local_NI(fid,path_obj):
# find number of scans!
#    print(path_obj)
     pathname=path_obj + "/how"
#    print(pathname)
     try:
        obj = fid[pathname]
        result = 0
        for item in obj.attrs.keys(): 
            if (item == "NI"):
               result = 1
     except:
         result = 0
     return result
#------------------------
def check_wmoid(fid,fname):
    try:
        antgainH = find_attr(fid,'antennagain')
    except:
        antgainH = 45.0
    try:
        pulsewidth = find_attr(fid,np.string_("pulsewidth"))
    except:
        pulsewidth = 1.0
    try:
        wavelength = find_attr(fid,np.string_("wavelength"))
        if (wavelength==0):
            wavelength=5.33
        elif ((wavelength>0) and (wavelength < 0.1)):
            wavelength=100.0*wavelength
    except:
        wavelength=5.33
    sensitivity = -108.0
    return antgainH,pulsewidth,wavelength,sensitivity
#------------------------
def copyattr(groupoid,name,groupood):
    global nbinsp,factorp,nrays,nbins
    if (name) not in groupood:
       group=groupood.create_group(name)

    for item in groupoid.attrs.keys():
        if (item == 'rscale'):
#          group.attrs[item] = newrscale
           group.attrs.create(item,newrscale,dtype=np.float64)
#           print("attribute type = ",group.attrs[item].dtype)
        elif (item == 'nbins'):
           factor = groupoid.attrs['rscale']/newrscale
           factorp = int(1/factor)
           nbins = groupoid.attrs[item]
           nbinsp = int(nbins*factor)
           group.attrs[item] = nbinsp
        elif (item == 'nrays'):
           nrays = groupoid.attrs[item]
           group.attrs[item] = groupoid.attrs[item]
        else:
           group.attrs[item] = groupoid.attrs[item]

def oldgrid(groupoid,name):
    global nrays,nbins,rscale
    for item in groupoid.attrs.keys():
        if (item == 'rscale'):
           rscale = groupoid.attrs['rscale']
        if (item == 'nbins'):
           nbins = groupoid.attrs[item]
        if (item == 'nrays'):
           nrays = groupoid.attrs[item]

def getgrid(groupoid,name):
    global nbinsp,factorp,nrays,nbins
    for item in groupoid.attrs.keys():
        if (item == 'nbins'):
           factor = groupoid.attrs['rscale']/newrscale
           factorp = int(1/factor)
           nbins = groupoid.attrs[item]
           nbinsp = int(nbins*factor)
        elif (item == 'nrays'):
           nrays = groupoid.attrs[item]

def newgrid(nbins,nrays,rscale,newrscale,newrayscale):
  # Calculates the new thinner grid to be used by all elevations
  binfactor = int(newrscale/rscale)
  newnbins  = int(nbins/binfactor)

  rayscale  = 360./nrays
  rayfactor = int(round(newrayscale/rayscale))

  newnrays  = int(nrays/rayfactor)

  return newnbins,newnrays,binfactor,rayfactor

def scan_nr(fid):
    nscans = 0
    groupid = fid['/']
    for item in groupid.keys():
        if "dataset" in item: nscans = nscans + 1

#    print("Number of datasets =",nscans)
    return nscans

def nr_quality(groupid):
    quality_fields = 0
    for item in groupid.keys():
        if "quality" in item: quality_fields = quality_fields + 1

    return quality_fields

def nr_datas(groupid):
    datas = 0
    for item in groupid.keys():
        if "data" in item: datas = datas + 1

    return datas

def controlled(fid,datasetnr):
    table = fid[datasetnr]
    for item in table.keys():
        if ('quality' in item): 
            qc_check = 1
            return qc_check 
        else:
            qc_check = 0
    return qc_check

def containsDBZH(fid,datasetnr):
    dbzh_check = 0
    table = fid[datasetnr]
    total_nr_data = nr_datas(table)
    for data_number in range(1,total_nr_data+1):
        datalinkname = datasetnr + "/data%d/what/" % (data_number)
        datalink = fid[datalinkname]
        if (datalink.attrs['quantity'].decode(encoding) == 'DBZH'):
            dbzh_check = 1
            return dbzh_check 
        else:
            dbzh_check = 0
    return dbzh_check


def empty_file(fid):
    table = fid['/']
    total_nr_scans = scan_nr(fid)
    if (total_nr_scans == 0):
        return total_nr_scans
    else:
        return 1

def corrupted(fid):
    table = fid['/']
    total_nr_scans = scan_nr(fid)
    for scan_number in range(1,total_nr_scans+1):
        scan_name = "/dataset%d/" % (scan_number)
        datasetlink = fid[scan_name]
        total_nr_data = nr_datas(datasetlink)
        for data_number in range(1,total_nr_data+1):
            data_name = "/dataset%d/data%d/" % (scan_number,data_number)
            try:
                data_link=fid[data_name]
                data = data_link['data'][()]
            except:
                print ("Unexpected error: dataset # data # is corrupt",scan_number,data_number)
                exit()

def usable_radwind(fid,scanlist):
    windlist = []
    for setnr in (scanlist):
        vrad_found = 0
        grps = "/dataset%d" % (setnr)
        grpi = fid[grps]
        nr_vars = nr_datas(grpi)
        for j in range(nr_vars):
            dsetnm = "/dataset%d/data%d" % ((setnr), (j+1))
            anm    = dsetnm+'/what'
            grp    = fid[anm]
            if (grp.attrs['quantity'].decode(encoding) == 'VRADH'): 
                dsetatt = "/dataset%d/how" % (setnr)
                attr_local = fid[dsetatt]
                attr_data_local = "dataset%d/data%d" % (setnr,j+1)
                vrad_found = 0
                if (find_local_NI(fid,grps) == 1):
                    if (attr_local.attrs["NI"] > NI_min):
#                      print("1: vrad found usable ",attr_local.attrs["NI"],setnr,j)
                       vrad_found = 1
                elif (find_local_NI(fid,attr_data_local) == 1): 
                   attr_data_local = "dataset%d/data%d/how" % (setnr,j+1)
                   attr_local = fid[attr_data_local]
                   if (attr_local.attrs["NI"] > NI_min):
#                     print("2: vrad found usable ",attr_local.attrs["NI"],setnr,j)
                      vrad_found = 1
                else: 
                   print("Did not find NI attribute! Trying to extract it from offset")
                   anm    = dsetnm+'/what'
                   grp    = fid[anm]
                   if (np.abs(grp.attrs["offset"]) >= NI_min):
#                      print("2: vrad found usable ",attr_local.attrs["NI"],setnr,j)
                       print("2: vrad found usable ",setnr,j)
                       vrad_found = 1
        windlist.append(vrad_found)
    return windlist

def overlapping(fid):
   table = fid['how']
   nscans = scan_nr(fid)
   beamwidth = find_attr(fid,'beamwidth')
   if (beamwidth == 0): # find a better check for this!!
       beamwidth = 1.0
   print("There are ",nscans," scans in file with a beamwidth = ",beamwidth)
   groupid=fid['/']
   lowest=90.0
   scanlist=[]
   anglelist=[]
   removedlist=[]

   while ((len(scanlist) < 1) and (len(removedlist) < nscans)): 
       lowest=90.0
       nextscan=100
       for i in range(nscans):
           j=i+1
           if (j not in removedlist):
               grpn = "/dataset%d" % (j) + "/where/"
               wgrpn = fid[grpn]
               elangle = wgrpn.attrs['elangle']
               if (elangle <= lowest): 
                   datasetnr="/dataset%d" % (j)
                   lowest=elangle
                   lscan=j
       qc_checked = 0
       qc_checked = controlled(fid,datasetnr)
       dbzh_checked = 0
       dbzh_checked = containsDBZH(fid,datasetnr)
       if (qc_checked == 0):
           print('Scan not quality controlled. Skipping.')
           removedlist.append(lscan)
       elif ((dbzh_checked == 0 ) and (qc_checked != 0)):    
           print('Scan contains quality field but no dbzh! Skipping')
           removedlist.append(lscan)
       else:
           scanlist.append(lscan)
           anglelist.append(lowest)
   print("scan nr:",lscan,"contains the elevation smallest angle = ",lowest)
   while (len(scanlist)+len(removedlist) < nscans):
       lowest=90.0
       nextscan=100
       qc_checked = 0
       for i in range(nscans):
           j=i+1
           if ((j not in scanlist) and (j not in removedlist)):
               grpn = "/dataset%d" % (j) + "/where/"
               wgrpn = fid[grpn]
               elangle = wgrpn.attrs['elangle']
               if (elangle <= lowest):
                   datasetnr="/dataset%d" % (j)
                   lowest=elangle
                   nextscan=j
       qc_checked = controlled(fid,datasetnr)
       if (qc_checked == 0):
           print('Scan not quality controlled. Skipping.')
           removedlist.append(nextscan)
       else:
           scanlist.append(nextscan)
           anglelist.append(lowest)
             
   print("Order of scans : ",scanlist)
   print("Angles of scans :",["%0.2f" % i for i in anglelist])
   windlist = usable_radwind(fid,scanlist)
   finallist=scanlist
   for i in reversed(range(len(scanlist))):
       if (i>0):
           if (anglelist[i]-anglelist[i-1] < 0.5*beamwidth):
               if (anglelist[i-1] >= beamwidth*0.5):
                   if (windlist[i] == 1):
                       del finallist[i-1]
                       del anglelist[i-1]
                       del windlist[i-1]
                   else: 
                       del finallist[i]
                       del anglelist[i]
                       del windlist[i]
               else:
                   break
   deadlist=finallist
   for i in reversed(range(len(deadlist))):
       if (beamwidth*0.5 > anglelist[i]):
           del finallist[i]
   deadlist=finallist
   for i in reversed(range(len(deadlist))):
       dbzh_checked = 0
       datasetnr="/dataset%d" % (deadlist[i])
       dbzh_checked = containsDBZH(fid,datasetnr)
       if (dbzh_checked == 0):
           del finallist[i]
   print("Final list after removal scheme ")
   print("Order of scans : ",finallist)
   print("Length of final list : ",len(finallist))
   return finallist
           
# ------------------

def createso(data,quality,gain,offset,nodata,undetect,qgain,qoffset,nbins,nrays,newnbins,newnrays,binfactor,rayfactor,restorethresh,clearsky,rscale,bias):
  # All observations below xx dBz is regarded as clear sky
  # For radial wind this is set to zero
#  clearlim = (clearsky - offset)/gain
  clearlim = clearsky

 # Create the new matrices initialized with "nodata"
  nodatatrue = nodata*gain+offset
  dataso    = np.zeros((newnrays,newnbins)) + nodata
  qualityso = np.zeros((newnrays,newnbins)) + nodata
  
  # Step through the rays in steps of rayfactor
  # and bins in steps of bin factor to craete a small tmp matrix
  raycount = -1

  for rray in range (0,nrays,rayfactor):
      if (rray+rayfactor > nrays):
          break

      raycount = raycount + 1
      bincount = -1
      for bbin in range (0,nbins,binfactor):
          if(bbin+binfactor > nbins):
              break

          bincount = bincount + 1

          radius = bbin * rscale
          arclength = (nrays/360)*rayfactor/360 * np.pi * 2 * radius
          if ((arclength//arclim) <= 1):
             rayfactor1 = rayfactor
          elif ((arclength//arclim) <= 2) and ((arclength//arclim) > 1):
             rayfactor1 = rayfactor-1
          else:      
             rayfactor1 = rayfactor-2
          if (rayfactor1 < 1):
              rayfactor = 1
          
          datatmp = data[rray:rray+rayfactor1,bbin:bbin+binfactor]
          qualtmp = quality[rray:rray+rayfactor1,bbin:bbin+binfactor]

          # Calculate how many obspoints gives i.e. 30% of the superobs
          # To be used with count1
          col,row = np.shape(datatmp)
          noobsinso = col*row
          percentage = int(noobsinso*sopercentlimit)          
          
      
          # Create a super obs from this smaller matrix if conditions are fulfilled 
          count1 = 0
          count2 = 0
          count3 = 0
          sotmp  = 0
          for raytmp in range (0,rayfactor1):
              for bintmp in range (0,binfactor):
#                  if (datatmp[raytmp,bintmp]*gain+offset > clearlim) and ((qualtmp[raytmp,bintmp]*qgain+qoffset > restorethresh) or (qualtmp[raytmp,bintmp] == 255)): 
                  if (datatmp[raytmp,bintmp]*gain+offset > clearlim) and (qualtmp[raytmp,bintmp]*qgain+qoffset > restorethresh) and (datatmp[raytmp,bintmp] != nodata): 
                      sotmp  = sotmp + datatmp[raytmp,bintmp]*gain+offset
                      count1 = count1 + 1
#                  elif (datatmp[raytmp,bintmp]*gain+offset <= clearlim) and ((qualtmp[raytmp,bintmp]*qgain+qoffset > restorethresh) or (qualtmp[raytmp,bintmp] == 255)): 
                  elif (datatmp[raytmp,bintmp]*gain+offset <= clearlim) and (qualtmp[raytmp,bintmp]*qgain+qoffset > restorethresh) and (datatmp[raytmp,bintmp] != nodata): 
                      count2 = count2 + 1
                  else:
                      count3 = count3 + 1

          if (count1>percentage):
              # Rainy pixels of good quality
              # need to know the filename and hence get the correct bias value
              # if name is in list then add/subtract scaling value from list using alpha*bias+beta
              if (bias !=0.0):
                  temp1=(sotmp/count1)+bias
#                  temp2=(temp1-offset)/gain
#                 print sotmp/count1,((sotmp/count1)*gain+offset),bias,temp2
#                  if (temp2<0):
#                      temp2=0
#                  elif (temp2>254):
#                      temp2=254
                  dataso[raycount,bincount] = int((temp1-offset)/gain)
#                 dataso[raycount,bincount] = sotmp/count1-(bias+32.0)/0.5
              else:
                  dataso[raycount,bincount] = int(((sotmp/count1)-offset)/gain)
              qualityso[raycount,bincount] = int((0.9-qoffset)/qgain)      # Arbitrary value above restorthresh
          elif(count2>0):
              # Clear pixels of good quality
              dataso[raycount,bincount] = undetect
              qualityso[raycount,bincount] = int((0.9-qoffset)/qgain)    # Arbitrary value above restorthresh
          else:
              # All pixels are of poor quality
              dataso[raycount,bincount] = 10     # Arbitrary value within observation range
              qualityso[raycount,bincount] = int((0.1-qoffset)/qgain)  # Arbitrary value below restorthresh
           
  return dataso,qualityso

# ------------------

def createso_dow(data,gain,offset,nodata,undetect,quality,qgain,qoffset,nbins,nrays,newnbins,newnrays,binfactor,rayfactor,restorethresh,rscale):

 # Create the new matrices initialized with "nodata"
  nodatatrue = nodata*gain+offset
  dataso    = np.zeros((newnrays,newnbins)) + nodata
  qualityso = np.zeros((newnrays,newnbins)) + nodata

  # Step through the rays in steps of rayfactor
  # and bins in steps of bin factor to create a small tmp matrix
  raycount = -1

  for rray in range (0,nrays,rayfactor):
      if (rray+rayfactor > nrays):
          break

      raycount = raycount + 1
      bincount = -1
      for bbin in range (0,nbins,binfactor):
          if(bbin+binfactor > nbins):
              break

          bincount = bincount + 1

          radius = bbin * rscale
          arclength = (nrays/360)*rayfactor/360 * np.pi * 2 * radius
          if ((arclength//arclim) <= 1):
             rayfactor1 = rayfactor
          elif ((arclength//arclim) <= 2) and ((arclength//arclim) > 1):
             rayfactor1 = rayfactor-1
          else:      
             rayfactor1 = rayfactor-2
          if (rayfactor1 < 1):
              rayfactor = 1
          
          datatmp = data[rray:rray+rayfactor1,bbin:bbin+binfactor]
          qualtmp = quality[rray:rray+rayfactor1,bbin:bbin+binfactor]

          # Calculate how many obspoints gives i.e. 30% of the superobs
          # To be used with count1
          col,row = np.shape(datatmp)
          noobsinso = col*row
          percentage = int(noobsinso*sopercentlimit)          
          
      
          # Create a super obs from this smaller matrix if conditions are fulfilled 
          count1 = 0
          count3 = 0
          sotmp  = 0
          for raytmp in range (0,rayfactor1):
              for bintmp in range (0,binfactor):
                  if (abs(datatmp[raytmp,bintmp]*gain+offset) < windlimit) and (datatmp[raytmp,bintmp] != nodata) and (datatmp[raytmp,bintmp] != undetect) and (qualtmp[raytmp,bintmp]*qgain+qoffset > restorethresh): 
                      count1 = count1 + 1
                      if (count1 == 1):
                          sotmp = np.array([datatmp[raytmp,bintmp]*gain+offset])
                      else:
                          sotmp = np.append(sotmp,datatmp[raytmp,bintmp]*gain+offset)
                  else:
                      count3 = count3 + 1
          if (count1>percentage):
              # Enough information to create a superobs
              # Check that the std is not too large
              tmpmean = np.mean(sotmp)
              stdtmp  = np.std(sotmp)
#              print "STD = ",stdtmp, tmpmean, count1
              if (stdtmp < windstd):
                  dataso[raycount,bincount] = int((tmpmean-offset)/gain)
                  qualityso[raycount,bincount] = int((0.9-qoffset)/qgain)      # Arbitrary value above restorethresh, not really used
              else:
                  # Too large std to be used
                  dataso[raycount,bincount] = nodata     # No data
                  qualityso[raycount,bincount] = int((0.1-qoffset)/qgain)  # Arbitrary value below restorethresh, not really used
          else:
              # Too few observations in this subset to be used
              dataso[raycount,bincount] = nodata     # No data
              qualityso[raycount,bincount] = int((0.1-qoffset)/qgain)  # Arbitrary value below restorethresh, not really used
           
  return dataso,qualityso

# ------------------

def convert(filename,dummy):
    print("Now working on file = ",filename)
    dummy_check = 1
    try:
        fid = h5py.File(filename, mode = "r",libver='earliest')
        corrupted(fid)
        fid.close()
    except:
        print("File or parts of file is corrupt and hence will not be processed!") 
        dummy_check = 0
        fid.close()
    if (dummy_check != 0):
       fid = h5py.File(filename, mode = "r",libver='earliest')
       if (empty_file(fid) == 0):
          print("File does not contain any data and hence will not be processed!") 
          dummy_check = 0

    if (dummy_check == 0):
       print("------------------------------------------")
    else:
        fid = h5py.File(filename, mode = "r",libver='earliest')
        if (datasource == "ode"):
            fname=filename[0:2]
        elif (datasource == "oifs"):
            fname=filename[31:33]
        elif (datasource == "nimbus"):
            fname=filename[31:33]
        print(fname)
        table = fid['how']
        if ('nscans' in table.keys()):
           nscans = table.attrs['nscans']
        else:
           nscans = scan_nr(fid)
#       if ( nscans > 1 ):
        finallist=overlapping(fid)
#       else:
#           finallist=[1]
#           print("Only one elevation => No elevation check made!!")

        pfilename=proroot+'/'+'r'+obstyperw+filename
        print('Output file =',pfilename)
        fod = h5py.File(pfilename, mode = "w")
# Reading and storing the root layer attributes
# HOW
# 
        copyattr(table,'how',fod)

# WHERE
        table = fid['where']
        copyattr(table,'where',fod)

# WHAT
        table = fid['what']
        copyattr(table,'what',fod)

# loop through the datasets and break if elevation angle to steep!
        setnr = 0
        count = 0
#       print "number of scans:",len(finallist)
        table = fod['/how']
        table.attrs['nscans']=len(finallist)

#       print len(finallist)
        for i in range(len(finallist)):
           varnr = 0
           setnr = setnr + 1
# new datasetX to outfile
           grpn = "/dataset%d" % (setnr)
           grpo = fod.create_group(grpn)
           grpow = grpn + '/where'
           grpohow = grpn +'/how'
# old datasetX from in file
           grps = "/dataset%d" % (finallist[setnr-1])
#          print "dataset = ",grps
           grpi = fid[grps]
           grpsw = grps + '/where'
           wherei = fid[grpsw]
           oldgrid(wherei,grpow)
           nbinsp,nraysp,binfactor,rayfactor=newgrid(nbins,nrays,rscale,newrscale,newrayscale) 
# here the newgrid parameters should be used!!
           copyattr(wherei,grpow,grpo)
           grpih = grps + '/how'
           grpoh = grpn + '/how'
           grpihh = 'how'
           if grpihh in grpi.keys():
              howi = fid[grpih]
              copyattr(howi,grpoh,grpo) 
              antgainH,pulsewidth,wavelength,sensitivity = check_wmoid(fid,fname)
              fod[grpoh].attrs['sensitivity']=sensitivity
              if "antgainH" not in grpo.attrs:
                  fod[grpoh].attrs['antgainH']=antgainH
              if "pulsewidth" not in grpo.attrs:
                  fod[grpoh].attrs['pulsewidth']=pulsewidth
              if "wavelength" not in grpo.attrs:
                  fod[grpoh].attrs['wavelength']=wavelength
              if "NI" not in howi.attrs:
                  NI = 1
              else:
                  NI = howi.attrs['NI']                  
           else:
               NI = 1
               grouphow = fod.create_group(grpoh)
               antgainH,pulsewidth,wavelength,sensitivity = check_wmoid(fid,fname)
               grouphow.attrs['antgainH'] = antgainH
               grouphow.attrs['pulsewidth']=pulsewidth
               grouphow.attrs['wavelength']=wavelength
               grouphow.attrs['sensitivity']=sensitivity

# put in check of wmoid here and direct to correct fix subroutine!
# datasetX/what
           grpiw = grps + '/what'
           grpow = grpn + '/what'
           whati = fid[grpiw]
           copyattr(whati,grpow,grpo)
           nr_vars = nr_datas(grpi)
           qdata  = np.zeros((int(nrays),int(nbins)))
           pqdata = np.zeros((int(nraysp),int(nbinsp)))
#        This is introduced for those that does not have the latest baltrad qc package!
#        New version of Baltrad QC puts in a qi_total for the total quality index and
#        that should be the default, otherwise it is going the old route with collective 
#        quality index in quality1
           datasetid="/dataset%d/" % (finallist[setnr-1])
           datasetgrp=fid[datasetid]
           quality_fields=nr_quality(datasetgrp)
           for q_nr in range(1,quality_fields+1):
               dsetnqh = "/dataset%d/quality%d/how" % (finallist[setnr-1],q_nr)
#              print 'dsetnqh =',dsetnqh
               qhow=fid[dsetnqh]
               if  np.string_("qi_total") in qhow.attrs['task']:
                   dsetnq = "/dataset%d/quality%d" % (finallist[setnr-1],q_nr)
#                   print("found qi_total in qualityset ",q_nr)
                   break
               else:
                   dsetnq = "/dataset%d/quality1" % (finallist[setnr-1])
#                   print("No qi_total found, using ",dsetnq)
#          print "Path to quality index ",dsetnq
#           dsetnq  = "/dataset%d/quality1" % (finallist[setnr-1])
           grpqno  = "/dataset%d/quality1" % (setnr)  #Will always be quality1
           grpqo   = fod.create_group(grpqno)
           grpqow  = grpqno + '/what'
           grpqoh  = grpqno + '/how'
           grpq    = fid[dsetnq]
           qdata   = grpq['data'][()]
           dsetnqw = dsetnq+"/what"
           dsetnqh = dsetnq+"/how"
           qwhat   = fid[dsetnqw]
           qhow    = fid[dsetnqh]
           copyattr(qwhat,grpqow,grpqo)
           copyattr(qhow,grpqoh,grpqo)
           qoffset = qwhat.attrs['offset']
           qgain   = qwhat.attrs['gain']

           dbzh_found=0
           vrad_found=0
           do_both= None
           for j in range(nr_vars):
              dsetnm = "/dataset%d/data%d" % ((finallist[setnr-1]), (j+1))
              anm    = dsetnm+'/what'
              grp    = fid[anm]
              if (grp.attrs['quantity'].decode(encoding) == 'DBZH'): 
                  dbzh_found=1 
                  any_dbzh_found=1
#              if ((grp.attrs['quantity'] == 'VRAD') and (filename[0:2]=='dk')):
#                  # Readup unambigousvelocity from /how/unambiguousvelocity to NI
#                  grpni=fid['/how']
#                  NI=grpni.attrs['unambiguousvelocity']
#                  alpha_vr=(NI/127.0)
#                  beta_vr=((-128.0/127.0)*NI)
#                  grp.attrs['gain']=alpha_vr
#                  grp.attrs['offset']=beta_vr
              quantity_name=grp.attrs['quantity']
              # Special treatment of radars with erroneous NI
              if fname == "es":
                  NI = float(grp.attrs['gain'])*127.0
#                 NI = 48
#             #
              if fname == "no":
                  NI = 64.0
              #
              if (quantity_name[0:5].decode(encoding) == 'VRADH'): 
                  vrad_found=1 
#                  print 'NI=',float(grp.attrs['gain'])*127.0
           if ((dbzh_found==1) and (vrad_found==1)): do_both=True

           for j in range(nr_vars):
              dsetnm = "/dataset%d/data%d" % ((finallist[setnr-1]), (j+1))
              anm = dsetnm+'/what'
              grp = fid[anm]
              quantity_name=grp.attrs['quantity']
              if ((quantity_name[0:4].decode(encoding) == 'DBZH') or ((quantity_name[0:5].decode(encoding) == 'VRADH') and (do_both))):
#             if (grp.attrs['quantity'] == 'DBZH') or ((grp.attrs['quantity'] == 'VRAD') and (do_both)):
                 varnr = varnr + 1
                 data   = np.zeros((int(nrays),int(nbins)))
                 pdata  = np.zeros((int(nraysp),int(nbinsp)))
                 newname="/dataset%d" % (setnr)
                 grpo=fod ["/"]
# Reading the attributes from the input file
                 anmo = newname+"/data%d" % (varnr) + '/what'
                 copyattr(grp,anmo,grpo)
                 offset=grp.attrs['offset']
                 gain=grp.attrs['gain']
                 nodata=grp.attrs['nodata']
                 undetect=grp.attrs['undetect']
# ----------------------
                 anw = newname + '/where'
                 fod[anw].attrs['nbins']=nbinsp 
                 fod[anw].attrs['nrays']=nraysp
# And now for the data 
                 grpd=fid[dsetnm]
                 data = grpd['data'][()]
                 quantity_name=grp.attrs['quantity']
             # Reflectivity
                 if (quantity_name[0:4].decode(encoding) == 'DBZH'):
#                if (grp.attrs['quantity'] == 'DBZH'):
                     clearsky = clearsky_dbz
                     bias = 0.0
                     pdata,pqdata = createso(data,qdata,gain,offset,nodata,undetect,qgain,qoffset,nbins,nrays,nbinsp,nraysp,binfactor,rayfactor,restorethresh,clearsky,rscale,bias)
             # Doppler Wind
                 if ((quantity_name[0:4].decode(encoding) == 'VRAD') and (do_both)):
#                if ((grp.attrs['quantity'] == 'VRAD') and (do_both)):
                     clearsky = offset
                     pdata,dummydata = createso_dow(data,gain,offset,nodata,undetect,qdata,qgain,qoffset,nbins,nrays,nbinsp,nraysp,binfactor,rayfactor,restorethresh,rscale)            

                 dataname="/dataset%d/data%d" % ((setnr), (varnr))+'/data'
                 fod.create_dataset(dataname,data=np.int_(pdata),chunks=(int(nraysp),int(nbinsp)),compression='gzip',compression_opts=6)

           qualityname="/dataset%d/quality1" % ((setnr))+'/data'
           fod.create_dataset(qualityname,data=np.int_(pqdata),chunks=(int(nraysp),int(nbinsp)),compression='gzip',compression_opts=6)
        fod.close() # buggfix added by mbs@dmi.dk    
        if ((len(finallist) == 0) or (any_dbzh_found==0)): 
            print("pfilename is empty",pfilename)
            os.remove(pfilename)
    print("finished")

if __name__ == "__main__":
    freeze_support()

    any_dbzh_found=0
# Check the input arguments
    newroot = ''
    proroot = ''
    obstyperw = ''
    argv = sys.argv[1:]
    try:
        opts, args = getopt.getopt(argv,"hd:i:o:t:",["date=","indir=","outdir=","obstype="])
    except getopt.GetoptError:
        print('prepopera.py -d <yyyymmdhh> -i <inputdir> -o <outputdir> -t <refl/wind/comb>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('prepopera.py -d <yyyymmddhh> -i <inputdir> -o <outputdir> -t <refl/wind/comb>')
            sys.exit()
        elif opt in ("-d", "--date"):
            sdate = arg
        elif opt in ("-i", "--indir"):
            newroot = arg
        elif opt in ("-o", "--outdir"):
            proroot = arg
        elif opt in ("-t", "--obstype"):
            obstyperw = arg
            if obstyperw != 'refl' and obstyperw != 'wind' and obstyperw != 'comb':
                print('Undefined obstype, must be refl, wind or comb (as in combined)', obstyperw)
                sys.exit()
    print('Date is       :', sdate)
    print('Input dir is  :', newroot)
    print('Output dir is :', proroot)
    print('Obstype is    :', obstyperw)


    if obstyperw == 'refl':
        newrscale      = newrscale_dbz
        arclim         = arclim_dbz
        newrayscale    = newrayscale_dbz
    elif obstyperw == 'wind':
        newrscale      = newrscale_dow
        arclim         = arclim_dow
        newrayscale    = newrayscale_dow
    elif obstyperw == 'comb':
        newrscale      = newrscale_comb
        arclim         = arclim_comb
        newrayscale    = newrayscale_comb

    print('newrscale   : ', newrscale)
    print('arclim      : ', arclim)
    print('newrayscale : ', newrayscale)


    print('Number of available threads: ',nocpu-1)
#    sdate=sys.argv[1]
    YY=sdate[:4] 
    MM=sdate[4:6]
    DD=sdate[6:8]
    HH=sdate[8:10]
    MI=sdate[10:12]
#    proroot= proroot + '/' + YY + '/' + MM + '/' + DD
    if not os.path.exists(proroot):
        os.makedirs(proroot)

    os.chdir(newroot)
    newrscale_store = newrscale

    fcountries=[]
    ffiles=[]
#    for co in countries:
#        fcountries.extend( co )

#    for cnames in fcountries:
    for cnames in countries:
        if (datasource == "ode"):
#            fnames=cnames + '_qcvol_pn129_*' + YY + MM + DD + "T" + HH + MI + '*Z*.h5'
            fnames=cnames + '*_qcvol_*' + YY + MM + DD + "T" + HH + '*Z*.h5'
        elif (datasource == "oifs"):
            fnames="T_PAZZ*_C_EUOC_" + YY + MM + DD + HH + "*00_" + cnames + '*.h5'
        elif (datasource == "nimbus"):
            fnames="T_PAZZ**_C_EUON_" + YY + MM + DD + HH + "*00_" + cnames + '*.h5'
#                   T_PAZZxx_C_EUON_20230601020000_iedub.hdf

        ffiles.extend( glob.glob(fnames) )
        print(fnames)
        length = len(glob.glob(fnames))
        if (length == 0):
            print('Radar ' + cnames + ' does not exist for this date.')

    print(ffiles)
    print(newroot)
    inumb=1
    TASKS = [(convert, (vfile,0)) for vfile in ffiles]
#   print TASKS
    # Create queues
    task_queue = Queue()
# Submit tasks
    for task in TASKS:
        task_queue.put(task)
# Start worker processes
    for i in range(NUMBER_OF_PROCESSES):
        Process(target=worker, args=(task_queue,)).start()
# Tell child processes to stop
    for i in range(NUMBER_OF_PROCESSES):
        task_queue.put('STOP')

