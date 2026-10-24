#!/usr/bin/env python3

import datetime
import os
import glob
import sys
from contextlib import contextmanager
import os

@contextmanager
def chdir(path):
    old = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old)

# get relevant ecma_full from odb_ecma.tar.gz for a period.
archpath="/home/a_si/oper/results/aosruc04ec/odb_assim"

# 3-hourly input, hourly output
startdate = "2025-07-01"
enddate = "2025-07-02"
odbrequestfile = "select_gnss.sql"
topdir=os.getcwd()
workdir=topdir + "/extracted_odb_data"
os.system("mkdir -p " + workdir)

def execute_sys(command):
  print("Executing: " + command)
  status = os.system(command)
  if status != 0:
    print("Command failed or not found.")
    sys.exit(1)

time = datetime.datetime.strptime(startdate,"%Y-%m-%d")
endtime = datetime.datetime.strptime(enddate,"%Y-%m-%d")

while time <= endtime:
  print(time)
  refdate=datetime.datetime.strftime(time,"%Y%m%d%H")
  ryear=datetime.datetime.strftime(time,"%Y")
  rmonth=datetime.datetime.strftime(time,"%m")
  rday=datetime.datetime.strftime(time,"%d")
  rys=datetime.datetime.strftime(time,"%y")
  rhour=datetime.datetime.strftime(time,"%H")
  julian=datetime.datetime.strftime(time,"%j")
  rminute=datetime.datetime.strftime(time,"%M")
  rpath=datetime.datetime.strftime(time,"/%Y/%m/%d/%H/")
  rhourmin=datetime.datetime.strftime(time,"%H%M")
  refdatum=datetime.datetime.strftime(time,"%Y%m%d")
  rdatetime=datetime.datetime.strftime(time,"%Y%m%d-%H%M")  

  execute_sys("cp " + archpath + rpath + "odb_ECMA_3dvar.tar.gz " + workdir + "/")

  with chdir(workdir):
    execute_sys("tar xvzf odb_ECMA_3dvar.tar.gz ECMA.gpssol")
    with chdir(workdir + "/ECMA.gpssol"):
      execute_sys("ln -sf " + topdir + "/IOASSIGN.ECMA IOASSIGN")
      sqlcommand = "odbsql -v " + topdir + "/" + odbrequestfile + " -k > ../ecma_" + refdate + ".dat"
      execute_sys(sqlcommand) 
    execute_sys("rm -r ./ECMA.gpssol")
    execute_sys("rm odb_ECMA_3dvar.tar.gz")

  time = time + datetime.timedelta(seconds=3600*3)
