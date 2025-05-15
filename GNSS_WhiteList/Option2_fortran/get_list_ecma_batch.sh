#!/bin/bash

yyyy1=2023
mm1=03
dd1=01
hh1=00

yyyy=$yyyy1
mm=$mm1
dd=$dd1

rm -f list_files_aust
for dd in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
do
for hh in 00 03 06 09 12 15 18 21 
do 
 cd /home/tothh/gnss_whitelist
 rm -rf ECMA*

 if [[ -f /home/tothh/ECMA_data/ECMA.${yyyy}${mm}${dd}${hh}.tar ]] 
  then
   rm -f myview${yyyy}${mm}${dd}${hh}.rpt
   cp /home/tothh/ECMA_data/ECMA.${yyyy}${mm}${dd}${hh}.tar .
   tar xvf ECMA.${yyyy}${mm}${dd}${hh}.tar
   cd ECMA
   cp -r ECMA.conv  ECMA_fake_1.conv
   use_odb
   odbdup -i ECMA_fake_1.conv -o ECMA_fake_1_cp.conv/
   cd ECMA_fake_1_cp.conv/
   #odbviewer -q 'select statid,date,time,an_sens_obs,lat,lon,stalt,obsvalue,fg_depar,an_depar,report_status.active@hdr,report_status.passive@hdr,report_status.blacklisted@hdr,report_status.rejected@hdr,datum_status.active@body, datum_status.blacklisted@body, datum_status.passive@body, datum_status.rejected@body from hdr,body,errstat,index' < "."
   odbviewer -q 'select statid,date,time,an_sens_obs,timeslot,lat,lon,stalt,orography,obsvalue,fg_depar,datum_status.active@body,datum_status.blacklisted@body,an_depar from hdr,body,errstat,index,modsurf WHERE (obstype == 1 && varno=128 )' < "."
   cp  myview.rpt /home/tothh/gnss_whitelist/myview${yyyy}${mm}${dd}${hh}.rpt
   cd /home/tothh/gnss_whitelist
   ls myview${yyyy}${mm}${dd}${hh}.rpt >> list_files_aust
  fi
done
done
exit
