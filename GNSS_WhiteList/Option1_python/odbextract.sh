#!/bin/bash
#SBATCH --job-name=gnssodb_untar_queue
#SBATCH --output=guq-%J.out
#SBATCH --error=guq-%J.err

#this script copy and decompress the odb of an experiment (code activates or deactivate with sw01)
#then make a sql queue to obtain the information of the stations that actually have been included
#finally concatenate the different queues keeping only one entry for each id and sort it at the same time


module load intel/2020
module load impi
module load odb
module load odb_api

wdir='/lustre/utmp/pns/odbgnss/exp1'
mkdir $wdir
ori='/MASIVO/pns/harmonie/AIBR_wl_gnss01'


sw01=no

if [ $sw01 == yes ] ; then

for yy in 2023
do
mkdir ${wdir}/${yy}/

for mm in 12
do
mkdir ${wdir}/${yy}/${mm}/


#for dd in 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
for dd in 01 02 03 04 05 06 07
do
mkdir ${wdir}/${yy}/${mm}/${dd}


for hh in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 
do 
mkdir ${wdir}/${yy}/${mm}/${dd}/${hh}

cd ${wdir}/${yy}/${mm}/${dd}/${hh}

cp ${ori}/${yy}/${mm}/${dd}/${hh}/odb_stuff.tar .
tar -xvf odb_stuff.tar
tar -xvf odbvar.tar



done

done
done 
done

fi


sw02=no

if [ $sw02 == yes ] ; then

for yy in 2023
do
#mkdir ${wdir}/${yy}/

for mm in 11 
do
#mkdir ${wdir}/${yy}/${mm}/


for dd in   27
#for dd in   15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 
#for dd in   01 02 03 04 05 06 07 
do
#mkdir ${wdir}/${yy}/${mm}/${dd}


for hh in 18 
#for hh in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 
do 
#mkdir ${wdir}/${yy}/${mm}/${dd}/${hh}

cd ${wdir}/${yy}/${mm}/${dd}/${hh}/odbvar/ECMA.conv/


export ODB_CMA=ECMA

#while [ ! -s "$wdir/gnss${yy}${mm}${dd}${hh}" ]; do


odbsql -v /pred/pns/bin/select_gnss.sql >> $wdir/gnss${yy}${mm}${dd}${hh}

#done

done

done
done 
done

fi


cd ${wdir}
tail -q -n +2  gnss* | awk '!seen[$1]++' | sort > listuniq #remove first line


