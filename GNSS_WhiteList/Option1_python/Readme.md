date final code version 3-12-2024
date update 01-04-2025

Code creted by Jana Sanchez AEMET and converted to python by Jose Miguel Perez de Garcia Delgado. AEMET


The white list has been generated following the next four steps:
-Generating a list with all the stations available in a period
scripts: readgnssdatas.py
inputs: a folder with ascii files like YYMMDDHH_gnssubh_v1_0020 (2023120722_gnssubh_v1_0020) of a selected period.
outputs: a file with a list of stations in ascii (list20)

-Running a experiment for a time that allows to perform statistics with the list. 

-Extracting all the observations that had been included in the run of an experiment, the first loop extract the odbvar files, and the loop below it makes queue to the odb. And makes a list containing only the different stations included in the period .
scripts: odbextract.sh select_gnss.sql
inputs: a folder with the database (for example /MASIVO/pns/harmonie/AIBR_wl_gnss01) 
outputs: a list with the stations that entered on a period (/lustre/utmp/pns/odbgnss/exp1/listuniq), and the result of the database queue (/lustre/utmp/pns/odbgnss/exp1/gnss2023111500)

-Calculate the sd and skewness for the innovations of the stations and making a list with the ones with lesser value.
scripts:  read_odbgnss_stats.py read_stat_stats.py sort_stats.py
pystat.sh call the python scripts
read_odbgnss_stats.py generates a file with the innovations of an id.
read_stat_stats.py reads the previous file and calculate the stats for each station
sort_stats.py make the final list, picking a station when there are more that one processing center.
inputs: a list with the stations that entered on a period (/lustre/utmp/pns/odbgnss/exp1/listuniq), and the result of the database queue (/lustre/utmp/pns/odbgnss/exp1/gnss2023111500)
outputs: the final white list (/lustre/utmp/pns/calc_wl_gnss/t02/wlist)


