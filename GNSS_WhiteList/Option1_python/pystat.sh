#!/bin/bash
#SBATCH --job-name=stat_calc
#SBATCH --output=guq-%J.out
#SBATCH --error=guq-%J.err

mkdir /lustre/utmp/pns/calc_wl_gnss/t02

python3 read_odbgnss_stats.py
python3 read_stat_stats.py
python3 sort_stats.py
#python3 correct_height_gnsss.py
