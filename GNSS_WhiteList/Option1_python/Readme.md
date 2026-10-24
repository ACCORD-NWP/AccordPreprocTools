date final code version 3-12-2024
date update 01-04-2025
update 2026-01-06

Code initially creted by Jana Sanchez (AEMET) and converted to python by Jose Miguel Perez de Garcia Delgado (AEMET).
Further unified by Benedikt Strajnar (ARSO).

The procedure currently includes 3 steps. The step 1 depends on the availablity of tools to query the ODB (mandalay, odbsql, odbviewer, ...)  

## 1) Produce an archive of extracted ODB using all GNSS ZTD observations (ECMA data base).
 
This step depends on the availablity of tools to query the ODB (mandalay, odbsql, odbviewer, ...) and the organization of ODB archive. 
An example using odbsql and ECMA ODB bases from ARSO archive is provided here.

scripts: extract_odb_data.py 
ODB query: select_gnss.sql
inputs: an archive with ODB databases (from a passive experiment)
outputs: a directory with ODB extracts in text format 

## 2) Create a list with metadata for all stations

This step depends on the use of OBSOUL or BUFR. For BUFR, it might be easier to construct the metadata file from the ODB archive directly.   

script: readgnssdatas.py 
inputs: a folder with ascii files like YYMMDDHH_gnssubh_v1_0020 (2023120722_gnssubh_v1_0020) of a selected period. 
outputs: a file with a list of stations in ascii (list20)

Instead, a simple way of doing this step is also to create is from the extrated ODB data
'''
awk 'NR > 1 && !seen[$1]++' <ODB_filename_prefix>* > unique_list.txt
'''
where <ODB_filename_prefix> is the initial string of the extracted ODB files.

## 3) Computation of the whitelist

This part includes gathering of departures per GNSS station, calculcation of statistics per station, data selection and optional thinning. 

script: compute_gnss_whitelist.py
inputs: ODB extracted data from step 1, or optionally, data per station if already available (e.g. to repeat the whitelist generation with different thinning), metadata file  
outputs: stationlist

