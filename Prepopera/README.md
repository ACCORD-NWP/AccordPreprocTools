Preprocessing of Radar data from NIMBUS 
may 2025

Mats Dahlbom DMI, Martin Ridal SMHI

Documentation: M. Ridal and M. Dahlbom: Assimilation of Multinational Radar Reflectivity Data in a Mesoscale Model: A Proof of Concept, doi: 10.1175/JAMC-D-16-0247.1

General Description of the Script:
Usage:
(python3) prepopera.py -d $YYYY$MM$DD$HH -i /path/to/your/input/directory -o /path/to/your/output/directory -t outtype (comb/wind/refl) 
or simply run prepopera.py -h and it will give the same result!


Prepopera was originally (> 10 years ago) written in order to make it possible to run opera volume radar files in hdf5 format from different countries and make them homogeneus and possible to pass the data ingestion in Harmonie-Arome.  


Prepopera does three different things:
1. Sanity check
Several checks are made in order to sort out corrupt, empty files and files that does not comply to the given standard, that would be problematic to pass to the observation ingestion of Harmonie-Arome.

2. Elevation overlap check
In order to avoid using the same volume of air twice in the data assimilation a check for overlapping elevations are made. The elevation separation is compared to the beamwith and if it is too small the elevation will be rejected.

3. Super observations
The creation of super observations (SO) is described in the article cited above as well as in the Salonen et al. (to be provided by Martin). The size of the SO is set by the user in the beginning of the script. 


There are three different outtypes (-t). These will be differently treated in the Harmonie-Arome data assimilation. The combined solution (comb) will use the winds anr reflectivity from the same hdf5-files with the same SO size. In order to be able to use different SO size and resolution of reflectivity and winds the refl and wind alternatives were introduced. For refl only reflecitity will be included in the resulting hdf5-files while for the wind option both reflectivity and winds will be included.


Three thinks to keep in mind! 

Number of "cpus" to use in the preprocessing! See nocpu!
The list of countries to be preprocessed! see countries!
Sizes of superobservations! look at the settings for various types of newrscale

