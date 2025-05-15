# Preprocessing of Radar data from NIMBUS 
may 2025

Mats Dahlbom DMI, Martin Ridal SMHI

## Documentation: 
M. Ridal and M. Dahlbom: Assimilation of Multinational Radar Reflectivity Data in a Mesoscale Model: A Proof of Concept, doi: 10.1175/JAMC-D-16-0247.1

## General Description of the Script:

**Usage:**

You need python3.

prepopera.py -d $YYYY$MM$DD$HH -i /path/to/your/input/directory -o /path/to/your/output/directory -t outtype (comb/wind/refl) 

or simply run prepopera.py -h and it will give the same result!

Prepopera was originally (> 10 years ago) written in order to make it possible to run opera volume radar files in hdf5 format from different countries and make them homogeneus and possible to pass the data ingestion in Harmonie-Arome.  

**Prepopera does three different things:**
1. Sanity check
Several checks are made in order to sort out corrupt, empty files and files that does not comply to the given standard, that would be problematic to pass to the observation ingestion of Harmonie-Arome.

2. Elevation overlap check
In order to avoid using the same volume of air twice in the data assimilation a check for overlapping elevations are made. The elevation separation is compared to the beamwith and if it is too small the elevation will be rejected.

3. Super observations
The creation of super observations (SO) is described in the article cited above as well as in the Salonen et al. (to be provided by Martin). The size of the SO is set by the user in the beginning of the script. 

There are three different outtypes (-t). These will be differently treated in the Harmonie-Arome data assimilation. The combined solution (comb) will use the winds and reflectivity from the same hdf5-files with the same SO size. In order to be able to use different SO size and resolution of reflectivity and winds the refl and wind alternatives were introduced. For refl only reflecitity will be included in the resulting hdf5-files while for the wind option both reflectivity and winds will be included.

**User defined parameters set in the beginning of the script:**

countries           = ["dk","se","fi","no","ee","fr","uk","ie","de","es","nl","be","pl","pt","is"] # Countries that will be included in the output (if available)

nocpu=cpu_count()              # Count the number of available cpus

NUMBER_OF_PROCESSES = 1        # Specify the maximum allowed parallel processes used. Note: No more than nocpu-2.

datasource          = "nimbus"   # Datasource ode/oifs/nimbus. Needed since the filenames are different

newrscale_dbz       = 6000     # Bin size in metres for the reflectivity SO (m)

arclim_dbz          = 6000     # Maximum size in ray direction for the reflectivity SO (m)

newrscale_dow       = 3000     # Bin size in metres for	the wind SO (m)

arclim_dow          = 3000     # Maximum size in ray direction for the wind SO (m)

newrscale_comb      = 6000     # Bin size in metres for	the combined SO (m)

arclim_comb         = 6000     # Maximum size in ray direction for the combined SO (m)

newrayscale_dbz     = 3        # Azimuth angle in degrees for the reflectivity SO

newrayscale_dow     = 2        # Azimuth angle in degrees for the wind SO

newrayscale_comb    = 3        # Azimuth angle in degrees for the combined SO

restorethresh       = 0.55     # Level of limit of the quality index to accept

clearsky_dbz        = 0        # Below this level, in dBz, is regarded as clear sky

sopercentlimit      = 0.3      # This fration of the SO need to contain good observations

windlimit           = 60       # Radial winds more than +-windlimit m/s is disregarded

windstd             = 5        # Maximum variation within the wind SO (m/s)

encoding=('utf-8')             # Encoding for p3 compatibility for conditional checks

NI_min=30.0                    # Minimum Nyquist velocity that is allowed!

