WHITE LIST PROCEDURE FOR GNSS ZTD

These scripts are to create a White List of GNSS observations taking into account the statistics using the innovations that come from odb files.


We have two options here:

Option1_python: 
It is a sofware developed in AEMET by Jana Sanchez and wrote into python by Jose Miguel Perez de Gracia . AEMET. april 2025

It chooses the best processing centre for each station according to these rules:
1)The obs with less % of rejections is chosen,
2)then the one whith lowest  SD and if necesarry because same value of SD,
3)then then the one with lowest skewness.
4) It does thinning

Option2_fortran: 
It is a software deleveloped by Pau Poli, MF, and it is being used by Florian Meier LACE, Geosphere Austria. April 2025.

It chooses the best processing centre for each station according to these rules:
1) The obs with stable position (coordinates) and little difference bewten height of model and station
2) tation with enough availability during the test period
3) the station whith lowest  SD 
4) Not too big bias
5) Gaussianity of the innovations
6)It does thinning

Documentation: Poli et al 2007, Yan et al 2008.

In both methods the input are the odb files results of running an experiment with gnss ztd passive (blacklisted:  (src/blacklist/ms_blacklist.b : if (VARIAB= apdds) then fail (EXPERIMENTAL) endif).
 and using all the combination STATPPCC (STAT=station, PPCC=processing centre) that are possible in your domain.
For this:
- you would have to comment in your exp where it ask to have a White list so to check if your station is there because you want all to enter this time so to have its statistics
- and would have to avoid any thinning  (src/arpifs/prep_obs/redgps.F90)
It would be good to have odb files from at least 15 days.

</pre>
