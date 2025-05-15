### How to use this package:

The input would be the odb files and then use mandalay.sql example to get the ascci files from them.
Example of input then would be data/myview2024020600.rpt 

Chose your settings in gpssol_mod.F90 

Compile the package using cmp.sh, it compiles:
parkind1.F90
gpssol_mod.F90
select_gpssol.F90


Then run the executable obtained called select_gpssol.x and you will get the white list.
