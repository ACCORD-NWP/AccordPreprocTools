#!/bin/ksh

echo "compil  parkind1"
gfortran -fPIC -fno-pie -mcmodel=large -c parkind1.F90
echo "compil  gpssol_mod"
gfortran -fPIC -fno-pie -mcmodel=large -c gpssol_mod.F90
echo "compil select_gpssol"
gfortran -fPIC -fno-pie -mcmodel=large -c select_gpssol.F90
echo "linking"
gfortran -fPIC -fno-pie -mcmodel=large -o select_gpssol.x select_gpssol.o gpssol_mod.o
