# AccordPreprocTools
Gathering various existing observation preprocessing tools and document them.


List of tools:

1) HOOF radar preprocessing software
   
2) prepopera radar preprocessing software
   
3) GNSS White List procedure


Candidates:

 Thinning tool for EMADDC Mode-S EHS



...

## Python package (HOOF and Prepopera)

HOOF and Prepopera are packaged as `accordpreproctools`, installable with e.g. poetry

The install is equivalent to running the original scripts directly:

```
hoof <namelist_file> <input_folder> <output_folder>
prepopera -d <yyyymmddhh> -i <input_dir> -o <output_dir> -t <refl/wind/comb>
```

See `HOOF/README.md` and `Prepopera/README.md` for details on each tool. `HOOF/HOOF_namelist.nam` remains as an example namelist to pass to `hoof`.
