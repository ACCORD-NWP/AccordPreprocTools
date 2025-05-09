# HOOF.py - Homogenization Of Opera Files

Tool for homogenizing Opera HDF5 radar files as preparation for usage in BATOR.

### **Usage:**

```
./HOOF.py namelist_file input_folder output_folder
```

### **Namelist:**
- **FileExtensions** - files with these extensions are homogenized.
- **SavedQuantities** - all possible names for DBZ, TH and VRAD quantities.
- **DbzQualityGroups** - quality groups (one or more of ROPO, SAT, BLOCK, TOTAL) that have to be present.
- **RadarAttributes common** - list of radar attributes which will be written to the output file and their default values.
Attributes belonging to a dataset, dataset/data or dataset/quality will be written for all datasets (elevations), data and quality groups.
A default value of `None` indicates that this attribute must be present in the input file, as it has no default value.
- **RadarAttributes NOD** - a list of radar attributes (same as above), specific to radar with the specified NOD (opera site identification)
- to comment a line in the namelist, use # as the first character in the line.

### **Actions:**
-  Opera files are 15 minutes aggregates, and most radars have measurements every 5 minutes, so there are usually 3 different measurements contained in one Opera file.
The homogenization tool separates measurements and either outputs one HDF file for each measurement or retains all measurements in a single HDF5 file.
- For each measurement, we rearrange the content according to the user specifications, given in a namelist.
- Validation of the required output meta-data is performed.
- The output file is written.
- The radar quantities (*/dataset/data/what/quantity*) that are written in the output are DBZ, TH and VRAD. 

