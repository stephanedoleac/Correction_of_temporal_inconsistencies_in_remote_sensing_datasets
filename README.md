This repository contains the codes used in our paper "North-South asymmetry in subtropical phytoplankton response to warming". The codes for completion are available at https://github.com/Sma6500/Chl_global_completion

- TGDM.py : Implement the Temporal Gap Detection Method (van Oostende et al. 2022)
- correction_calibration_inconsistencies.py : Correct the calibration inconsistencies based on SeaWiFS and MODIS-AQUA time series

The corrected datasets analysed in the study are available at https://doi.org/10.5281/zenodo.17122287

### Usage of TGDM.py
```
from TGDM import TGDMProcessor

# Initialize processor
tgdm = TGDMProcessor("input.nc", "varname", window=27)

# Apply TGDM and get corrected DataArray
corrected_da = tgdm.apply_tgdm()

# Or save directly to NetCDF
tgdm.save_corrected("corrected.nc", overwrite=True)
```


### Usage of correction_calibration_inconsistencies.py
```
from correction_calibration_inconsistencies import BiasCorrector

# Initialize with neighborhood width (default=5)
corrector = BiasCorrector(width=5)

# Correct a dataset
corrector.correct_dataset(
    merged_path="/path/to/merged.nc",
    modis_path="/path/to/modis.nc",
    seawifs_path="/path/to/seawifs.nc",
    output_path="/path/to/output.nc"
)
```




### References
van Oostende, M., Hieronymi, M., Krasemann, H., Baschek, B. & Röttgers, R. Correction of inter-mission inconsistencies in merged ocean colour satellite data. Front. Remote Sens. 3, (2022). 


