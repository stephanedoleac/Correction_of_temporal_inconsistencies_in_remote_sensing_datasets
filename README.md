This repository contains the codes used in our paper "North-South asymmetry in subtropical phytoplankton response to recent warming". The codes for completion are available at https://github.com/Sma6500/Chl_global_completion

- TGDM.py : Implement the Temporal Gap Detection Method (van Oostende et al. 2022)
- correction_calibration_inconsistencies.py : Correct the calibration inconsistencies based on SeaWiFS and MODIS-AQUA time series

# Installation
Required modules :
- numpy (https://numpy.org/)
- math (https://docs.python.org/3/library/math.html)
- xarray (https://docs.xarray.dev/en/stable/)
- os (https://docs.python.org/3/library/os.html)
- tqdm (https://tqdm.github.io/)
- netCDF4 (https://unidata.github.io/netcdf4-python/)
- pathlib (https://docs.python.org/3/library/pathlib.html)

Tested versions :
- v3.9.7

You can install those packages in a few minutes either with ```pip3``` or with ```conda```.

# Usage
### Imput data
Both TGDM.py and correction_calibration_inconsistencies.py work with chlorophyll-a remote sensing datasets with a daily frequency. The datasets used in our study "North-South asymmetry in subtropical phytoplankton response to recent warming" can be downloaded at the following adresses :
- OC-CCI : https://data.marine.copernicus.eu/product/OCEANCOLOUR_GLO_BGC_L3_MY_009_107
- GlobColour-CMEMS : https://data.marine.copernicus.eu/product/OCEANCOLOUR_GLO_BGC_L3_MY_009_103
- GlobColour-AVW : ftp://ftp.hermes.acri.fr

Additionnally, correction_calibration_inconsistencies.py requires chlorophyll-a data from SeaWiFs and MODIS-AQUA with a daily frequency. These datasets can be downloaded from the ACRI-ST ftp server (ftp://ftp.hermes.acri.fr)


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

### Expected output
The corrected datasets analysed in the study are available at https://doi.org/10.5281/zenodo.17122287

TGDM.py and correction_calibration_inconsistencies.py both take a few minutes to run. The completion takes a few hours.


# References
van Oostende, M., Hieronymi, M., Krasemann, H., Baschek, B. & Röttgers, R. Correction of inter-mission inconsistencies in merged ocean colour satellite data. Front. Remote Sens. 3, (2022). 


