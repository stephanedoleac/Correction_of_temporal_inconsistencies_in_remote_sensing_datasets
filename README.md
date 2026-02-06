# Overview

This repository contains the codes used in our paper "North-South asymmetry in subtropical phytoplankton response to recent warming". The codes for completion are available at https://github.com/Sma6500/Chl_global_completion

- TGDM.py : Implement the Temporal Gap Detection Method (van Oostende et al. 2022)
- correction_calibration_inconsistencies_monthly.py : Correct the calibration inconsistencies based on SeaWiFS and MODIS-AQUA time series for monthly datasets
- correction_calibration_inconsistencies_daily.py : Correct the calibration inconsistencies based on SeaWiFS and MODIS-AQUA time series for daily datasets
- data_formatting_for_LFCA.py : Format datasets to perform a low-frequency component analysis (LFCA) with the NCSTAT software

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

Information about NCSTAT installation can be found on the following page : https://pagesperso.locean-ipsl.upmc.fr/terray/ncstat/index.html

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

### Usage of correction_calibration_inconsistencies_daily.py
```
from correction_calibration_inconsistencies import BiasCorrector_daily

# Initialize with neighborhood width (default=5)
corrector = BiasCorrector_daily(width=5)

# Correct a dataset
corrector.correct_dataset(
    merged_path="/path/to/merged.nc",
    modis_path="/path/to/modis.nc",
    seawifs_path="/path/to/seawifs.nc",
    output_path="/path/to/output.nc"
)
```

### Usage of correction_calibration_inconsistencies_monthly.py
```
from correction_calibration_inconsistencies import BiasCorrector_monthly

# Initialize with neighborhood width (default=5)
corrector = BiasCorrector_monthly(width=5)

# Correct a dataset
corrector.correct_dataset(
    merged_path="/path/to/merged.nc",
    modis_path="/path/to/modis.nc",
    seawifs_path="/path/to/seawifs.nc",
    output_path="/path/to/output.nc"
)
```

### Usage of data_formatting_for_LFCA.py
The code data_formatting_for_LFCA.py can be used to format datasets and build a meshmask to perform a LFCA with NCSTAT. The LFCA can then be done with the two following bash commands :
```
comp_eof_miss_3d -f=path_to_formatted_data -v=var_name -m=path_to_meshmask -n=30 -o=output_path_for_eof
comp_loess_rot_eof_3d -f=path_to_eof -v=var_name -m=path_to_meshmask -nt=96 -itdeg=1 -se=1:10 -o=output_path_for_lfc
```

### Expected output
The corrected datasets analysed in the study are available at https://doi.org/10.5281/zenodo.17122287

TGDM.py, correction_calibration_inconsistencies_daily.py and correction_calibration_inconsistencies_monthly.py take a few minutes to run. The completion takes a few hours.

data_formatting_for_LFCA.py takes a few seconds to run, and the subsequent LFCA can take up to a few hours.


# References
van Oostende, M., Hieronymi, M., Krasemann, H., Baschek, B. & Röttgers, R. Correction of inter-mission inconsistencies in merged ocean colour satellite data. Front. Remote Sens. 3, (2022). 

# License

This software is released under the MIT License. See the [LICENSE](LICENSE) file for details.

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

