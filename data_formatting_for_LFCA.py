# SPDX-License-Identifier: MIT

                            # -------------------------------------- #
                            #   Formatting data before performing    #
                            #        an LFCA with NCSTAT             #
                            # -------------------------------------- #

import numpy as np
import xarray as xr
import netCDF4 as nc
from tqdm import tqdm

### Paths
path_input_data = ""           # Path to input dataset
path_output_data = ""          # Path to output dataset
path_mask_ocean_case1 = ""     # Path to mask of case 1 waters
path_output_meshmask = ""      # Path to output meshmask

###############################################################################################
### Data formatting 
###############################################################################################

### Load data
dataset = xr.open_dataset(path_input_data)
dim_time = len(dataset.time.data)

mask_ocean_case1 = np.load(path_mask_ocean_case1)  # Load mask of case 1 waters

### Compute anomalies
data_noseasonal = dataset.CHL.groupby("time.month") - dataset.CHL.groupby("time.month").mean(dim="time")

### Remove pixels with less than 70 % of values
data_filtered = np.empty(np.shape(data_noseasonal.data))

for t in tqdm(range(dim_time)) :
    data_filtered[t,:,:] = np.where((data_noseasonal.count(dim="time").data < 0.7*dim_time) | (mask_ocean_case1 == False), np.nan, data_noseasonal.data[t,:,:])

### Remove outliers        
data_nooutliers_1 = np.where(data_filtered > np.nanquantile(data_filtered, q=0.999), np.nan, data_filtered)
data_nooutliers = np.where(data_filtered < np.nanquantile(data_filtered, q=0.001), np.nan, data_nooutliers_1)

### Put missing data to the correct format
data_missingvalue = np.where(np.isnan(data_nooutliers), np.float32(-999), data_nooutliers) 

### Make xr.dataset and Netcdf file
new_nc_file = xr.Dataset(data_vars=dict(CHL=(["time", "lat", "lon"], np.array(data_missingvalue, dtype='f'), dict(missing_value=np.float32(-999),
                                                                                                                  _FillValue=np.float32(-999)))),
                         coords=dict(time=np.array(dataset.time.data, dtype='f'),
                                     lat=np.array(np.arange(0,180,1), dtype='f'),
                                     lon=np.array(np.arange(0,360,1), dtype='f'))
                        )

new_nc_file.to_netcdf(path_output_data)


###############################################################################################
### Mesh mask creation 
###############################################################################################

lon = np.arange(-179.5,180,1)
lat = np.arange(-89.5,90,1)
width_gridcell = 40075000/360

e1n = np.ones([len(lat), len(lon)])*width_gridcell

e2n = []
for j in range(360) :
    e2n.append(np.cos(np.deg2rad(lat))*width_gridcell)

e2n = np.transpose(e2n)

new_mask = np.where(np.sum(np.isnan(data_nooutliers), axis=0) == dim_time, 0, 1)

mesh_mask = xr.Dataset(data_vars=dict(nmask=(["lat", "lon"], np.array(new_mask, dtype='f'), dict(valid_min=np.float32(0),
                                                                                                 valid_max=np.float32(1),
                                                                                                 lon1_Eastern_limit=np.float32(0),
                                                                                                 lon2_Western_limit=np.float32(360),
                                                                                                 lat1_Southern_limit=np.float32(0),
                                                                                                 lat2_Northern_limit=np.float32(180))),
                                      e1n=(["lat", "lon"], np.array(e1n, dtype="f"), dict(units="meters",
                                                                                          missing_value=np.float32(1))),
                                      e2n=(["lat", "lon"], np.array(e2n, dtype="f"), dict(units="meters",
                                                                                          missing_value=np.float32(1)))),
                       coords=dict(lat=np.array(np.arange(0,180,1), dtype='f'),
                                   lon=np.array(np.arange(0,360,1), dtype='f'))
)
                      
mesh_mask.to_netcdf(path_output_meshmask)
