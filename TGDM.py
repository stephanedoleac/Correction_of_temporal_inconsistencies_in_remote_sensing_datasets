# SPDX-License-Identifier: MIT

import numpy as np
import math
import xarray as xr
import os


class TGDMProcessor:
    """
    TGDM applied to a 3D NetCDF variable with dimensions (time, lat, lon).

    Attributes
    ----------
    dataset : xarray.Dataset
        Input dataset containing the variable to process.
    varname : str
        Name of the variable to which TGDM will be applied.
    window : int
        Window size (in days) for TGDM moving-window computation.
    """

    def __init__(self, input_nc, varname, window=27):
        """
        Initialize the TGDM processor.

        Parameters
        ----------
        input_nc : str or xarray.Dataset
            Path to NetCDF file or an already opened xarray.Dataset.
        varname : str
            Name of the variable to apply TGDM.
        window : int, optional
            Window size in days (default is 27).
        """
        # Load dataset
        if isinstance(input_nc, xr.Dataset):
            self.dataset = input_nc
        else:
            self.dataset = xr.open_dataset(input_nc, decode_times=True)

        self.varname = varname
        self.window = window

        # Ensure the variable exists
        if varname not in self.dataset:
            raise KeyError(f"Variable '{varname}' not found in dataset.")

        self.dataarray = self.dataset[varname]
        self._check_dimensions()
        self._prepare_data()

    def _check_dimensions(self):
        """
        Check that the variable has the expected dimensions.

        Raises
        ------
        ValueError
            If variable does not have 'time' dimension or is not 3D (time, lat, lon).
        """
        if 'time' not in self.dataarray.dims:
            raise ValueError("Input variable must have a 'time' dimension.")

        spatial_dims = [d for d in self.dataarray.dims if d != 'time']
        if len(spatial_dims) != 2:
            raise ValueError(f"Variable must be 3D with dims (time, lat, lon). Found: {spatial_dims}")

        self.lat_dim, self.lon_dim = spatial_dims

    def _prepare_data(self):
        """
        Load data into memory, reshape for computation, and check assumptions.

        Raises
        ------
        ValueError
            If time dimension length is not a multiple of 365 (full years required).
        """
        self.data = self.dataarray.values
        self.T, self.nlat, self.nlon = self.data.shape

        # Ensure daily data without leap years
        if self.T % 365 != 0:
            raise ValueError(
                f"Time length T={self.T} is not a multiple of 365. "
                "TGDM requires full years of 365 days (no leap days)."
            )
        self.n_years = self.T // 365
        self.npix = self.nlat * self.nlon

        # Create a binary mask of valid data (1 = valid, 0 = NaN)
        self.mask_flat = (~np.isnan(self.data)).astype(np.int8).reshape(self.T, self.npix)
        self.data_flat = self.data.reshape(self.T, self.npix)

    def compute_tgdm_mask(self):
        """
        Compute the TGDM mask based on moving-window minimum observations.

        Returns
        -------
        yearly_mask_tiled : np.ndarray
            Array of shape (T, npix) with 1 for valid data and NaN for insufficient observations.
        """
        w = int(self.window)
        if w < 1 or w > self.T:
            raise ValueError("window must be >=1 and <= number of time steps.")

        half_window = math.floor(w / 2)

        # Compute cumulative sum for efficient moving-window sum
        cs = np.vstack([np.zeros((1, self.npix), dtype=np.int32),
                        np.cumsum(self.mask_flat, axis=0, dtype=np.int32)])
        ws = cs[w:, :] - cs[:-w, :]

        # Place summed window into full-length array
        summed_full = np.full((self.T, self.npix), np.nan, dtype=float)
        summed_full[half_window:half_window + ws.shape[0], :] = ws

        # Reshape to (n_years, 365, npix) and compute minimum obs per day-of-year
        summed_reshaped = summed_full.reshape(self.n_years, 365, self.npix)
        min_obs = np.nanmin(summed_reshaped, axis=0)  # shape (365, npix)

        # Build yearly mask: NaN where observations are insufficient
        yearly_mask = np.where(min_obs == 0, np.nan, 1.0)
        yearly_mask_tiled = np.tile(yearly_mask, (self.n_years, 1))

        return yearly_mask_tiled

    def apply_tgdm(self, output_varname=None):
        """
        Apply TGDM to the variable and return a corrected xarray.DataArray.

        Parameters
        ----------
        output_varname : str, optional
            Name of the corrected variable. Defaults to varname + "_tgdm".

        Returns
        -------
        corrected_da : xarray.DataArray
            Corrected data array with TGDM applied.
        """
        mask = self.compute_tgdm_mask()
        corrected_flat = self.data_flat * mask
        corrected = corrected_flat.reshape(self.T, self.nlat, self.nlon)

        outname = output_varname or (self.varname + "_tgdm")
        corrected_da = xr.DataArray(
            corrected,
            coords=self.dataarray.coords,
            dims=self.dataarray.dims,
            attrs=self.dataarray.attrs,
            name=outname
        )
        return corrected_da

    def save_corrected(self, output_nc, overwrite=False, output_varname=None):
        """
        Apply TGDM and save the corrected dataset to a NetCDF file.

        Parameters
        ----------
        output_nc : str
            Path to save the corrected NetCDF file.
        overwrite : bool, optional
            If True, overwrite existing file (default False).
        output_varname : str, optional
            Name of the corrected variable in output dataset.

        Returns
        -------
        ds_out : xarray.Dataset
            Dataset containing the corrected variable.
        """
        corrected_da = self.apply_tgdm(output_varname=output_varname)
        ds_out = self.dataset.copy()
        ds_out[corrected_da.name] = corrected_da

        if os.path.exists(output_nc) and overwrite:
            os.remove(output_nc)

        ds_out.to_netcdf(output_nc)
        print(f"Wrote corrected variable '{corrected_da.name}' to {output_nc}")
        return ds_out

