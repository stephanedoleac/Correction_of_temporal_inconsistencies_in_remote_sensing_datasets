# SPDX-License-Identifier: MIT

import numpy as np
import xarray as xr
from tqdm import tqdm
import netCDF4 as nc
from pathlib import Path


class BiasCorrector:
    def __init__(self, width=5):
        """
        Initialize the BiasCorrector.

        Parameters
        ----------
        width : int
            Neighborhood size for local averaging when computing local means.
        """
        self.width = width

    def get_local_mean(self, data, i, j, lat):
        """
        Compute the local mean time series in a square neighborhood centered at (i,j).

        Parameters
        ----------
        data : ndarray (time, lat, lon)
            Input dataset.
        i : int
            Latitude index.
        j : int
            Longitude index.
        lat : ndarray
            Array of latitudes.

        Returns
        -------
        ndarray
            Local mean time series (length = time dimension).
            Returns NaN if the region is dominated by land or very high latitudes (> 80°).
        """
        n = len(data[:, 0, 0])
        if np.abs(lat[i]) > 80:
            return np.empty([n]) * np.nan  # Ignore polar regions

        n_continents = 0
        temp = []
        width = self.width

        # Handle longitude wrapping near 0° and 360°
        if j < width // 2:
            # Wrap around at left edge
            for ii in range(i - width // 2, i + width // 2 + 1):
                for jj in range(360 - width // 2 + j, 360):
                    temp.append(data[:, ii, jj])
                    if np.isnan(np.nanmean(data[:, ii, jj])):
                        n_continents += 1
                for jj in range(0, width - width // 2 + j):
                    temp.append(data[:, ii, jj])
                    if np.isnan(np.nanmean(data[:, ii, jj])):
                        n_continents += 1
        elif j >= 360 - width // 2:
            # Wrap around at right edge
            for ii in range(i - width // 2, i + width // 2 + 1):
                for jj in range(j - width // 2, 360):
                    temp.append(data[:, ii, jj])
                    if np.isnan(np.nanmean(data[:, ii, jj])):
                        n_continents += 1
                for jj in range(0, width // 2 - (359 - j)):
                    temp.append(data[:, ii, jj])
                    if np.isnan(np.nanmean(data[:, ii, jj])):
                        n_continents += 1
        else:
            # Standard case (no longitude wrapping)
            for ii in range(i - width // 2, i + width // 2 + 1):
                for jj in range(j - width // 2, j + width // 2 + 1):
                    temp.append(data[:, ii, jj])
                    if np.isnan(np.nanmean(data[:, ii, jj])):
                        n_continents += 1

        # If too much land in the neighborhood → return NaN
        if n_continents > width**2 // 2:
            return np.empty([n]) * np.nan
        else:
            return np.nanmean(temp, axis=0)

    def compute_bias_from_SeaWiFS(self, data_merged, data_SeaWiFS, lat):
        """
        Compute the bias field by comparing merged dataset with SeaWiFS reference.

        Parameters
        ----------
        data_merged : ndarray (time, lat, lon)
            Merged dataset.
        data_SeaWiFS : ndarray (time, lat, lon)
            SeaWiFS reference dataset.
        lat : ndarray
            Latitude array.

        Returns
        -------
        ndarray (lat, lon)
            Bias field (relative changes between merged and SeaWiFS).
        """
        # Restrict to SeaWiFS period (1998–2010)
        data_period = data_merged[0 : 13 * 12, :, :]
        mask_commun = np.abs(np.sign(data_period)) * np.abs(np.sign(data_SeaWiFS))

        # Mask both datasets consistently
        data_masked = data_period * mask_commun
        data_SeaWiFS_masked = data_SeaWiFS * mask_commun

        # Two reference periods: early vs later SeaWiFS
        reference_periods = np.array([[0, 4 * 12], [5 * 12, 13 * 12]])
        change_medianes_data = np.empty(np.shape(data_masked[0, :, :]))
        change_medianes_SeaWiFS = np.empty(np.shape(data_masked[0, :, :]))

        # Compute relative median change for each grid cell
        for i in tqdm(range(np.shape(data_masked)[1])):
            for j in range(np.shape(data_masked)[2]):
                # For merged dataset
                local_mean = self.get_local_mean(data_masked, i, j, lat)
                change_medianes_data[i, j] = (
                    np.nanmedian(local_mean[reference_periods[1, 0]:reference_periods[1, 1]])
                    - np.nanmedian(local_mean[reference_periods[0, 0]:reference_periods[0, 1]])
                ) / np.nanmedian(local_mean[reference_periods[0, 0]:reference_periods[0, 1]])

                # For SeaWiFS
                local_mean = self.get_local_mean(data_SeaWiFS_masked, i, j, lat)
                change_medianes_SeaWiFS[i, j] = (
                    np.nanmedian(local_mean[reference_periods[1, 0]:reference_periods[1, 1]])
                    - np.nanmedian(local_mean[reference_periods[0, 0]:reference_periods[0, 1]])
                ) / np.nanmedian(local_mean[reference_periods[0, 0]:reference_periods[0, 1]])

        return change_medianes_data - change_medianes_SeaWiFS

    def compute_bias_from_MODIS(self, data_merged, data_MODIS, merged_path, lat):
        """
        Compute the bias field by comparing merged dataset with MODIS reference.

        Parameters
        ----------
        data_merged : ndarray (time, lat, lon)
            Merged dataset.
        data_MODIS : ndarray (time, lat, lon)
            MODIS reference dataset.
        merged_path : Path
            Used to infer product type (affects reference periods).
        lat : ndarray
            Latitude array.

        Returns
        -------
        ndarray (2, lat, lon)
            Bias field (two steps of change detected with MODIS).
        """
        product = Path(merged_path).parent.name

        # MODIS overlap period (2002–2022)
        data_period = data_merged[5 * 12 : 25 * 12, :, :]
        mask_commun = np.abs(np.sign(data_period)) * np.abs(np.sign(data_MODIS))

        data_masked = data_period * mask_commun
        data_MODIS_masked = data_MODIS * mask_commun

        # Reference periods depend on product type
        if product == "Yu2023_v1.0":
            reference_periods = np.array([[0, 8 * 12], [10 * 12, 13 * 12], [18 * 12, 21 * 12]])
        else:
            reference_periods = np.array([[0, 8 * 12], [10 * 12, 13 * 12], [16 * 12, 21 * 12]])

        change_medianes_data = np.empty(np.concatenate([[2], np.shape(data_masked[0, :, :])]))
        change_medianes_MODIS = np.empty(np.concatenate([[2], np.shape(data_masked[0, :, :])]))

        # Compute relative median changes for each grid cell (two steps)
        for i in tqdm(range(np.shape(data_masked)[1])):
            for j in range(np.shape(data_masked)[2]):
                # For merged dataset
                local_mean = self.get_local_mean(data_masked, i, j, lat)
                change_medianes_data[0, i, j] = (
                    np.nanmedian(local_mean[reference_periods[1, 0]:reference_periods[1, 1]])
                    - np.nanmedian(local_mean[reference_periods[0, 0]:reference_periods[0, 1]])
                ) / np.nanmedian(local_mean[reference_periods[0, 0]:reference_periods[0, 1]])
                change_medianes_data[1, i, j] = (
                    np.nanmedian(local_mean[reference_periods[2, 0]:reference_periods[2, 1]])
                    - np.nanmedian(local_mean[reference_periods[1, 0]:reference_periods[1, 1]])
                ) / np.nanmedian(local_mean[reference_periods[1, 0]:reference_periods[1, 1]])

                # For MODIS
                local_mean = self.get_local_mean(data_MODIS_masked, i, j, lat)
                change_medianes_MODIS[0, i, j] = (
                    np.nanmedian(local_mean[reference_periods[1, 0]:reference_periods[1, 1]])
                    - np.nanmedian(local_mean[reference_periods[0, 0]:reference_periods[0, 1]])
                ) / np.nanmedian(local_mean[reference_periods[0, 0]:reference_periods[0, 1]])
                change_medianes_MODIS[1, i, j] = (
                    np.nanmedian(local_mean[reference_periods[2, 0]:reference_periods[2, 1]])
                    - np.nanmedian(local_mean[reference_periods[1, 0]:reference_periods[1, 1]])
                ) / np.nanmedian(local_mean[reference_periods[1, 0]:reference_periods[1, 1]])

        return change_medianes_data - change_medianes_MODIS

    def make_total_bias(self, bias_seawifs, bias_modis, merged_path, length):
        """
        Construct a continuous time series of total bias from SeaWiFS and MODIS contributions.

        Parameters
        ----------
        bias_seawifs : float
            Bias derived from SeaWiFS.
        bias_modis : ndarray (2,)
            Two-step bias derived from MODIS.
        merged_path : Path
            Used to infer product type (affects bias application dates).
        length : int
            Total length of the time series.

        Returns
        -------
        ndarray (length,)
            Time series of bias to subtract from the merged dataset.
        """
        product = Path(merged_path).parent.name

        # Bias application dates depend on product
        if product == "Yu2023_v1.0":
            date_biases = [4 * 12 + 6, 14 * 12, 20 * 12 + 6]
        else:
            date_biases = [4 * 12 + 6, 14 * 12, 19 * 12 + 6]

        # Build total bias stepwise
        total_bias = np.concatenate(
            [
                np.zeros([date_biases[0]]),
                bias_seawifs * np.ones([date_biases[1] - date_biases[0]]),
                (bias_modis[0] + bias_seawifs) * np.ones([date_biases[2] - date_biases[1]]),
                (bias_modis[1] + bias_modis[0] + bias_seawifs) * np.ones([length - date_biases[2]]),
            ]
        )
        return total_bias - np.nanmean(total_bias)  # Center to zero mean

    def correct_dataset(self, merged_path, modis_path, seawifs_path, output_path):
        """
        Apply bias correction to a merged dataset using SeaWiFS and MODIS.

        Parameters
        ----------
        merged_path : str or Path
            Path to merged dataset NetCDF file.
        modis_path : str or Path
            Path to MODIS dataset NetCDF file.
        seawifs_path : str or Path
            Path to SeaWiFS dataset NetCDF file.
        output_path : str or Path
            Path for corrected NetCDF output.
        """
        # Load datasets
        dataset_merged = xr.open_dataset(merged_path)
        data_merged = dataset_merged.CHL.data

        dataset = xr.open_dataset(modis_path)
        data_modis = dataset.chlor_a.isel(lat=slice(None, None, -1), time=slice(5, 21 * 12 - 7)).data

        dataset = xr.open_dataset(seawifs_path)
        data_seawifs = dataset.CHL.isel(lat=slice(None, None, -1)).data

        lat = dataset_merged.lat.data
        lon = dataset_merged.lon.data

        # Compute biases
        biases_SeaWiFS = self.compute_bias_from_SeaWiFS(data_merged, data_seawifs, lat)
        biases_MODIS = self.compute_bias_from_MODIS(data_merged, data_modis, merged_path, lat)

        # Apply correction
        data_debiased = np.empty(np.shape(data_merged)) * np.nan
        mask_ocean_case1 = ~np.isnan(np.nanmean(data_merged, axis=0))  # Mask land

        for i in tqdm(range(len(lat))):
            for j in range(len(lon)):
                if mask_ocean_case1[i, j]:
                    total_bias = self.make_total_bias(
                        biases_SeaWiFS[i, j], biases_MODIS[:, i, j], merged_path, len(data_merged[:, 0, 0])
                    )
                    # Subtract bias scaled by local median
                    data_debiased[:, i, j] = data_merged[:, i, j] - total_bias * np.nanmedian(
                        data_merged[:, i, j]
                    )

        # Save corrected dataset
        nc_file = xr.Dataset(
            data_vars=dict(CHL=(["time", "lat", "lon"], data_debiased)),
            coords=dict(time=dataset_merged.time.data, lat=lat, lon=lon),
        )
        nc_file.to_netcdf(output_path)

