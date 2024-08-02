import os
import numpy as np
import pandas as pd
import geopandas as gpd
import netCDF4 as nc
import xarray as xs
from functools import reduce
from datetime import datetime
import tempfile

class GSDESoil:
    def __init__(self, directory, input_basin, output_shapefile, nc_filename):
        self.directory = directory
        self.input_basin = input_basin
        self.output_shapefile = output_shapefile
        self.nc_filename = nc_filename
        self.file_paths = []
        self.gsde_df = pd.DataFrame()
        self.merged_gdf = gpd.GeoDataFrame()
        self.weights_used = []
        self.mesh_intervals = []
        self.lon = []
        self.lat = []
        self.segid = []

    def load_data(self, file_names):
        """
        Load data from multiple CSV files and merge them into a single DataFrame.
        """
        self.file_paths = [os.path.join(self.directory, filename) for filename in file_names]
        self.gsde_df = self.load_and_merge_files(self.file_paths)

    @staticmethod
    def load_and_merge_files(file_list, key='COMID'):
        """
        Load and merge multiple CSV files on the specified key.
        """
        dfs = [pd.read_csv(fp) for fp in file_list]
        return reduce(lambda left, right: pd.merge(left, right, on=key, how='outer'), dfs)

    def clean_column_names(self, search_list, replace_list):
        """
        Clean column names by replacing search strings with corresponding replace strings.
        """
        new_columns = []
        for column in self.gsde_df.columns:
            for search, replace in zip(search_list, replace_list):
                if search in column:
                    column = column.replace(search, replace)
            column = column.replace('.', '')  # Remove periods from column names
            new_columns.append(column)
        self.gsde_df.columns = new_columns

    def fill_and_clean_data(self, exclude_cols=['COMID'], exclude_patterns=['OC', 'BD'], max_val=100):
        """
        Fill and clean data, replacing values greater than max_val with NaN and forward/backward filling.
        """
        for col in self.gsde_df.columns:
            if col not in exclude_cols and not any(pattern in col for pattern in exclude_patterns):
                self.gsde_df.loc[self.gsde_df[col] > max_val, col] = np.nan
        self.gsde_df.sort_values('COMID', inplace=True)
        self.gsde_df.fillna(method='ffill', inplace=True)
        self.gsde_df.fillna(method='bfill', inplace=True)

    def calculate_weights(self, gsde_intervals, mesh_intervals):
        """
        Calculate weights for different soil intervals.
        """
        self.mesh_intervals = mesh_intervals
        weights_used = []
        for mesh_interval in mesh_intervals:
            start, end = mesh_interval
            weights = []
            for gsde_interval in gsde_intervals:
                gsde_start, gsde_end = gsde_interval
                overlap_start = max(start, gsde_start)
                overlap_end = min(end, gsde_end)
                weight = (overlap_end - overlap_start) / (end - start) if overlap_start < overlap_end else 0
                weights.append(weight)
            weights_used.append([w / sum(weights) for w in weights if sum(weights) > 0])
        self.weights_used = weights_used

    def calculate_mesh_values(self, column_names):
        """
        Calculate mesh values by applying weights to soil property columns.
        """
        for prop, cols in column_names.items():
            extracted_data = self.gsde_df[cols]
            weights_array = np.array(self.weights_used).T
            mesh_values = np.dot(extracted_data, weights_array)
            if prop == 'OC':
                mesh_values *= 0.01 * 1.72
            for i, mesh_col in enumerate([f'mesh{prop}{j+1}' for j in range(len(self.mesh_intervals))]):
                self.gsde_df[mesh_col] = mesh_values[:, i]

    def merge_and_save_shapefile(self):
        """
        Merge soil data with a shapefile and save the result.
        """
        gdf = gpd.read_file(self.input_basin).to_crs(epsg=4326)
        gdf['COMID'] = gdf['COMID'].astype(int)
        self.merged_gdf = gdf.merge(self.gsde_df, on='COMID', how='left')
        self.merged_gdf.to_file(self.output_shapefile)

    def set_coordinates(self, input_ddb):
        """
        Set longitude and latitude values from a NetCDF drainage database.
        """
        db = xs.open_dataset(input_ddb)
        self.lon = db.variables['lon'].values
        self.lat = db.variables['lat'].values
        self.segid = db.variables['subbasin'].values
        db.close()

    def create_netcdf_file(self, ncname, properties, lon, lat, ind):
        """
        Create a NetCDF file with the processed soil data.
        """
        try:
            rootgrp = nc.Dataset(self.nc_filename, "w", format="NETCDF4")
        except PermissionError:
            temp_dir = tempfile.gettempdir()
            self.nc_filename = os.path.join(temp_dir, f"MESH_parameters_{ncname}.nc")
            rootgrp = nc.Dataset(self.nc_filename, "w", format="NETCDF4")

        # Calculate the indices of the matching COMID values
        ind = []
        for i in range(len(self.segid)):
            fid = np.where(np.int32(self.merged_gdf['COMID'].values) == self.segid[i])[0]
            ind = np.append(ind, fid)
        ind = np.int32(ind)
        
        MaxNumGRUs = 16
        num_soil_lyrs = 4
        subbasin_dim = rootgrp.createDimension("subbasin", len(lon))
        ngru_dim = rootgrp.createDimension("ngru", MaxNumGRUs)
        nsol_dim = rootgrp.createDimension("nsol", num_soil_lyrs)

        lon_var = rootgrp.createVariable("lon", "f4", ("subbasin",), fill_value=-1.0)
        lat_var = rootgrp.createVariable("lat", "f4", ("subbasin",), fill_value=-1.0)
        time_var = rootgrp.createVariable("time", "f4", ("subbasin",), fill_value=-1.0)
        lon_var.units = "degrees_east"
        lat_var.units = "degrees_north"
        time_var.units = "days since 1980-10-01 00:00:00.0 -0:00"

        lon_var[:] = np.array(lon[ind])
        lat_var[:] = np.array(lat[ind])
        time_var[:] = np.zeros(len(lon))

        for prop in properties:
            for i, depth in enumerate(['1', '2', '3', '4']):
                data_var = rootgrp.createVariable(f'{prop.lower()}{depth}', "f4", ("nsol", "subbasin"), fill_value=-1.0)
                data_var.long_name = f"{prop} Content of Soil Layer {depth}"
                data_var[:] = np.array(self.merged_gdf[f'mesh{prop}{depth}'].values[ind])

        rootgrp.Conventions = "CF-1.0"
        rootgrp.source = "MERIT geogabrics and GSDE soil"
        rootgrp.institution = "ECCC"
        rootgrp.references = "xx et al. (xxxx) journal xx:xx-xx"
        rootgrp.history = f"Fuad Yassin, {datetime.now().strftime('%Y-%m-%d')}"
        rootgrp.featureType = "point"
        rootgrp.close()

