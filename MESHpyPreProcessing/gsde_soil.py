import os
import numpy as np
import pandas as pd
import geopandas as gpd
from functools import reduce
import xarray as xs

class GSDESoil:
    def __init__(self, directory, input_basin, output_shapefile):
        self.directory = directory
        self.input_basin = input_basin
        self.output_shapefile = output_shapefile
        self.file_paths = []
        self.gsde_df = pd.DataFrame()
        self.merged_gdf = gpd.GeoDataFrame()
        self.weights_used = []
        self.mesh_intervals = []
        self.lon = []
        self.lat = []
        self.segid = []
        self.num_soil_lyrs = 0

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

    def fill_and_clean_data(self, exclude_cols=['COMID'], exclude_patterns=['OC', 'BD', 'BDRICM', 'BDTICM'], max_val=100):
        """
        Fill and clean data, replacing values greater than max_val with NaN and forward/backward filling.
        """
        for col in self.gsde_df.columns:
            if col not in exclude_cols:
                if not any(pattern in col for pattern in exclude_patterns):
                    # Replace values greater than max_val with NaN
                    self.gsde_df.loc[self.gsde_df[col] > max_val, col] = np.nan
            
            # Setting the possible minimum and maximum values for the depth to bedrock.
            if 'BDRICM' in col:
                # Replace values of 9999.0 with NaN
                self.gsde_df.loc[self.gsde_df[col] == 9999.0, col] = np.nan
                # Convert from cm into meters by dividing the column values by 100.0
                self.gsde_df[col] = self.gsde_df[col] / 100.0
                # Cap values greater than 2.0 to 2.0
                self.gsde_df.loc[self.gsde_df[col] > 2.0, col] = 2.0
                self.gsde_df.loc[self.gsde_df[col] < 0.1, col] = 0.1
                
            if 'BDTICM' in col:
                # Replace values of 9999.0 with NaN
                self.gsde_df.loc[self.gsde_df[col] == 9999.0, col] = np.nan
                # Convert from cm into meters by dividing the column values by 100.0
                self.gsde_df[col] = self.gsde_df[col] / 100.0
                # Cap values greater than 4.1 to 4.1
                self.gsde_df.loc[self.gsde_df[col] > 4.1, col] = 4.1
                self.gsde_df.loc[self.gsde_df[col] < 0.1, col] = 0.1
                
        self.gsde_df.sort_values('COMID', inplace=True)
        self.gsde_df.fillna(method='ffill', inplace=True)
        self.gsde_df.fillna(method='bfill', inplace=True)

    def calculate_weights(self, gsde_intervals, mesh_intervals):
        """
        Calculate weights for different soil intervals.
        """
        self.mesh_intervals = mesh_intervals
        self.num_soil_lyrs = len(mesh_intervals)  # Set the number of soil layers based on mesh intervals
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
