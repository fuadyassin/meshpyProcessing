import os
import numpy as np
import geopandas as gpd
import netCDF4 as nc
from datetime import datetime
import tempfile
import xarray as xs

class NetCDFWriter:
    def __init__(self, nc_filename, shapefile_path, input_ddb_path):
        self.nc_filename = nc_filename
        self.shapefile_path = shapefile_path
        self.input_ddb_path = input_ddb_path
        self.merged_gdf = gpd.GeoDataFrame()
        self.lon = []
        self.lat = []
        self.segid = []
        self.num_soil_lyrs = 0

    def read_shapefile(self):
        """
        Read the shapefile and set the merged GeoDataFrame.
        """
        self.merged_gdf = gpd.read_file(self.shapefile_path).to_crs(epsg=4326)

    def set_coordinates(self):
        """
        Set longitude and latitude values from a NetCDF drainage database.
        """
        db = xs.open_dataset(self.input_ddb_path)
        self.lon = db.variables['lon'].values
        self.lat = db.variables['lat'].values
        self.segid = db.variables['subbasin'].values
        db.close()

    def set_num_soil_layers(self, num_layers):
        """
        Set the number of soil layers.
        """
        self.num_soil_lyrs = num_layers

    def write_netcdf(self, properties):
        """
        Create a NetCDF file with the processed soil data.
        properties: dict
            A dictionary with two keys:
            - 'layer_dependent': List of property names tied to the number of soil layers.
            - 'layer_independent': List of property names dependent only on the subbasin.
        """
        try:
            rootgrp = nc.Dataset(self.nc_filename, "w", format="NETCDF4")
        except PermissionError:
            temp_dir = tempfile.gettempdir()
            self.nc_filename = os.path.join(temp_dir, f"MESH_parameters.nc")
            rootgrp = nc.Dataset(self.nc_filename, "w", format="NETCDF4")

        # Calculate the indices of the matching COMID values
        ind = []
        for i in range(len(self.segid)):
            fid = np.where(np.int32(self.merged_gdf['COMID'].values) == self.segid[i])[0]
            ind = np.append(ind, fid)
        ind = np.int32(ind)

        subbasin_dim = rootgrp.createDimension("subbasin", len(self.lon))
        nsol_dim = rootgrp.createDimension("nsol", self.num_soil_lyrs)

        lon_var = rootgrp.createVariable("lon", "f4", ("subbasin",), fill_value=-1.0)
        lat_var = rootgrp.createVariable("lat", "f4", ("subbasin",), fill_value=-1.0)
        time_var = rootgrp.createVariable("time", "f4", ("subbasin",), fill_value=-1.0)
        lon_var.units = "degrees_east"
        lat_var.units = "degrees_north"
        time_var.units = "days since 1980-10-01 00:00:00.0 -0:00"

        lon_var[:] = np.array(self.lon[ind])
        lat_var[:] = np.array(self.lat[ind])
        time_var[:] = np.zeros(len(self.lon))

        # Handle properties tied to number of soil layers
        if 'layer_dependent' in properties:
            for prop in properties['layer_dependent']:
                data_var = rootgrp.createVariable(f'{prop.lower()}', "f4", ("subbasin", "nsol"), fill_value=-1.0)
                data_var.long_name = f"{prop} Content of Soil Layer"
                for i in range(self.num_soil_lyrs):
                    data_var[:, i] = np.array(self.merged_gdf[f'mesh{prop}{i+1}'].values[ind])

        # Handle properties dependent only on subbasin
        if 'layer_independent' in properties:
            for prop in properties['layer_independent']:
                data_var = rootgrp.createVariable(f'{prop.lower()}', "f4", ("subbasin",), fill_value=-1.0)
                data_var.long_name = f"{prop} Content per Subbasin"
                data_var[:] = np.array(self.merged_gdf[prop].values[ind])

        rootgrp.Conventions = "CF-1.0"
        rootgrp.source = "MERIT geogabrics and GSDE soil"
        rootgrp.institution = "ECCC"
        rootgrp.references = "xx et al. (xxxx) journal xx:xx-xx"
        rootgrp.history = f"Fuad Yassin, {datetime.now().strftime('%Y-%m-%d')}"
        rootgrp.featureType = "point"
        rootgrp.close()
