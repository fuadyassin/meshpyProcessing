# MESHpyPreProcessing/remap_rdrs_climate_data.py

import os
import numpy as np
import xarray as xs
import geopandas as gpd
import glob

def remap_rdrs_climate_data(input_directory, output_directory, input_basin, input_ddb, start_year, end_year):
    # Existing implementation for multiple years
    os.makedirs(output_directory, exist_ok=True)
    basin = gpd.read_file(input_basin)
    db = xs.open_dataset(input_ddb)
    lon = db.variables['lon'].values
    lat = db.variables['lat'].values
    segid = db.variables['subbasin'].values
    db.close()

    files = []
    for year in range(start_year, end_year + 1):
        files.extend(glob.glob(os.path.join(input_directory, f"remapped_remapped_ncrb_model_{str(year)}*.nc")))

    for file_path in files:
        process_file(file_path, segid, lon, lat, output_directory)

def remap_rdrs_climate_data_single_year(input_directory, output_directory, input_basin, input_ddb, year):
    # Implementation for a single year
    os.makedirs(output_directory, exist_ok=True)
    basin = gpd.read_file(input_basin)
    db = xs.open_dataset(input_ddb)
    lon = db.variables['lon'].values
    lat = db.variables['lat'].values
    segid = db.variables['subbasin'].values
    db.close()

    files = glob.glob(os.path.join(input_directory, f"remapped_remapped_ncrb_model_{str(year)}*.nc"))
    
    for file_path in files:
        process_file(file_path, segid, lon, lat, output_directory)

def process_file(file_path, segid, lon, lat, output_directory):
    # Helper function to process a single file
    forc = xs.open_dataset(file_path)
    ind = []
    for i in range(len(segid)):
        fid = np.where(np.int32(forc['COMID'].values) == segid[i])[0]
        ind = np.append(ind, fid)
    ind = np.int32(ind)

    forc_vec = xs.Dataset()
    variables_to_process = [
        'RDRS_v2.1_A_PR0_SFC', 'RDRS_v2.1_P_P0_SFC', 'RDRS_v2.1_P_HU_09944',
        'RDRS_v2.1_P_TT_09944', 'RDRS_v2.1_P_FB_SFC', 'RDRS_v2.1_P_FI_SFC', 'RDRS_v2.1_P_UVC_09944'
    ]
    for var in variables_to_process:
        data = forc[var].values[:, ind]
        forc_vec[var] = (("time", "subbasin"), data)
        forc_vec[var].attrs = forc[var].attrs

    forc_vec = forc_vec.assign_coords(
        time=forc['time'].values,
        lon=(['subbasin'], lon[ind]),
        lat=(['subbasin'], lat[ind])
    )
    forc_vec['lon'].attrs = {
        'long_name': 'longitude',
        'units': 'degrees_east'
    }
    forc_vec['lat'].attrs = {
        'long_name': 'latitude',
        'units': 'degrees_north'
    }

    forc_vec.attrs.update({
        'Conventions': 'CF-1.6',
        'history': 'Processed on Aug 06, 2024',
        'License': 'The data were written by Fuad Yassin.',
        'featureType': 'timeSeries'
    })

    forc_vec['crs'] = xs.DataArray(np.int32(1), attrs={
        'grid_mapping_name': 'latitude_longitude',
        'longitude_of_prime_meridian': 0.0,
        'semi_major_axis': 6378137.0,
        'inverse_flattening': 298.257223563
    })

    forc_vec['subbasin'] = xs.DataArray(segid, dims=['subbasin'])
    forc_vec['subbasin'].attrs.update({
        'long_name': 'shape_id',
        'units': '1',
        'cf_role': 'timeseries_id'
    })

    output_path = os.path.join(output_directory, os.path.basename(file_path).replace('.nc', '_modified.nc'))
    encoding = {var: {'zlib': True, 'complevel': 6} for var in forc_vec.data_vars}
    forc_vec.to_netcdf(output_path, encoding=encoding)
    print(f"Processed and saved: {output_path}")

    forc.close()
