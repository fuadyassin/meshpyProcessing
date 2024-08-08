import os
import argparse
import numpy as np
import xarray as xs
import geopandas as gpd
import glob
from natsort import natsorted

def remap_rdrs_climate_data(input_directory, output_directory, input_basin, input_ddb, start_year, end_year):
    # Ensure output directory exists
    os.makedirs(output_directory, exist_ok=True)

    # Read basin and drainage database files
    basin = gpd.read_file(input_basin)
    db = xs.open_dataset(input_ddb)
    lon = db.variables['lon'].values
    lat = db.variables['lat'].values
    segid = db.variables['subbasin'].values
    db.close()

    # Print some basin and database information
    print("Basin Info:")
    print(basin.head())  # Print the first few rows of the basin GeoDataFrame
    print("Longitude:", lon[:5])  # Print the first few longitude values
    print("Latitude:", lat[:5])  # Print the first few latitude values
    print("Subbasin IDs:", segid[:5])  # Print the first few subbasin IDs

    # List files based on year range
    files = []
    for year in range(start_year, end_year + 1):
        year_files = glob.glob(os.path.join(input_directory, f"*_{str(year)}*.nc"))
        files.extend(year_files)
        print(f"Files for year {year}: {year_files}")

    # Sort files in natural order
    files = natsorted(files)
    print("Sorted files:", files)

    # Process each file
    for file_path in files:
        print(f"Processing file: {file_path}")
        process_file(file_path, segid, lon, lat, output_directory)

def remap_rdrs_climate_data_single_year(input_directory, output_directory, input_basin, input_ddb, year):
    remap_rdrs_climate_data(input_directory, output_directory, input_basin, input_ddb, year, year)

def process_file(file_path, segid, lon, lat, output_directory):
    # Placeholder for the actual file processing logic
    print(f"Started processing file: {file_path}")
    # Open dataset
    forc = xs.open_dataset(file_path)

    # Extract indices of forcing IDs based on the drainage database
    ind = []
    for i in range(len(segid)):
        fid = np.where(np.int32(forc['COMID'].values) == segid[i])[0]
        ind = np.append(ind, fid)
    ind = np.int32(ind)

    # Create a new dataset with data ordered as needed
    forc_vec = xs.Dataset()
    variables_to_process = ['RDRS_v2.1_A_PR0_SFC', 'RDRS_v2.1_P_P0_SFC', 'RDRS_v2.1_P_HU_09944', 
                            'RDRS_v2.1_P_TT_09944', 'RDRS_v2.1_P_FB_SFC', 'RDRS_v2.1_P_FI_SFC', 'RDRS_v2.1_P_UVC_09944']
    for var in variables_to_process:
        data = forc[var].values[:, ind]
        forc_vec[var] = (("time", "subbasin"), data)
        forc_vec[var].attrs = forc[var].attrs

    # Correctly setting coordinates:
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

    # Metadata and attributes
    forc_vec.attrs.update({
        'Conventions': 'CF-1.6',
        'history': 'Processed on Apr 06, 2024',
        'License': 'The data were written by Fuad Yassin.',
        'featureType': 'timeSeries'
    })

    # Define coordinate system
    forc_vec['crs'] = xs.DataArray(np.int32(1), attrs={
        'grid_mapping_name': 'latitude_longitude',
        'longitude_of_prime_meridian': 0.0,
        'semi_major_axis': 6378137.0,
        'inverse_flattening': 298.257223563
    })

    # Define a variable for the points and set the 'timeseries_id'
    forc_vec['subbasin'] = xs.DataArray(segid, dims=['subbasin'])
    forc_vec['subbasin'].attrs.update({
        'long_name': 'shape_id',
        'units': '1',
        'cf_role': 'timeseries_id'
    })

    # Save to netCDF
    output_path = os.path.join(output_directory, os.path.basename(file_path).replace('.nc', '_modified.nc'))
    encoding = {var: {'zlib': True, 'complevel': 6} for var in forc_vec.data_vars}
    forc_vec.to_netcdf(output_path, encoding=encoding)
    print(f"Processed and saved: {output_path}")

    # Close the dataset
    forc.close()
    print(f"Finished processing file: {file_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process RDRS climate data.")
    subparsers = parser.add_subparsers(dest="command")

    all_years_parser = subparsers.add_parser("all_years")
    all_years_parser.add_argument("--input_directory", required=True)
    all_years_parser.add_argument("--output_directory", required=True)
    all_years_parser.add_argument("--input_basin", required=True)
    all_years_parser.add_argument("--input_ddb", required=True)
    all_years_parser.add_argument("--start_year", type=int, required=True)
    all_years_parser.add_argument("--end_year", type=int, required=True)

    single_year_parser = subparsers.add_parser("single_year")
    single_year_parser.add_argument("--input_directory", required=True)
    single_year_parser.add_argument("--output_directory", required=True)
    single_year_parser.add_argument("--input_basin", required=True)
    single_year_parser.add_argument("--input_ddb", required=True)
    single_year_parser.add_argument("--year", type=int, required=True)

    args = parser.parse_args()

    if args.command == "all_years":
        remap_rdrs_climate_data(
            args.input_directory,
            args.output_directory,
            args.input_basin,
            args.input_ddb,
            args.start_year,
            args.end_year
        )
    elif args.command == "single_year":
        remap_rdrs_climate_data_single_year(
            args.input_directory,
            args.output_directory,
            args.input_basin,
            args.input_ddb,
            args.year
        )
