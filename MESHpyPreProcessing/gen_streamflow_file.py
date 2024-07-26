# gen_streamflow_file.py

import pandas as pd
import numpy as np
import requests
from owslib.ogcapi.features import Features
from datetime import datetime
import time

class GenStreamflowFile:
    def __init__(self):
        self.oafeat = Features("https://api.weather.gc.ca/")

    def create_date_range(self, start_date, end_date):
        return pd.date_range(start=start_date, end=end_date)

    def extract_flow_data_us(self, station_list, start_date, end_date):
        dates = self.create_date_range(start_date, end_date)
        data_dict = {'Date': dates}
        date_index_dict = {str(date.date()): idx for idx, date in enumerate(dates)}
        station_info = []

        for station in station_list:
            data_dict[station] = [np.nan] * len(dates)
        
        for station in station_list:
            start_time_station = time.time()
            url = f"https://waterservices.usgs.gov/nwis/dv/?format=json&sites={station}&startDT={start_date}&endDT={end_date}&parameterCd=00060&statCd=00003"
            response = requests.get(url)
            if response.status_code == 200:
                data = response.json()
                if 'timeSeries' in data['value']:
                    time_series = data['value']['timeSeries'][0]
                    variable_info = time_series['variable']
                    unit = variable_info.get('unit', {}).get('unitCode', None)
                    parameter_units = variable_info.get('variableDescription', None)

                    records = time_series['values'][0]['value']
                    flow_data = pd.DataFrame(records)
                    flow_data['value'] = pd.to_numeric(flow_data['value'], errors='coerce')

                    # Convert flow data to cms if the unit is cfs
                    if unit == 'ft3/s' or parameter_units == 'Cubic Feet per Second':
                        flow_data['value'] = flow_data['value'] * 0.0283168
                    
                    flow_data['dateTime'] = pd.to_datetime(flow_data['dateTime']).dt.date.astype(str)

                    start_time_loop = time.time()
                    for date, flow in zip(flow_data['dateTime'], flow_data['value']):
                        if date in date_index_dict:
                            date_index = date_index_dict[date]
                            data_dict[station][date_index] = flow
                    end_time_loop = time.time()
                    print(f"Time taken in the for loop for station {station}: {end_time_loop - start_time_loop} seconds")
                    
                    site_info = time_series['sourceInfo']
                    station_info.append({
                        'Station_Number': site_info['siteCode'][0]['value'],
                        'Station_Name': site_info['siteName'],
                        'Latitude': site_info['geoLocation']['geogLocation']['latitude'],
                        'Longitude': site_info['geoLocation']['geogLocation']['longitude'],
                        'Drainage_Area': next((prop['value'] for prop in site_info['siteProperty'] if prop['name'] == 'drain_area_va'), None),
                        'Unit': unit,
                        'Parameter_Units': parameter_units
                    })
                else:
                    print(f"Flow data column not found for station: {station}")
            else:
                print(f"Failed to retrieve data for station: {station}")
            
            end_time_station = time.time()
            print(f"Time taken to retrieve data for station {station}: {end_time_station - start_time_station} seconds")

        combined_df = pd.DataFrame(data_dict)
        return combined_df, station_info

    def fetch_hydrometric_data_ca(self, station_numbers, start_date, end_date, limit=1000):
        dates = self.create_date_range(start_date, end_date)
        data_dict = {'Date': dates}
        date_index_dict = {str(date.date()): idx for idx, date in enumerate(dates)}
        station_info = []

        for station_number in station_numbers:
            data_dict[station_number] = [np.nan] * len(dates)
            
            offset = 0
            full_data = []
            
            start_time_station = time.time()
            
            while True:
                url = f"https://api.weather.gc.ca/collections/hydrometric-daily-mean/items"
                params = {
                    'STATION_NUMBER': station_number,
                    'datetime': f"{start_date}/{end_date}",
                    'limit': limit,
                    'offset': offset,
                    'f': 'json'
                }

                response = requests.get(url, params=params)
                response_data = response.json()

                if 'features' in response_data and response_data['features']:
                    full_data.extend(response_data['features'])
                    offset += limit
                    if len(response_data['features']) < limit:
                        break
                else:
                    break
            
            if full_data:
                data_list = [
                    {
                        'Date': feature['properties']['DATE'],
                        'value': feature['properties']['DISCHARGE'] if feature['properties']['DISCHARGE'] is not None else -1000
                    }
                    for feature in full_data
                ]
                
                flow_data = pd.DataFrame(data_list)
                flow_data['value'] = pd.to_numeric(flow_data['value'], errors='coerce')
                flow_data['Date'] = pd.to_datetime(flow_data['Date']).dt.date.astype(str)

                start_time_loop = time.time()
                for date, flow in zip(flow_data['Date'], flow_data['value']):
                    if date in date_index_dict:
                        date_index = date_index_dict[date]
                        data_dict[station_number][date_index] = flow
                end_time_loop = time.time()
                print(f"Time taken in the for loop for station {station_number}: {end_time_loop - start_time_loop} seconds")

                first_feature = full_data[0]['properties']
                geometry = full_data[0]['geometry']
                station_info.append({
                    'Station_Number': first_feature['STATION_NUMBER'],
                    'Station_Name': first_feature['STATION_NAME'],
                    'Latitude': geometry['coordinates'][1],
                    'Longitude': geometry['coordinates'][0],
                    'Drainage_Area': first_feature.get('DRAINAGE_AREA_GROSS', None)
                })
            else:
                print(f"Flow data not found for station: {station_number}")
            
            end_time_station = time.time()
            print(f"Time taken to retrieve data for station {station_number}: {end_time_station - start_time_station} seconds")

        combined_df = pd.DataFrame(data_dict)
        return combined_df, station_info

    def write_flow_data_to_file_obstxt(self, file_path, flow_data, site_details):
        flow_data = flow_data.fillna(-1.000)
        flow_data['Date'] = pd.to_datetime(flow_data['Date'])

        with open(file_path, "w") as file_conn:
            start_date = flow_data['Date'].min()
            end_date = flow_data['Date'].max()
            file_conn.write(f"Observedstreamflow\t{start_date.strftime('%Y/%m/%d')}\t{end_date.strftime('%Y/%m/%d')}\n")
            num_stations = flow_data.shape[1] - 1
            num_days = flow_data.shape[0]
            start_year = start_date.strftime('%Y')
            start_day_of_year = start_date.timetuple().tm_yday
            file_conn.write(f"{num_stations}  {num_days}  {num_days}  24 {start_year}  {start_day_of_year} 00\n")
            
            for station_id in flow_data.columns[1:]:
                station_info = next((item for item in site_details if item["Station_Number"] == station_id), None)
                if station_info:
                    lat = station_info['Latitude']
                    lon = station_info['Longitude']
                    drainage_area = station_info['Drainage_Area']
                    if drainage_area is None:
                        drainage_area = -1.0
                    station_name = station_info['Station_Name']
                    file_conn.write(f"{int(lat * 60):4d} {int(lon * 60):4d} {lat:12.6f} {lon:12.6f} {station_id:12s} {float(drainage_area):12.3f} {station_name}\n")
            
            for i in range(num_days):
                flow_values = flow_data.iloc[i, 1:].values
                formatted_flow_values = " ".join(f"{x:12.4f}" for x in flow_values)
                file_conn.write(f"{formatted_flow_values}\n")

    def write_flow_data_to_file_ensim(self, file_path, flow_data, site_details):
        flow_data = flow_data.fillna(-1.000)

        header = [
            "########################################",
            ":FileType tb0  ASCII  EnSim 1.0",
            "#",
            "# DataType               Time Series",
            "#",
            ":Application             EnSimHydrologic",
            ":Version                 2.1.23",
            ":WrittenBy          PythonScript",
            f":CreationDate       {datetime.now().strftime('%Y-%m-%d')}",
            "#",
            "#---------------------------------------",
            ":SourceFile                   flow_data",
            "#",
            ":Name               streamflow",
            "#",
            ":Projection         LATLONG",
            ":Ellipsoid          WGS84",
            "#",
            f":StartTime          {flow_data['Date'].iloc[0]} 00:00:00.00000",
            "#",
            ":AttributeUnits            1.0000000",
            ":DeltaT               24",
            ":RoutingDeltaT         1",
            "#",
            ":ColumnMetaData",
            f"   :ColumnUnits             {' '.join(['m3/s' for _ in range(flow_data.shape[1] - 1)])}",
            f"   :ColumnType             {' '.join(['float' for _ in range(flow_data.shape[1] - 1)])}",
            f"   :ColumnName           {' '.join(flow_data.columns[1:])}",
            "   :ColumnLocationX    " + ' '.join([f"{site['Longitude'] * 60:.6f}" for site in site_details]),
            "   :ColumnLocationY    " + ' '.join([f"{site['Latitude'] * 60:.6f}" for site in site_details]),
            f"   :coeff1            {' '.join(['0.0000E+00' for _ in range(flow_data.shape[1] - 1)])}",
            f"   :coeff2            {' '.join(['0.0000E+00' for _ in range(flow_data.shape[1] - 1)])}",
            f"   :coeff3            {' '.join(['0.0000E+00' for _ in range(flow_data.shape[1] - 1)])}",
            f"   :coeff4            {' '.join(['0.0000E+00' for _ in range(flow_data.shape[1] - 1)])}",
            f"   :Value1            {' '.join(['1' for _ in range(flow_data.shape[1] - 1)])}",
            ":EndColumnMetaData",
            ":endHeader"
        ]
        
        with open(file_path, "w") as file_conn:
            file_conn.write("\n".join(header) + "\n")
            
            for _, row in flow_data.iterrows():
                flows = row[1:].values
                flow_string = " ".join(f"{flow:12.4f}" for flow in flows)
                file_conn.write(f"{' ' * 22}{flow_string}\n")
