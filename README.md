# meshpyProcessing



## Installation
You can install the package using pip:

```bash
pip install git+https://github.com/fuadyassin/meshpyProcessing.git
```
# Streamflow file preparation

from MESHpyPreProcessing.gen_streamflow_file import GenStreamflowFile

gen_flow = GenStreamflowFile()

## Define sample lists of station IDs
station_ca = ["05GG001", "05AC012"]

station_us = ["06132200", "05020500"]

## Define date range

start_date = "1980-03-01"

end_date = "2018-01-10"

combined_data_ca, station_info_ca = gen_flow.fetch_hydrometric_data_ca(station_ca, start_date, end_date)

combined_data_us, station_info_us = gen_flow.extract_flow_data_us(station_us, start_date, end_date)

## Combine data into a single DataFrame
combined_data = pd.merge(combined_data_ca, combined_data_us, on='Date', how='outer')

## Combine station info
combined_station_info = station_info_ca + station_info_us

## Write the data to a file
gen_flow.write_flow_data_to_file_obstxt('output.txt', combined_data, combined_station_info)

gen_flow.write_flow_data_to_file_ensim('output_ensim.txt', combined_data, combined_station_info)
