# meshpyProcessing



## Installation

You can install the package using pip:

```bash
pip install git+https://github.com/fuadyassin/meshpyProcessing.git


from MESHpyPreProcessing.gen_streamflow_file import GenStreamflowFile

gen_flow = GenStreamflowFile()
combined_df, station_info = gen_flow.extract_flow_data_us(station_list, start_date, end_date)
gen_flow.write_flow_data_to_file_obstxt('output.txt', combined_df, station_info)
