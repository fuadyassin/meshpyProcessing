import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
import netCDF4 as nc
import matplotlib.colors as mcolors
import math

def plot_grus_ddb(output_basin_path, ddbnetcdf_path, save_path='landuse_gru_plots_with_basin_boundary.png'):
    """
    Plots GRU data from a NetCDF file overlaid with the dissolved boundary of a basin shapefile.

    Parameters:
    - output_basin_path (str): Path to the shapefile containing the basin data.
    - ddbnetcdf_path (str): Path to the NetCDF file containing the GRU data.
    - save_path (str): File path where the plot will be saved. Defaults to 'landuse_gru_plots_with_basin_boundary.png'.

    Returns:
    - None: The function saves the plot to the specified path.
    """
    
    # Load and dissolve the shapefile to get a single boundary
    sub_agg_gdf = gpd.read_file(output_basin_path).to_crs(epsg=4326)
    dissolved_basin = sub_agg_gdf.dissolve().boundary

    # Load NetCDF data
    with nc.Dataset(ddbnetcdf_path) as dataset:
        # Extract variables directly into a DataFrame
        drainage_df = pd.DataFrame({
            'subbasin': dataset.variables['subbasin'][:].astype(int)
        })

        # Extract land use classes and append as separate columns
        landuse_classes = dataset.variables['LandUse'][:]
        gru_class = dataset.variables['GRU'][:]
        for i, landuse in enumerate(landuse_classes):
            drainage_df[landuse] = gru_class[:, i]

    # Merge subbasin agg with drainage database
    sub_agg_ddb_merged_gdf = sub_agg_gdf.merge(drainage_df, left_on='COMID', right_on='subbasin', how='left')

    # Identify the columns corresponding to landuse classes dynamically
    gru_columns = [col for col in sub_agg_ddb_merged_gdf.columns if col in landuse_classes]

    # Pre-calculate the percentage values and store them along with column names
    percentage_pairs = [(col, round(sub_agg_ddb_merged_gdf[col].sum() / len(sub_agg_ddb_merged_gdf) * 100, 3)) for col in gru_columns]

    # Sort the columns by their percentage values in descending order
    percentage_pairs_sorted = sorted(percentage_pairs, key=lambda x: x[1], reverse=True)

    # Unpack the sorted columns and percentages
    sorted_gru_columns, sorted_rounded_percentages = zip(*percentage_pairs_sorted)

    # Determine the layout of the subplots
    num_plots = len(sorted_gru_columns)
    num_rows = math.floor(math.sqrt(num_plots))
    num_cols = math.ceil(num_plots / num_rows)
    
    # Set the figure size based on the number of rows and columns
    fig_width = 2 * num_cols
    fig_height = 2 * num_rows

    # Create a copy of the colormap and modify it
    cmap = plt.cm.get_cmap("gnuplot2_r").copy()
    cmap.set_under('white', alpha=0)

    # Function to split title into two lines
    def split_title(title):
        words = title.split()
        split_point = len(words) // 2
        return ' '.join(words[:split_point]) + '\n' + ' '.join(words[split_point:])

    # Set up the figure and axes for the calculated grid
    fig, axes = plt.subplots(nrows=num_rows, ncols=num_cols, figsize=(fig_width, fig_height), sharex=True, sharey=True)
    fig.subplots_adjust(hspace=0.6, wspace=0.1)

    # Create an Image for colormap scaling
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=mcolors.Normalize(vmin=0.01, vmax=1))
    sm._A = []

    # Loop through and plot each column in the new order
    for i, (col, rounded_percentage) in enumerate(zip(sorted_gru_columns, sorted_rounded_percentages)):
        ax = axes.flatten()[i]
        sub_agg_ddb_merged_gdf.plot(column=col, ax=ax, cmap=cmap, vmin=0.01, vmax=1)

        # Plot the dissolved basin boundary on top
        dissolved_basin.plot(ax=ax, edgecolor='black', linewidth=1)

        # Convert title to two lines
        title = split_title(col)
        ax.set_title(title, fontsize=12, ha='center')
        ax.set_xticks([])
        ax.set_yticks([])

        # Display the pre-calculated percentage with a percentage sign in the top left corner of the subplot
        ax.text(0.05, 0.95, f"{rounded_percentage}%", transform=ax.transAxes, fontsize=12,
                verticalalignment='top', bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="none", alpha=0.5))

    # Hide any unused subplots
    for i in range(len(sorted_gru_columns), axes.flatten().shape[0]):
        axes.flatten()[i].set_visible(False)

    # Adjust layout and place colorbar
    fig.tight_layout(pad=1.0)
    cbar_ax = fig.add_axes([0.99, 0.15, 0.03, 0.7])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=16)

    # Save the figure
    plt.savefig(save_path)

    # Show the plot
    plt.show()
