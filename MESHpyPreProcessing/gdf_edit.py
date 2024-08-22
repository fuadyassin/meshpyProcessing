import geopandas as gpd
import pandas as pd

def flag_ncaalg(
    gdf1: gpd.GeoDataFrame,
    gdf2: gpd.GeoDataFrame,
    threshold: float = 0.1,  # Threshold set to 10% by default
    output_path: str = None
) -> gpd.GeoDataFrame:
    """
    Improve spatial analysis with `assign_intersection_flag_gdf_sindex` function.

    Introduce the `assign_intersection_flag_gdf_sindex` function to efficiently identify
    hydrological polygons in a GeoDataFrame that intersect with non-contributing area polygons
    beyond a specified threshold percentage. This enhancement uses spatial indexing for optimized
    performance, significantly reducing computation time for large datasets. The function adds
    an "ncontr" column to the primary GeoDataFrame, initially marking all polygons with a value of 1.
    Polygons representing non-contributing areas that meet the intersection threshold are then
    updated to a value of 2. This update facilitates the distinction of critical hydrological features
    based on their spatial relationship with non-contributing zones, enabling more nuanced environmental
    and hydrological analyses.

    Parameters
    ----------
    gdf1 : gpd.GeoDataFrame
        The first GeoDataFrame.
    gdf2 : gpd.GeoDataFrame
        The second GeoDataFrame.
    threshold : float, optional
        The threshold for considering an intersection significant, as a fraction of
        the first GeoDataFrame's polygon area (default is 0.1 for 10%).
    output_path : str, optional
        Path where the modified first GeoDataFrame should be saved. If None, the file is not saved.

    Returns
    -------
    gpd.GeoDataFrame
        The modified GeoDataFrame of the first GeoDataFrame with the 'ncontr' column added.
    """
    # Initialize the 'ncontr' column to 1 for all rows in gdf1 as default
    gdf1['ncontr'] = 1
    
    # Create spatial index for gdf2
    spatial_index = gdf2.sindex
    
    # Iterate over gdf1 using spatial indexing to find potential intersections
    for index, row in gdf1.iterrows():
        # Use spatial index to find potential intersections
        possible_matches_index = list(spatial_index.query(row['geometry'], predicate='intersects'))
        if not possible_matches_index:
            continue  # No intersections, move to next row
        
        # Filter possible matches for actual intersection
        possible_matches = gdf2.iloc[possible_matches_index]
        actual_intersections = possible_matches[possible_matches.intersects(row['geometry'])]
        
        # Calculate area fractions for actual intersections
        for _, match in actual_intersections.iterrows():
            intersection_area = row['geometry'].intersection(match['geometry']).area
            area_fraction = intersection_area / row['geometry'].area
            if area_fraction > threshold:
                gdf1.at[index, 'ncontr'] = 2  # Update 'ncontr' to 2 for significant intersections
                break  # Break after finding the first significant intersection
                
    # Save the modified gdf1 to a new shapefile if an output path is provided
    if output_path is not None:
        gdf1.to_file(output_path)
    
    return gdf1

def flag_ncaalg_from_files(
    shapefile1: str,
    shapefile2: str,
    threshold: float = 0.1,  # Threshold set to 10% by default
    output_path: str = None
) -> gpd.GeoDataFrame:
    """
    Read two shapefiles, set their CRS to EPSG:4326, and apply the `flag_ncaalg` function.

    Parameters
    ----------
    shapefile1 : str
        Path to the first shapefile.
    shapefile2 : str
        Path to the second shapefile.
    threshold : float, optional
        The threshold for considering an intersection significant, as a fraction of
        the first GeoDataFrame's polygon area (default is 0.1 for 10%).
    output_path : str, optional
        Path where the modified first GeoDataFrame should be saved. If None, the file is not saved.

    Returns
    -------
    gpd.GeoDataFrame
        The modified GeoDataFrame of the first GeoDataFrame with the 'ncontr' column added.
    """
    # Read the shapefiles into GeoDataFrames
    gdf1 = gpd.read_file(shapefile1)
    gdf2 = gpd.read_file(shapefile2)

    # Set the CRS to EPSG:4326 in place
    gdf1.to_crs(epsg=4326, inplace=True)
    gdf2.to_crs(epsg=4326, inplace=True)

    # Call the original flag_ncaalg function and will return the output of the original fun
    return flag_ncaalg(gdf1, gdf2, threshold, output_path)
