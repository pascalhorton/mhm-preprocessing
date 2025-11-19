import rasterio
import numpy as np
import geopandas as gpd
from pathlib import Path

CATCHMENT_PATH = R'C:\Users\Pascal Horton\Documents\Data\Projects\2024 Broye mHM\Data\GIS\Catchment\Broye_total_v2.shp'
DEM_PATH = R'C:\Users\Pascal Horton\Documents\Data\Projects\2024 Broye mHM\Data\Analyses\Data preprocessing\Python_50m\corrected_dem_final.tif'
FLOWDIR_PATH = R'C:\Users\Pascal Horton\Documents\Data\Projects\2024 Broye mHM\Data\Analyses\Data preprocessing\Python_50m\flow_direction.tif'
FLOWACC_PATH = R'C:\Users\Pascal Horton\Documents\Data\Projects\2024 Broye mHM\Data\Analyses\Data preprocessing\Python_50m\flow_accumulation.tif'
OUTPUT_DIR = R'C:\Users\Pascal Horton\Documents\Data\Projects\2024 Broye mHM\Data\Analyses\Data preprocessing'


def main():
    output_dir = Path(OUTPUT_DIR)

    # Load the catchment shape file
    catchment = gpd.read_file(CATCHMENT_PATH)

    # Export existing rasters as ASCII
    export_as_ascii(DEM_PATH, output_dir, 'dem.asc')
    export_as_ascii(FLOWDIR_PATH, output_dir, 'fdir.asc')

    # Flow accumulation needs to be adjusted (pyshed flow acc - 1)
    flowacc = rasterio.open(FLOWACC_PATH)
    flowacc_data = flowacc.read(1)
    flowacc_data -= 1
    save_as_ascii(flowacc_data, output_dir, 'facc.asc', flowacc)

    # Load the DEM
    dem_raster = rasterio.open(DEM_PATH)

    # Compute the slope
    slope = compute_slope_percent_rise(dem_raster)
    save_as_ascii(slope, output_dir, 'slope.asc', dem_raster)

    # Compute the aspect
    aspect = compute_aspect_degrees(dem_raster)
    save_as_ascii(aspect, output_dir, 'aspect.asc', dem_raster)

    print('Done!')


def compute_aspect_degrees(dem_raster):
    dem = dem_raster.read(1)
    transform = dem_raster.transform
    x, y = transform[0], -transform[4]
    aspect = np.arctan2(-np.gradient(dem, x, axis=1), np.gradient(dem, y, axis=0))
    aspect = np.degrees(aspect) % 360

    return aspect


def compute_slope_percent_rise(dem_raster):
    dem = dem_raster.read(1)
    transform = dem_raster.transform
    x, y = transform[0], -transform[4]
    slope_percent_rise = np.arctan(np.sqrt(np.gradient(dem, x, axis=1) ** 2 + np.gradient(dem, y, axis=0) ** 2))

    # Compute the gradients along the x and y axes
    gradient_x = np.gradient(dem, x, axis=1)
    gradient_y = np.gradient(dem, y, axis=0)

    # Compute the slope as the square root of the sum of the squares of the gradients
    slope_percent_rise = np.sqrt(gradient_x**2 + gradient_y**2) * 100

    return slope_percent_rise


def save_as_ascii(array, output_dir, output_name, ref_raster):
    profile = ref_raster.profile
    profile.update(driver='AAIGrid')
    output_path = output_dir / output_name
    with rasterio.open(output_path, 'w', **profile) as dst:
        dst.write(array, 1)


def export_as_ascii(raster_path, output_dir, output_name):
    raster = rasterio.open(raster_path)
    data = raster.read(1)
    profile = raster.profile
    profile.update(driver='AAIGrid')
    output_path = output_dir / output_name
    with rasterio.open(output_path, 'w', **profile) as dst:
        dst.write(data, 1)


if __name__ == '__main__':
    main()
