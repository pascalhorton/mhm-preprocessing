import math
import rasterio
import rasterio.features
from rasterio.transform import rowcol
from shapely.geometry import Point
from shapely.geometry import mapping
from pysheds.grid import Grid
import geopandas as gpd
import numpy as np
from pathlib import Path

# Paths
DEM_PATH = R'C:\Data\Projects\2024 Broye mHM\Data\GIS\DEM\DEM_Broye_50m.tif'
STREAMS_SHP = R'C:\Data\Projects\2024 Broye mHM\Data\GIS\River network\River_network_tlm3d_v2.shp'
OUTPUT_DIR = R'C:\Data\Projects\2024 Broye mHM\Analyses\Data preprocessing'
OUTLET_SHP = R'C:\Data\Projects\2024 Broye mHM\Data\GIS\Catchment\Outlet.shp'
DELTA = 0.01  # Elevation difference to lower the next pixel when correcting

output_dir = Path(OUTPUT_DIR)
output_dir.mkdir(exist_ok=True)

dem = rasterio.open(DEM_PATH)

# Correct the DEM using Pysheds
pysheds_grid = Grid.from_raster(str(DEM_PATH))
pysheds_dem = pysheds_grid.read_raster(str(DEM_PATH))
pit_filled_dem = pysheds_grid.fill_pits(pysheds_dem)
flooded_dem = pysheds_grid.fill_depressions(pit_filled_dem)
inflated_dem = pysheds_grid.resolve_flats(flooded_dem)

# Compute flow accumulation
fdir = pysheds_grid.flowdir(inflated_dem, nodata_out=np.int64(0))
acc = pysheds_grid.accumulation(fdir, nodata_out=np.float64(-9999))

if OUTLET_SHP:
    # Load the outlet
    outlet = gpd.read_file(OUTLET_SHP)
    (x, y) = (outlet.geometry.x[0], outlet.geometry.y[0])

    # Snap the outlet to the nearest cell with a high flow accumulation
    x_snap, y_snap = pysheds_grid.snap_to_mask(acc > 10000, (x, y))

    # Compute the catchment
    catchment = pysheds_grid.catchment(x=x_snap, y=y_snap, fdir=fdir)

    # Save the catchment
    output_catchment_path = output_dir / 'catchment.tif'
    with rasterio.open(output_catchment_path, 'w', **dem.profile) as dst:
        dst.write(catchment, 1)

# Save the flow direction
output_fdir_path = output_dir / 'flow_direction.tif'
with rasterio.open(output_fdir_path, 'w', **dem.profile) as dst:
    dst.write(fdir, 1)

# Save the flow accumulation
output_acc_path = output_dir / 'flow_accumulation.tif'
with rasterio.open(output_acc_path, 'w', **dem.profile) as dst:
    dst.write(acc, 1)

dem.close()
