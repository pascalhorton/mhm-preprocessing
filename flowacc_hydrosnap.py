from hydro_snap import recondition_dem

# Paths
DEM_PATH = R'C:\Data\Projects\2024 Broye mHM\Data\GIS\DEM\DEM_Broye_50m.tif'
STREAMS_SHP = R'C:\Data\Projects\2024 Broye mHM\Data\GIS\River network\River_network_tlm3d_v2.shp'
OUTPUT_DIR = R'C:\Data\Projects\2024 Broye mHM\Analyses\Data preprocessing'
OUTLET_SHP = R'C:\Data\Projects\2024 Broye mHM\Data\GIS\Catchment\Outlet.shp'
CATCHMENT_SHP = R'C:\Data\Projects\2024 Broye mHM\Data\GIS\Catchment\Broye_total_v2.shp'
BREACHES_SHP = R'C:\Data\Projects\2024 Broye mHM\Data\GIS\River network\Breaches.shp'
DELTA = 0.01  # Elevation difference to lower the next pixel when correcting

# Recondition the DEM
recondition_dem(DEM_PATH, STREAMS_SHP, OUTPUT_DIR, delta=DELTA, outlet_shp=OUTLET_SHP, catchment_shp=CATCHMENT_SHP,
                breaches_shp=BREACHES_SHP)
