# xarray and dask are required
import xarray as xr
import rioxarray as rxr

BASE_DIR = R'\\hydroshare.giub.unibe.ch\data\Meteorology\Switzerland\MeteoSwiss_gridded_products'
DIR_PRECIP = f'{BASE_DIR}/RhiresD_v2.0_swiss.lv95'
DIR_TMEAN = f'{BASE_DIR}/TabsD_v2.0_swiss.lv95'
DIR_TMAX = f'{BASE_DIR}/TmaxD_v2.0_swiss.lv95'
DIR_TMIN = f'{BASE_DIR}/TminD_v2.0_swiss.lv95'
OUTPUT_DIR = R'C:\Data\Projects\2024 Broye mHM\Data'
# MASK_PATH = R'C:\Data\Projects\2024 Broye mHM\Data\GIS\Mask_L2.tif'
MASK_PATH = None

Y_START = 1993
Y_END = 2022

DOMAIN = [2542000, 1148000, 2571000, 1199000]  # [x_min, y_min, x_max, y_max]


def main():
    # Export precip
    export_mch_data(DIR_PRECIP, OUTPUT_DIR, DOMAIN, Y_START, Y_END, MASK_PATH,
                    'RhiresD_ch01h.swiss.lv95', 'precip',
                    'RhiresD', 'pre')
    print('Precipitation data exported')

    # Export temperature
    export_mch_data(DIR_TMEAN, OUTPUT_DIR, DOMAIN, Y_START, Y_END, MASK_PATH,
                    'TabsD_ch01r.swiss.lv95', 'tmean',
                    'TabsD', 'tavg')
    print('Temperature data exported')

    # Export temperature max
    export_mch_data(DIR_TMAX, OUTPUT_DIR, DOMAIN, Y_START, Y_END, MASK_PATH,
                    'TmaxD_ch01r.swiss.lv95', 'tmax',
                    'TmaxD', 'tmax')
    print('Temperature max data exported')

    # Export temperature min
    export_mch_data(DIR_TMIN, OUTPUT_DIR, DOMAIN, Y_START, Y_END, MASK_PATH,
                    'TminD_ch01r.swiss.lv95', 'tmin',
                    'TminD', 'tmin')
    print('Temperature min data exported')

    print('All data exported')


def export_mch_data(dir_data, dir_output, domain, year_start, year_end, mask_path,
                    file_name_input, file_name_output, var_name_input, var_name_output):
    # Load data
    file_paths = [
        f"{dir_data}/{file_name_input}_{year}01010000_{year}12310000.nc" for year
        in range(year_start, year_end + 1)]
    ds = xr.open_mfdataset(file_paths, decode_times=False)

    # Select domain
    ds = ds.sel(E=slice(domain[0], domain[2]),
                N=slice(domain[1], domain[3]))

    # Set variable to double precision
    ds[var_name_input] = ds[var_name_input].astype('float64')

    # Set time as integer
    ds['time'] = ds['time'].astype('int64')

    # Drop unnecessary variables (lat, lon, swiss_lv95_coordinates)
    ds = ds.drop_vars(['lat', 'lon', 'swiss_lv95_coordinates'])

    # Rename variables
    ds = ds.rename({var_name_input: var_name_output, 'E': 'x', 'N': 'y'})

    # Reverse y-axis to match expected orientation
    ds = ds.reindex(y=list(reversed(ds.y)))

    # Load mask
    if mask_path:
        mask = rxr.open_rasterio(mask_path)
        mask = mask.squeeze().drop_vars('band')

        # Mask data
        ds = ds.where(mask == 1)

    # Replace NaN with -9999 for the main variable
    ds[var_name_output] = ds[var_name_output].where(ds[var_name_output].notnull(), -9999)

    ds.to_netcdf(f'{dir_output}/{file_name_output}_{year_start}-{year_end}.nc')


if __name__ == '__main__':
    main()
