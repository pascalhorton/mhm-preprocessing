# Collection of scripts to prepare the data for the mHM hydrological model.

These were designed for Swiss datasets, but can be adapted to other contexts.

The provided scripts are:
- extract_meteo_data.py : Extraction of the meteorological forcing data.
- flowacc_hydrosnap.py : Computation of flow accumulation using hydro-snap (https://github.com/pascalhorton/hydro-snap) that proceeds to a DEM reconditioning prior to the calculation of the flow accumulation.
- flowacc_pysheds.py : Computation of flow accumulation using pysheds, without DEM reconditioning.
- prepare_morpho_grids.py : Preparation of morphological grids required by mHM.

