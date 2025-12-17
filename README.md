# Generic Gridding Routine
One single Julia code `gridL2_Dates.jl` to grid all kinds of satellite data onto rectangular grid at arbitrary spacing (and determined spatial and temporal resolution). Input defined just via json files, which makes it generic as long as you have corner coordinates for your satellite dataset (you can grid XCO2, methane CO, etc from TROPOMI or SIF from OCO-2, OCO-3, TROPOMI, GOME-2, etc)..

The main program, gridL2_Dates.jl, can compute gridded averages and save everything into a netCDF4 file that can be read via tools such as python, Julia, Panoply, etc ( 3D dataset with time/lat/lon). Most importantly, it can oversample, i.e. the gridding takes the actual footprint overlap with the final grid into account. This is done by splitting a footprint via its corner coordinates into a set of points nxn within that footprint (typicall n=10, i.e. 100 points). For these, a simple "in-the-box" gridding routine will be applied, i.e. tha algorithm finds in which grid boxe each sub-pixel of a footprint falls. This results in fractional "samples" per grid box.

The program was written as a first step into Julia (as everything else was too slow and C++ development time just too long). So, some things are still clumsy and not really "Julian". Feel free to make things more elegant and create a pull request if you have done so. 

## Install Julia

The current version works well with Julia versions > 1.0 but I didn't test all of them. Download it for your platform from [Julia's old
releases](https://julialang.org/downloads/oldreleases/#v131_dec_30_2019).

You might have to install a couple of Julia packages, this can done like in Julia using the built-in package manager (press `]` at the Julia prompt):

```julia
julia> ]
(v1.3) pkg> add ArgParse, Statistics, Glob, JSON, Dates, Printf
```

## How to run the program

Hint for first time users, you need to use the console/terminal to use it properly. On a Mac, it is good to add an alias, e.g. 

```
alias julia /Applications/Julia-1.3.app/Contents/Resources/julia/bin/julia
```

in your .bash_profile or whatever you are using.

Once this is done, you can test out what options there are (I often do this myself):

```
$ julia ./gridL2_Dates.jl --help
usage: gridL2_Dates.jl [--Dict DICT] [-o OUTFILE] [--monthly]
                       [--compSTD] [--latMin LATMIN] [--latMax LATMAX]
                       [--lonMin LONMIN] [--lonMax LONMAX]
                       [--dLat DLAT] [--dLon DLON]
                       [--startDate STARTDATE] [--stopDate STOPDATE]
                       [--dDays DDAYS] [-h]

optional arguments:
  --Dict DICT           JSON dictionary file to use (default:
                        "/home/cfranken/code/gitHub/Gridding/gridding/tropomi_all.json")
  -o, --outFile OUTFILE
                        output filename (default OCO2_SIF_map.nc)
                        (default: "OCO2_SIF_map.nc")
  --monthly             Use time-steps in terms of months (not days)
  --compSTD             compute standard deviation within dataset
  --latMin LATMIN       Lower latitude bound (type: Float32, default:
                        -90.0)
  --latMax LATMAX       Upper latitude bound (type: Float32, default:
                        90.0)
  --lonMin LONMIN       Lower longitude bound (type: Float32, default:
                        -180.0)
  --lonMax LONMAX       Upper longitude bound (type: Float32, default:
                        180.0)
  --dLat DLAT           latitude resolution (type: Float32, default:
                        1.0)
  --dLon DLON           longitude resolution (type: Float32, default:
                        1.0)
  --startDate STARTDATE
                        Start Date (in YYYY-MM-DD) (default:
                        "2018-03-07")
  --stopDate STOPDATE   Stop Date (in YYYY-MM-DD) (default:
                        "2018-10-31")
  --dDays DDAYS         Time steps in days (or months if --monthly is
                        set) (type: Int64, default: 8)
  -h, --help            show this help message and exit
```

The most important part of this script is `--Dict` as most description about the files you want to grid are within the json file you will use here. You can also change the spatial resolution (e.g. `--dLat`), the lat/lon boundaries if you don't want the globe and the date ranges and temporal resolution you are interested in (`--startDate`, `--stopDate`, `--dDays` defines the number of days you want to aggregate together (if you set `--monthly`, this will be the step in months). If you set `--compSTD`, the standard deviation of the data within a grid box will be computed (the original variable name will be kept with an `_std` added to it and currently discarding all attributes). This is not necessarily an uncertainty estimate but the true spread of data within a grid box. 

An example on how to run this is

```
julia ./gridL2_Dates.jl --dLat 1 --dLon 1  --dDays 16 --startDate 2019-02-01 --stopDate 2020-01-30  --Dict oco2_all.json  -o ~/oco2_16day_1deg_scf.nc
```

### JSON files
Now to the most important part, the json structure makes the code very general (check out https://www.json.org/json-en.html).
An example file content is below:
```
{
    "basic":{
    "lat": "latitude",
        "lon": "longitude",
        "lat_bnd": "Latitude_Corners",
        "lon_bnd": "Longitude_Corners",
        "time": "time"
},
    "grid":{
        "sif_757": "Daily_SIF_757nm",
        "sif_771": "Daily_SIF_771nm",
        "sif_757_relative": "Science/SIF_Relative_757nm",
        "sif_771_relative": "Science/SIF_Relative_771nm",
        "dcCorr": "Science/daily_correction_factor"
    },
    "filePattern": "YYYY/MM/oco2_LtSIF_??MMDD_*.nc4",
    "folder": "/net/fluo/data2/groupMembers/cfranken/data/kurosu/test_data/daily/oco2/new/"
}
```

The structure will be used to define dictionaries in Julia, with a key (the left side, pointing to internal variables in the code) and a value (this is actually the path to the respetive dataset within the files you want to grid). What is needed for sure is the key/value pair for `lat_bnd` and `lon_bnd` as these corner coordinates in the files are used to perform proper gridding and over-sampling (I think time is not even used right now). The other part is in the `grid` group, eveything you will add here will be gridded, i.e. all variables on the left side of the key/value pairs (and saved how you call it on the left). 

Last but not least, we need a `filePattern`, which defines where to look for the data (`YYYY`, `MM` and `DD` are keywords, which will be internally used to find the right matching years, months and days). Wildcards can be used (`?`, `*`, etc). Then you need to provide the main folder where the data is located (the full path is `folder/filePattern`).

### Filter criteria:
You can add filter criteria as well, as we sometimes want to make sure that quality filters are applied, angles are within a specific range, and so forth. Within the json file, this can be done by adding groups called `filter_eq`, `filter_gt`, or `filter_lt`). These test for equalities, greater than or lower than. An example is here:
```
{
    "basic":{
    "lat": "PRODUCT/latitude",
    "lon": "PRODUCT/longitude",
    "lat_bnd": "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/latitude_bounds",
    "lon_bnd": "PRODUCT/SUPPORT_DATA/GEOLOCATIONS/longitude_bounds",
    "time": "PRODUCT/time_utc"
},
    "grid":{
        "ch4": "PRODUCT/methane_mixing_ratio_bias_corrected",
        "altitude": "PRODUCT/SUPPORT_DATA/INPUT_DATA/surface_altitude"
    },
    "filter_gt":{
        "PRODUCT/qa_value": 0.3,
        "PRODUCT/methane_mixing_ratio_bias_corrected": 1600
    },
    "filter_lt":{
        "PRODUCT/methane_mixing_ratio_bias_corrected": 2200
    },
    "filePattern": "S5P_*_L2__CH4____YYYYMMDD*.zip",
    "folder": "/net/fluo/data2/projects/TROPOMI_GASES/ch4/"
}
```
here, we want to make sure that the `qa_value` is larger than 0.3 and that the methane mixing ratio is within a reasonable range (be careful doing things like this, just added for demonstration).

## How to run the modularized script from within python
Loading the modularized julia script exposes two functions: run_gridding (rectangular grid) and run_gridding_km (sinusoidal km-scale grid). The JSON dictionaries can have the same structure, however, the "folder" key is not used by the two functions as a list of files is expected as input. 
The functions can be used from python using the juliacall package:

```
from juliacall import Main as jl
jl.include("grid_S5L2_Dates_module.jl")

files_in = ["/path/to/your/file_1.nc", "/path/to/your/file_2.nc"]

res = jl.S5Grid.run_gridding(files_in,
    json_dict="jsonFiles/variable_dict.json",
    monthly=False, 
    compSTD=True,
    latMin=35., 
    latMax=60.,  
    lonMin=0., 
    lonMax=45., 
    dLat=0.1, 
    dLon=0.1, 
    startDate_str="2025-12-14", 
    stopDate_str="2025-12-15", 
    dDays=2, 
    oversample_temporal=1.0,
    FillValue=-999
)
```
The output consists of lat, lon, time, averaged_data [time, lon, lat, variables],  n (number of averaged soundings per grid cell) [time, lon, lat], and optionally: std_data [time, lon, lat, variables] (when compSTD is set to true).
       

The difference in the run_gridding_km (sinusoidal km-scale grid) function with respect to original version is that dLat/dLon are not needed as input anymore as both are now controlled by target_km.

```
from juliacall import Main as jl
jl.include("grid_S5L2_Dates_module.jl")

files_in = ["/path/to/your/file_1.nc", "/path/to/your/file_2.nc"]

res = jl.S5Grid.run_gridding_km(files_in,
    json_dict="jsonFiles/variable_dict.json",
    target_km=10.0,
    monthly=False, 
    compSTD=True,
    latMin=35., 
    latMax=60.,  
    lonMin=0., 
    lonMax=45., 
    startDate_str="2015-11-24", 
    stopDate_str="2015-11-26", 
    dDays=2, 
    oversample_temporal=1.0,
    FillValue=-999
)
```

The output is a little bit more complex for this function as we do not have a rectangular grid anymore. Instead of a matrix array[lon, lat], ragged arrays (a vector of vectors) are used. In total the output consists of: lat_centers, n_lons_per_row, time, averaged_data [time][lat_row][lon_bin, var_idx], n (number of averaged soundings per grid cell) [time][lat_row][lon_bin], and optionally: std_data [time][lat_row][lon_bin, var_idx] (when compSTD is set to true). The longitudes are returned in a normalized format and can be converted back by:
```
# Assuming column (c) in a specific row (r)
d_lon = (lon_max - lon_min) / n_lons_per_row[r]
lon_min + (c + 0.5) * d_lon # Center Longitudes 
lon_min + c * d_lon # Longitude boundary 1
lon_min + (c + 1) * d_lon # Longitude boundary 2
````

A full (re-)construction of coordinates for a given variable for plotting would look like this:
``` python
var_idx = 0
patches = []
values  = []
    
lats = res.lat_centers
d_lat = lats[1] - lats[0]
    
for r in range(len(lats)):
    n_bins = res.n_lons_per_row[r]
    d_lon = (lon_max - lon_min) / n_bins
    lat_b = lats[r] - d_lat/2
    lat_t = lats[r] + d_lat/2
        
    for c in range(n_bins):
        val = res.averaged_data[time_step_idx][r][c, var_idx]
        val = np.nan if val == -999 else val
         
        lon_l = lon_min + c * d_lon
        lon_r = lon_min + (c+1) * d_lon
        
        # Define the 4 corners of the bin
        vertices = [(lon_l, lat_b), (lon_r, lat_b), (lon_r, lat_t), (lon_l, lat_t)]
        patches.append(vertices)
        values.append(val)
```


## Code of Conduct:
Please feel free to use this tool but make sure that you help the community if you find bugs, improve it etc. Any modifications that are useful should be made publicly available, you can fork and create a pull request. Also, let us know if you find bugs. On top of that, please acknowledge the tool if you use it in publications.

## MODIS files
We use these as well but I have no time to document all of them right now. Contact us if you want to use them (cfranken"at"caltech.edu).
