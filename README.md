# Global Climate Statistical Analysis Library (GCSAL)
GCSAL is a software package in MATLAB that allows the user to
- automatically download the Integrated Global Radiosonde Archive (IGRA) raw
text data from the NOAA website
- efficiently process and save the data in a h5 hierarchical file
- quickly access the data, aggregate statistics, and generate plots

## Examples
The following example files are provided in the root directory.  
* **GCSAL_ex1.m**, **GCSAL_ex2.m**, etc. go through the most commonly used functions
of the Global Climate Statistical Analysis Library class. These examples include
analysis of the probability distributions of different atmospheric data based on
location and time of day or time of year.
* **IGRA_to_h5_example.m** takes the GCSAL process from downloading the source
data from NOAA website through to creating a GCSAL Matlab object.  These steps
are only necessary if you want to download and re-build the GCSAL library from
the original source data. Typically this is not necessary as you can use the
provided .h5 file on the website. Reasons you might want to do this are 1) want
to update the data with the latest measurements or 2) want to make a change to the
way the data is stored in the .h5 file.

## Requirements
GCSAL requires MATLAB 2017 and runs on Mac, Linux, or Windows OS

## Building GCSAL
No building required.

## Installing GCSAL
No installation is required.  The user needs to download two data files,
*gcsal.h5* and *gcsal.h5.info.mat*, and put them in the *h5_data* directory.
- gcsal.h5: https://www.dropbox.com/s/2m3glr0drhds33l/gcsal.h5?dl=0
- gcsal.h5.info.mat: https://www.dropbox.com/s/ks9fs3xombb9xqs/gcsal.h5.info.mat?dl=0

## How GCSAL works
GCSAL consists of 2 main groups of functions.  The first group is for processing
the text files efficiently and store the data in the h5 hierarchical file
format.  The second group of functions is for query the H5 data file and create
maps and statistical plots.

## Full documentation
The Global Climate Statistical Analysis Library (GCSAL) allows one to view
climate statistics formulated from over 60 years of radiosonde data from weather
balloons launched at more than 3000 locations around the world! The GCSAL efficiently
processes and compresses the 80 GB of raw text data into 17 GB of data in the h5
hierarchical file format. It provides a simple MATLAB interface to access the
vast quantities of climate data. One can view statistical distributions and
perform statistical operations on the following quantities from sea level to
30 km altitude: wind speed, wind direction, temperature, pressure, dewpoint
depression, and relative humidity.

### List of functions
- Text to Mat to H5
  - IGRA
    - format_definitions
    - datafile2mat
    - datafile2mat_dir
    - mat2h5
    - mat2h5_dir
    - Param
    - Methods
      - read_columns
      - h5write
      - h5write_param
      - Static
        - pad_left
        - convert_to_min_int
        - txt2data
        - compare_bytes_unique
        - str2int
        - char2numerals
        - str2float
        - bits2ints
        - ints2bits
        - compress_idx
        - uncompress_idx
      - Static (Private)
        - unique_inverse
        - compress_txt
        - remove_bad_vals
  - H5
    - create_and_write
    - fullpath
    - load
    - recursive_load
- Query H5 and histograms
  - GCSAL
    - Methods
      - counts
      - counts2
      - countsN
      - query
      - stations_search
      - plot_world_map
      - find_countries
      - find_stations
      - find_headers
      - find_def
      - clear_entries_cache
    - Methods (Private)
      - add_header_params_to_entries
      - stations_from_latlong
      - stations_from_countries
      - stations_from_regex
      - load_all_headers
      - load_from_stations
      - load_group
      - load_param
      - read_from_cached_entries
      - cache_param
      - add_entry_idx_to_headers
    - Static
      - filter_data_by_range
      - histcounts
      - histcounts2
      - histcountsN
      - counts2pdf
      - get_bin_centers
      - find_keys
      - h5info_find_children
      - plot_stations
      - count_and_plot_entries
      - default_bin_edges
      - get_label
      - description_from_filters
      - stations_intersect
      - stations_union
      - stations_setxor
    - Static (Private)
      - struct_set_operation
      - initialize_stations
      - get_entry_idx_in_range
      - header_to_entry_idx
      - station_id_str
  - Map
    - find_in_lat_long_range
    - find_nearest
    - map_stations_by_country
    - multipatch
    - world_map
    - inpolygon


## Join the GCSAL community
* POC: Greg Katz (<gregbkatz@gmail.com>) and David Liu (<zliu@fb.com>)

## License
GCSAL is MIT-licensed.

NOAA Integrated Global Radiosonde Archive (IGRA) data is licensed under the
World Meteorological Organization (WMO) Resolution 40 NOAA Policy NCEI data and
products.
