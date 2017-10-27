#!/bin/bash
# Use wget to download all data
wget -m ftp://ftp.ncdc.noaa.gov/pub/data/igra/data/data-por/

# cd to download directory
cd ftp.ncdc.noaa.gov/pub/data/igra/data/data-por || exit

# unzip all
find ./ -name \*.zip -exec unzip -n {} \;
