#!/bin/bash

# Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved

# Use wget to download all data
wget -m ftp://ftp.ncdc.noaa.gov/pub/data/igra/data/data-por/

# cd to download directory
cd ftp.ncdc.noaa.gov/pub/data/igra/data/data-por || exit

# unzip all
find ./ -name \*.zip -exec unzip -n {} \;
