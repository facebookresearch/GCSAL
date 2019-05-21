% Copyright (c) Facebook, Inc. and its affiliates.
%
% This source code is licensed under the MIT license found in the
% LICENSE file in the root directory of this source tree.
%
% This example file goes through the most commonly used function of the
% Global Climate Statistical Analysis Library class. The examples show you
% can do analysis of the probability distributions of different atmospheric
% data based on location and time of day or time of year.


%% Step 1 is to download the source .h5 file and set up your paths correctly

clear;
close all;
clc;

if (~exist('g','var') || ~isa(g,'GCSAL.GCSAL'))
    % The gcsal.h5 and gcsal.h5.info.mat files are available for download from the
    % website and should be placed in the h5_data directory.

    %%%%%%%%%%%%%% CHANGE THESE %%%%%%%%%%%%%%%%%%
    % Set this to wherever you put the gcsal.h5 file and gcsal.h5.info.mat
    % files downloaded from dropbox
    h5_dir = './h5_data/';

    % Directory to code. The folder +GCSAL which contains this file should be
    % in this directory
    codebase_dir = './';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Full path to .mat file with h5 info
    h5_file = fullfile(h5_dir, 'gcsal.h5');
    h5_mat_file = [h5_file '.info.mat'];

    % Set up Matlab path
    addpath(genpath(codebase_dir))

    %% Load GCSAL object from .mat file

    % This requires about 6 gb of RAM just to hold all of the header
    % information in memory

    % Normally you should load the GCSAL object from the .mat file but if it
    % doesn't exist on your path you can use the .h5 file. After using the .h5
    % file a .mat file will be created automatically for subsequent use
    if ~exist(h5_mat_file, 'file')
        g = GCSAL.GCSAL(h5_file);
    else
        g = GCSAL.GCSAL(h5_mat_file);
        g.h5_fname = h5_file;
    end
end

% Introduction: Printout the list of variables that are in the header data
% and the entries data.
% header data is for a single balloon launch - things like time, date, and location
% entries data is the measurements of the baloon - wind speed, pressure, etc.
g.defs.header.params
g.defs.entries.params

% For more details look at a single parameter
g.defs.entries.params.wspd
g.defs.entries.params.gph


%% Stations_search

% Find stations within 25 degrees of the equator
stations1 = g.station_search('LatLong', [-25 25 -180 180]);

% Find stations from a bunch of countries
country_names = {'Mexico', 'Brazil','Algeria', 'Burkina Faso', 'Ghana', 'Niger',...
           'Nigeria',  'Egypt', 'Sudan', 'Ethiopia', 'Uganda',...
           'Kenya', 'Tanzania', 'Madagascar', 'India', 'Sri Lanka',...
           'Nepal', 'Bangladesh','Myanmar', 'Thailand', 'Vietnam', 'Cambodia',...
           'Ukraine', 'Uzbekistan', 'Turkey',...
           'Indonesia', 'Philippines'};

stations2 = g.station_search('Countries', country_names);

% Find stations in countries AND within 25 degrees of equator
stations3 = g.station_search('Countries', country_names, ...
    'LatLong', [-25 25 -180 180]);

% Find stations in countries OR within 25 degrees of equator
stations4 = GCSAL.GCSAL.stations_union(stations1, stations2);

% Plot stations 4 on wolrd map
figure; hold all; g.plot_world_map();
GCSAL.GCSAL.plot_stations(stations4, 'r+');

% Find stations with IDs beginnign with the letter A.
% Note that in regex ^ means beginning of the line
stations5 = g.station_search('IDRegex', '^A');

% Find stations in Brazil or India AND within 25 degrees of
% the equator AND with station IDs ending in 5
% Note that in regex $ means end of the line
stations6 = g.station_search('Countries', {'Brazil', 'India'}, ...
    'LatLong', [-25 25 -180 180], ...
    'IDRegex', '5$');

% Find stations within +/-2.5 deg latitude about Guatemala
stations7 =  g.station_search('Lat', [14.583323, -90.527309], 'Range', 2.5);

%% Data query

% Get stations located in Botswana
stations = g.station_search('Countries', {'Botswana'});

% Get all geopotential height and  windspeed data along with hour, month,
% and year data
entries1 = g.query(stations, {'gph', 'wspd', 'hour', 'month', 'year'});

% Plot distribution of hours and years for the data in entries1
figure; histogram(vertcat(entries1.hour)); xlabel('Hour'); ylabel('# of occurences');
figure; histogram(vertcat(entries1.year)); xlabel('Year'); ylabel('# of occurences');

% Get gph and wspd data measured between 6 and 4 pm
entries2 = g.query(stations, {'gph', 'wspd'}, 'hour', [6 16]);

% Plot distribution of hours for the data in entries2
figure; histogram(vertcat(entries2.hour)); xlabel('Hour'); ylabel('# of occurences');

% Get data corresponding only to measuresments taken in August between
% 4am and Noon and in the years 1990 to 1999
entries3 = g.query(stations, {'gph', 'wspd'}, ...
                  {'month', 'hour', 'year'}, ...
                  {8, [4 12], [1990 1999]});

% Plot distribution of years for the data in entries3
figure; histogram(vertcat(entries3.year)); xlabel('Year'); ylabel('# of occurences');

% To save RAM clear out the entries cache. The entries cache holds data for
% any data that has been loaded so far which makes it faster to the access
% the data subsequently but also uses up RAM;
g.clear_entries_cache();

%% One dimensional histograms

% Get stations in Brazil
stations = g.station_search('Countries', 'Brazil');

% Histogram for all windspeeds
[N, entries] = g.counts(stations, 'wspd');

% Define custom bin_edges
bin_edges = 0:1:80;

% Do counts for windspeed filtered on geopotential
% altitude between 20 and 30 km and with custom bin edges
[N, entries] = g.counts(stations, 'wspd', 'FilterFields', {'gph'}, ...
     'FilterRanges', {[20 30]}, 'Edges', bin_edges);

% Now additionally filter on measurements taken in August
% between 4 and 10 am
[N, entries] = g.counts(stations, 'wspd', 'FilterFields', {'gph', 'month', 'hour'}, ...
    'FilterRanges', {[20 30], 8, [4 10]}, 'Edges', bin_edges);

%% Two dimensional histograms

% Get some stations
stations = g.station_search('Countries', 'Brazil');

% Do counts between gph and various other parameters
[N, entries] = g.counts2(stations, 'gph', 'wspd');
[N, entries] = g.counts2(stations, 'gph', 'press');
[N, entries] = g.counts2(stations, 'gph', 'temp');
[N, entries] = g.counts2(stations, 'gph', 'rh');
[N, entries] = g.counts2(stations, 'gph', 'dpdp');
[N, entries] = g.counts2(stations, 'gph', 'wdir');
[N, entries] = g.counts2(stations, 'gph', 'wspd');

% Do counts between gph and pressure with custom bin edges
[N, entries] = g.counts2(stations, 'gph', 'wspd', ...
   'XEdges', 0:0.5:40, 'YEdges', 0:1:80);

% Do counts for data measured between 6am and 4pm and in August
[N, entries] = g.counts2(stations, 'gph', 'wspd', ...
'FilterFields', {'hour', 'month'}, 'FilterRanges', {[6 16], 8});


%% N dimensional counts

% Get some stations
stations = g.station_search('Countries', 'Brazil');

% Make 5-dimensional count matrix with default bin edges and no filtering
resolutions = struct('fld', {'lat', 'lon', 'gph', 'month', 'wspd'});
N = g.countsN(stations, resolutions);

% Add custom bin edges to gph field and limit data to only data between 6
% and 10 am
resolutions(3).edges = 0:1:80;
N = g.countsN(stations, resolutions, 'FilterFields', {'hour'}, 'FilterRanges', [6 10]);

% Clear the cache to save RAM
g.clear_entries_cache();

%% N dimensional counts on full library
% This section has been commented out because it takes ~35 minutes and 30
% gb to complete

% % Get all stations
% stations = g.stations;
%
% % Cache data for lat, lon, gph, and month
% % This may take a few minutes and requires another 18 gb of RAM for a total
% % of 24 gb of RAM including the header data that is loaded when the GCSAL
% % object is initialized.
% g.query(stations, {'lat', 'lon', 'gph', 'month'});
%
% % Turn off cache so we don't use any more RAM going forward
% g.do_cache_entries = false;
%
% % Make 5-dimensional count matrix with default bin edges and no filtering
% % This may take about 5 minutes for each countsN call.
% % Also temporarily uses another 6 gb of RAM for each countsN call for a
% % grand total of ~30gb of RAM. As long as do_cache_entries is false this
% % last 6 gb RAM is cleared between each call to countsN
% resolutions = struct('fld', {'lat', 'lon', 'gph', 'month', 'wdir'});
% N1 = g.countsN(stations, resolutions); % This may take about 5 minutes
%
% resolutions = struct('fld', {'lat', 'lon', 'gph', 'month', 'wspd'});
% N2 = g.countsN(stations, resolutions); % This may take about 5 minutes
%
% resolutions = struct('fld', {'lat', 'lon', 'gph', 'month', 'temp'});
% N3 = g.countsN(stations, resolutions); % This may take about 5 minutes
%
% resolutions = struct('fld', {'lat', 'lon', 'gph', 'month', 'press'});
% N4 = g.countsN(stations, resolutions); % This may take about 5 minutes
%
% resolutions = struct('fld', {'lat', 'lon', 'gph', 'month', 'dpdp'});
% N5 = g.countsN(stations, resolutions); % This may take about 5 minutes
%
% resolutions = struct('fld', {'lat', 'lon', 'gph', 'month', 'rh'});
% N6 = g.countsN(stations, resolutions); % This may take about 5 minutes
