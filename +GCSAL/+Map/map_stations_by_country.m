function countries = map_stations_by_country(all_stations, countries_data)
% countries_with_station = map_stations_by_country(all_stations, countries_data)
%   For each country in countries_data, finds stations in all_stations that
%   are within the borders of the country based on the latitude and
%   longitude of the stations and countries.
%
%   Returns a struct array countries where each element contains the
%   country name and a list of stations ids that were found to be contained
%   by that country.
%
%   If countries_data is not input, then the function will try to load
%   countries_data from ne_110m_admin_0_countries.mat if the file exists on
%   the path
%
% INPUTS
%     all_stations - struct array, each element should contain lat and long
%                    of the station in degrees and id for the identifier of
%                    the station
%   countries_data - struct array, each element should contain Lat and Lon of
%                    the borders of the country in degrees as well as name to
%                    identify the country
%
% OUTPUTS
%    countries - struct array, each element contains name, Lat, Lon, and
%                stations. name is a string. Lat and Lon are the borders
%                of the country in degrees. stations is a string matrix
%                with each row an id for a station contained in that
%                country.

% Copyright (c) 2017-present, Facebook, Inc.
% All rights reserved.
%
% This source code is licensed under the BSD-style license found in the
% LICENSE file in the root directory of this source tree. An additional grant
% of patent rights can be found in the PATENTS file in the same directory.


% If countries_data not given, try to load it form default .mat file
if ~exist('countries_data', 'var')
    fname = 'ne_110m_admin_0_countries.mat';
    if exist(fname, 'file')
        fprintf('Loading countries data from %s\n', which(fname))
        countries_data = load(fname);
        countries_data = countries_data.worldData;
    else
        fprintf('Could not load country data. %s file not found\n', fname)
        countries = struct([]);
        return
    end
end


% Some data sets use name and others use NAME
if isfield(countries_data, 'name')
    names = {countries_data.name};
elseif isfield(countries_data, 'NAME')
    names = {countries_data.NAME};
else
    error('countries_data is missing the name field')
end

% Initialize output with Lat, Lon, and name fields from countries_data
countries = struct('Lat', {countries_data.Lat}, ...
                   'Lon', {countries_data.Lon}, ...
                   'name', names, ...
                   'stations', '');


% Extract vector from all_stations struct array
stations_lat = [all_stations.lat];
stations_long = [all_stations.lon];
stations_id = vertcat(all_stations.id);

% Loop through each country in countries_data
for i = 1:length(countries)

    % Get logical array corresponding to stations found in country
    in_idx = GCSAL.Map.inpolygon2(stations_long, stations_lat, ...
        countries(i).Lon, countries(i).Lat);

    % Record station ids found in the country
    countries(i).stations = stations_id(in_idx,:);

end


end
