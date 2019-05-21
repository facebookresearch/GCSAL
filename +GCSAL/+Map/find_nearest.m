function [stations, arclen] = find_nearest(all_stations, lat, lon, n)
% Copyright (c) Facebook, Inc. and its affiliates.
%
% This source code is licensed under the MIT license found in the
% LICENSE file in the root directory of this source tree.
%
% [stations, arclen] = find_nearest(all_stations, lat, lon, n)
%
%   Returns an array of station structs of the nearest n stations
%   relative to the specified lat/lon
%
%   lat and lon must be single values each in degrees.
%
% INPUTS
%   all_stations - struct array, each elemetn contains lat, long, id
%            lat - scalar defining reference latitude in degrees.
%            lon - scalar defining reference longitude in degrees.
%              n - number of nearest stations
%
% OUTPUTS
%       stations - struct array, subset of all_stations located within the
%                  lat/long search box
%         arclen - vector of distances in meters


    lats = [all_stations(:).lat];
    lons = [all_stations(:).lon];
    E = referenceEllipsoid('wgs84');
    [arclen, ~] = distance(lats, lons, lat, lon, E);
    [~,idx] = sort(arclen);
    stations = all_stations(idx(1:n));
    arclen = arclen(idx(1:n));

end
