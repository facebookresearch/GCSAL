function [stations, latbox, longbox] = ...
    find_in_lat_long_range(all_stations, latrange, longrange)
% [stations, longbox, latbox] = ...
%     find_in_lat_long_range(all_stations, latrange, longrange)
%
%   Returns an array of station structs that fall within the
%   box defined by latrange and longrange from the list of stations in
%   all_stations. all_stations is a struct array with each element
%   containing lat, long, and id fields.
%
%   Additionally returns longbox and latbox which can be used to plot the
%   the searchbox that was used.
%
%   latrange and longrange must be two element vectors and are
%   in degrees. This function accounts for angle wraparound. So latrange
%   could be [-45 45] to find stations in latitudes between -45 and 45
%   degrees or it could be [45 -45] to find stations with latitude above 45
%   deg or below -45.  The same goes for longrange
%
% INPUTS
%   all_stations - struct array, each elemetn contains lat, long, id
%       latrange - two element vector defining latitude limits in degrees.
%                  Angle wraparound is OK.
%      longrange - two element vector defining longitude limits in degrees.
%                  Angle wraparound is OK.
%
% OUTPUTS
%       stations - struct array, subset of all_stations located within the
%                  lat/long search box
%         latbox - vector of latitude values in degrees that can be used to
%                  make a plot representing the search box used. This
%                  handles wraparound nicely for a map plot by inserting
%                  NaNs for a discontinuous plot line.
%        longbox - see latbox

% Copyright (c) 2017-present, Facebook, Inc.
% All rights reserved.
%
% This source code is licensed under the BSD-style license found in the
% LICENSE file in the root directory of this source tree. An additional grant
% of patent rights can be found in the PATENTS file in the same directory.


% Handle angle wraparound
% After this function lat/long ranges may have multiple rows for multiple
% boxes to search. Multiple boxes occur when a lat/long range spans across
% an edge of the map where the angle wraps around.
[latranges, longranges] = latlongwrap(latrange, longrange);

% Extract lat/long values from all_stations struct array
all_stations_lat = [all_stations.lat];
all_stations_long = [all_stations.lon];

% Initialize logical array to false
i_stations = false(size(all_stations));

% Initialize latbox and longbox output vectors
latbox = [];
longbox = [];

% Loop through the wrapped lat/long ranges
for i = 1:size(latranges, 1)

    % Get idx and lat/long box for the current range
    [idx_curr, latbox_curr, longbox_curr] = evaluate_range(...
        all_stations_lat, all_stations_long, latranges(i,:), longranges(i,:));

    % Append idx_curr to i_stations
    i_stations = i_stations | idx_curr;

    % Append lat/long box
    latbox = [latbox latbox_curr]; %#ok<AGROW>
    longbox = [longbox longbox_curr]; %#ok<AGROW>
end

% Index into all_stations
stations = all_stations(i_stations);

end

function [i_stations, latbox, longbox] = evaluate_range(...
    all_stations_lat, all_stations_long, latrange, longrange)

    % Find stations in latitude range
    ilat = find_in_range(all_stations_lat, latrange);

    % Find stations in longitude range
    ilong = find_in_range(all_stations_long, longrange);

    % Find stations in both lat and long ranges
    i_stations = (ilat & ilong);

    % Process lat/long range for pretty plotting
    [latbox, longbox] = latlongbox(latrange, longrange);

end

function in_range = find_in_range(val, range)
% return logical index for values betweeen range(1) and range(2) inclusive

in_range = val >= range(1) & val <= range(2);

end

function [lat, long] = latlongwrap(lat, long)

% Error check that lat/long are vectors
if ~isvector(lat) || ~isvector(long)
    error('lat and long must be vectors')
end

% Enforce lat/long are row vectors
lat = lat(:)';
long = long(:)';

% Erroc check that lat/long are length 2
if length(lat) ~= 2 || length(long) ~= 2
    error('lat and long must be length 2')
end

% Ensure lat/long are betwee -180 and 180
lat =  dmodpi(lat);
long = dmodpi(long);

% Error check that latitude is between -90 and 90
if any(lat < -90 | lat > 90)
    error('lat must be between -90 and 90')
end

% Get booleans for whether lat/long ranges wrap around edges of map
latwrap = lat(2) <= lat(1);
longwrap = long(2) <= long(1);

% If both lat and long wrap, then we need four search boxes going to the
% edge of the map at all four courners
if latwrap && longwrap
    lat = [lat(1) 90;
        lat(1) 90;
        -90    lat(2);
        -90    lat(2)];

    long = [long(1) 180;
        -180    long(2);
        long(1) 180;
        -180    long(2)];

% If only lat is wrapped then we need two search boxes going to the top and
% bottom of the map
elseif latwrap
    lat = [lat(1) 90;
        -90    lat(2)];

    long = [long; long];

% If only lat is wrapped then we need two search boxes going to the left
% and right edges of the map
elseif longwrap
    long = [long(1) 180;
        -180    long(2)];

    lat = [lat; lat];

% Neither is wrapped
else
    % In this case lat and long are good as is, just a single search box
end


end

function ang_wrapped = dmodpi(ang_deg)
% Ensure ang_deg is between -180 and 180

ang_wrapped = mod(ang_deg + 180, 360) - 180;

end

function [lat_box, long_box] = latlongbox(latrange, longrange)
% Return vectors that can be used to plot the search box made by latrange
% and longrange but with lines on the edge of the map hidden to show how
% the box is wrapped around to the other side of the map.
%
% To achieve this each of the lines of the box are constructed one at a
% time and if the line is found to be on the edge of the map it is replaced
% with NaN

% These are the points of the search box with the first point repeated at
% the end to complete the circuit
long_pts = [longrange(1) longrange(2) longrange(2) longrange(1) longrange(1)];
lat_pts =  [latrange(1)  latrange(1)  latrange(2)  latrange(2)  latrange(1)];

% initalize outputs
long_box = [];
lat_box = [];

% Loop throught he 4 edges of the square
for i = 1:4


    % For this edge pull the next two points form long_pts and then check
    % those points to see if they lie at the extreme and NaN if so.
    curr_long = NaN_if_all_same_extreme(long_pts(i:i+1), [-180 180]);
    curr_lat =  NaN_if_all_same_extreme( lat_pts(i:i+1), [-90 90]);

    % Append
    long_box = [long_box curr_long  NaN]; %#ok<AGROW>
    lat_box =  [lat_box  curr_lat   NaN]; %#ok<AGROW>
end

end

function vals = NaN_if_all_same_extreme(vals, extremes)
% If vals are equal to each other and equal to some value in extremes then
% return NaN

L = length(vals);
if vals(1) == vals(2:L)
    if any(vals(1) == extremes)
        vals(1:L) = NaN;
    end
end

end
