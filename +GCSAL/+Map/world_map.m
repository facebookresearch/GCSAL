function world_map(countries)
% Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved
%
% world_map(countries)
%   Creates a cartesion map of the world by creating patches for each
%   country in countries. countries should be a struct array with each
%   element containing the fields Lat and Lon for the latitude/longitude of
%   the country borders in degrees.


% Background patch for ocean in light blue
patch([-180 -180 180 180], [-90 90 90 -90], [0 1 1])

% Loop through all countries and make a yellow patch with a grey border
% multipatch handles the fact that the country borders may have NaN for
% separating non-continuous borders
for i = 1:length(countries)
    GCSAL.Map.multipatch(countries(i).Lon, countries(i).Lat, [1 1 0], ...
        'EdgeColor', [.7 .7 .7]);
end

% axis equal so map does not distort
axis equal

% Limit axes by longitude and latitude min/max values
axis([-180 180 -90 90])

end
