function defs = format_definitions( )
% Copyright (c) Facebook, Inc. and its affiliates.
%
% This source code is licensed under the MIT license found in the
% LICENSE file in the root directory of this source tree.
%
% [ defs ] = igra_format_definitions( )
% Returns a struct containing the format definitions for igra data and
% stations list text files
%
% See the following references for more information on the IGRA text file
% format.
% Ref: https://www1.ncdc.noaa.gov/pub/data/igra/data/igra2-data-format.txt
% Ref: https://www1.ncdc.noaa.gov/pub/data/igra/igra2-list-format.txt


% Load up definitions in cell matrix format
headers_cell =  headers_format_definition_as_cell_matrix;
entries_cell =  entries_format_definition_as_cell_matrix;
stations_cell = stations_format_definition_as_cell_matrix;

% Convert definitions to structs
defs.header    = cell2struct(headers_cell);
defs.entries   = cell2struct(entries_cell);
defs.stations  = cell2struct(stations_cell);


end

function def = cell2struct(cell_matrix)

%%%% Convert cell arrays to structs
% field names for each column in the above cell matrices
columns = {'varname', 'type', 'col_idx', 'bad_vals', ...
           'function_handle', 'units', 'description'};

% Loop through the rows in the cell array
def.row_width = 0; % initialize row_width
for i_var = 1:size(cell_matrix, 1)
    curr_varname = cell_matrix{i_var, 1};

    % Loop through the columns
    for i_col = 1:length(columns)
        curr_col_name = columns{i_col};

        % Assign struct field to cell array element
        def.params.(curr_varname).(curr_col_name) = cell_matrix{i_var, i_col};
    end

    % Get row_widht by finding the maximum col_idx value
    def.row_width = max(def.row_width, def.params.(curr_varname).col_idx(end));
end

end

function out = scale_func(x, scale_factor)
% Helper function for applying a scale_factor and converting to single
% This is defined her for use in the function_handle column of the
% definitions below

    out = single(x)*scale_factor;

end

function def = headers_format_definition_as_cell_matrix()
% Source: https://www1.ncdc.noaa.gov/pub/data/igra/data/igra2-data-format.txt
% ---------------------------------
% Variable   Columns  Type
% ---------------------------------
% HEADREC       1-  1  Character
% ID            2- 12  Character
% YEAR         14- 17  Integer
% MONTH        19- 20  Integer
% DAY          22- 23  Integer
% HOUR         25- 26  Integer
% RELTIME      28- 31  Integer
% NUMLEV       33- 36  Integer
% P_SRC        38- 45  Character
% NP_SRC       47- 54  Character
% LAT          56- 62  Integer
% LON          64- 71  Integer
% ---------------------------------
def = {
    'id',         'char',   2:12,  {},     [], '', 'Station ID';
    'year',       'uint16', 14:17, {},     [], '', 'Year';
    'month',      'uint8',  19:20, {},     [], '', 'Month';
    'day',        'uint8',  22:23, {},     [], '', 'Day';
    'hour',       'uint8',  25:26, {},     [], '', 'Hour';
    'reltime_hr', 'uint8',  28:29, {'99'}, [], '', 'Release Time Hour';
    'reltime_min','uint8',  30:31, {'99'}, [], '', 'Release Time Minute';
    'numlevs',    'uint32', 33:36, {},     [], '', '# of levels in the sounding';
    'p_src',      'char',   38:45, {''},   [], '', 'Data Source Code for Pressure Levels';
    'np_src',     'char'    47:54, {''},   [], '', 'Data Source Code for Non-pressure Levels';
    'lat',        'int32',  56:62, {},     @(x)scale_func(x,1e-4), 'deg', 'Latittude';
    'lon',        'int32',  64:71, {},     @(x)scale_func(x,1e-4), 'deg', 'Longitude';
    };
end

function def = entries_format_definition_as_cell_matrix()
% Source: https://www1.ncdc.noaa.gov/pub/data/igra/data/igra2-data-format.txt
% -------------------------------
% Variable        Columns Type
% -------------------------------
% LVLTYP1         1-  1   Integer
% LVLTYP2         2-  2   Integer
% ETIME           4-  8   Integer
% PRESS          10- 15   Integer
% PFLAG          16- 16   Character
% GPH            17- 21   Integer
% ZFLAG          22- 22   Character
% TEMP           23- 27   Integer
% TFLAG          28- 28   Character
% RH             29- 33   Integer
% DPDP           35- 39   Integer
% WDIR           41- 45   Integer
% WSPD           47- 51   Integer
% -------------------------------
% defs.entries = {
%     'lvltyp1',   'uint8',  1,     {},                 [], '';
%     'lvltyp2',   'uint8',  2,     {},                 [], '';
%     'etime',     'int32',  4:8,   {'-8888', '-9999'}, [], '';
%     'press',     'uint32', 10:15, {'-8888', '-9999'}, [], 'Pa';
%     'pflag',     'char',   16,    {''},               [], '';
%     'gph',       'uint16', 17:21, {'-8888', '-9999'}, @(x)scale_func(x, 1e-3), 'km';
%     'zflag',     'char',   22,    {''},               [], '';
%     'temp',      'uint16', 23:27, {'-8888', '-9999'}, @(x)scale_func(x, 0.1), 'deg C';
%     'tflag',     'char',   28,    {''},               [], '';
%     'rh',        'uint16', 29:33, {'-8888', '-9999'}, @(x)scale_func(x, 0.1), '%';
%     'dpdp',      'uint16', 35:39, {'-8888', '-9999'}, @(x)scale_func(x, 0.1), 'deg C'
%     'wdir',      'uint16', 41:45, {'-8888', '-9999'}, [], 'deg from North';
%     'wspd',      'uint16', 47:51, {'-8888', '-9999'}, @(x)scale_func(x,0.1), 'm/s';
%     };

def = {
    'lvltyp1',   'uint8',  1,     {},                 [], '', '';
    'lvltyp2',   'uint8',  2,     {},                 [], '', '';
    'etime_min', 'int32',  4:6,   {'-88', '-99'},     [], '', '';
    'etime_sec', 'int32',  7:8,   { '88',  '99'},     [], '', '';
    'press',     'int32',  10:15, {'-8888', '-9999'}, [], 'PA', 'Pressure';
    'pflag',     'char',   16,    {''},               [], '', 'Pressure Flag';
    'gph',       'int32',  17:21, {'-8888', '-9999'}, @(x)scale_func(x, 1e-3), 'km', 'Altitude';
    'zflag',     'char',   22,    {''},               [], '', 'Altitude Flag';
    'temp',      'int16',  23:27, {'-8888', '-9999'}, @(x)scale_func(x, 0.1), 'deg C', 'Temperature';
    'tflag',     'char',   28,    {''},               [], '', 'Temperature Flag';
    'rh',        'int16',  29:33, {'-8888', '-9999'}, @(x)scale_func(x, 0.1), '%', 'Relative Humidity';
    'dpdp',      'int16',  35:39, {'-8888', '-9999'}, @(x)scale_func(x, 0.1), 'deg C', 'Dewpoint Depresesion';
    'wdir',      'int16',  41:45, {'-8888', '-9999'}, [], 'deg from North', 'Wind direction (90 = East)';
    'wspd',      'int16',  47:51, {'-8888', '-9999'}, @(x)scale_func(x,0.1), 'm/s', 'Wind Speed';
    };
end

function def = stations_format_definition_as_cell_matrix()

% Source: https://www1.ncdc.noaa.gov/pub/data/igra/igra2-list-format.txt
% ------------------------------
% Variable   Columns   Type
% ------------------------------
% ID            1-11   Character
% LATITUDE     13-20   Real
% LONGITUDE    22-30   Real
% ELEVATION    32-37   Real
% STATE        39-40   Character
% NAME         42-71   Character
% FSTYEAR      73-76   Integer
% LSTYEAR      78-81   Integer
% NOBS         83-88   Integer
% ------------------------------
def = {
    'id',   'char',    1:11,   {},       [], '', 'Station ID';
    'lat',  'single',  13:20,  {'-98.8888'},  [], '', 'Latitude';
    'long', 'single',  22:30,  {'-998.8888'}, [], '', 'Longitude';
    'elev', 'single',  32:37,  {'-998.8'},    [], '', 'Elevation';
    'state',     'char',    39:40,  {}, [], '', 'U.S. State';
    'name',      'char',    42:71,  {}, [], '', 'Name';
    'fstyear',   'int16',   73:76,  {}, [], '', 'First Year';
    'lstyear',   'int16',   78:81,  {}, [], '', 'Last Year';
    'nobs',      'int32',   83:88,  {}, [], '', '# of Observations';
    };
end
