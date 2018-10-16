classdef GCSAL < handle
    % Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved
    %
    % GCSAL - Global Climate Statical Analysis Library
    %   Inherits handle class making objects of this class pointers


    properties
        h5_fname % path to .h5 source file for loading data
        h5_info  % struct returned by h5info()
        headers  % struct array, each element containing header data from an IGRA station
        stations % struct array, each element containing id, latitude, and longitude information for an IGRA station
        defs     % parameter format definitions
        countries  % struct array, each element countaining infomration for a country including its name, lat/long of its borders, and which staitons are contained by that country
        entries % struct of entries data that has been cached
        do_cache_entries  % boolean whether to cache entries are not. Normally true but if you are running out of RAM you can turn this off
        quiet_mode  % Suppress status text messages and waitbars
        plot_mode  % Suppress plots
    end

    methods
        function obj = GCSAL(in_file)
            % obj = GCSAL(in_file)
            %   Create a GCSAL object from in_file. in_file should be the
            %   path to either a .h5 or .mat file. The .h5 file could have
            %   been created with GCSAL.IGRA.mat2h5_dir for example. A .mat
            %   file would have been created by an earlier call to this
            %   constructor function

            % Default to output all comments
            obj.quiet_mode = false;
            
            % Default to make all plots
            obj.plot_mode = true;
            
            % Set do_cache_entries true
            obj.do_cache_entries = true;

            % Load format definitions
            obj.defs = GCSAL.IGRA.format_definitions();

            % Expected fields for loading from/saving to .mat file
            flds = {'h5_info', 'h5_fname', 'headers', 'countries', 'stations'};

            % Initialize obj.entries to empty struct
            obj.clear_entries_cache();

            % Extract extension from in_file
            [~,~, ext] = fileparts(in_file);

            % Switch on file extension of in_file
            switch ext

                % For .h5 file, load basic info from .h5 and save to .mat
                case '.h5'

                    % in_file was .h5 so assign h5_fname
                    obj.h5_fname = in_file;

                    % Get h5info
                    tic; fprintf('Parsing info from h5 file...');
                    obj.h5_info = h5info(obj.h5_fname);
                    fprintf(' Complete in %.1f seconds\n', toc);

                    % Load all headers to RAM for quicker data access
                    tic; fprintf('Loading headers...');
                    obj.headers = obj.load_all_headers();
                    fprintf(' Complete in %.1f seconds\n', toc);

                    % Create stations struct array from headers and create
                    % countries from stations struct. These data
                    % structures offer quick access to finding structures
                    % based on their lat/long location or which country
                    % they are in
                    tic; fprintf('Initiazling stations struct and country map ...');
                    obj.stations = GCSAL.GCSAL.initialize_stations(obj.headers);
                    obj.countries = GCSAL.Map.map_stations_by_country(obj.stations);
                    fprintf(' Complete in %.1f seconds\n', toc);

                    % Save the workspace to a .mat file for quicker loading
                    % in the future
                    tic; fprintf('Saving h5 info to .mat file for faster loading...');
                    for i = 1:length(flds)
                        to_mat.(flds{i}) = obj.(flds{i}); %#ok<STRNU>
                    end
                    save([in_file '.info.mat'], '-struct', 'to_mat');
                    fprintf(' Complete in %.1f seconds\n', toc);


                    % For .mat file, just load and check expected variables are
                    % present
                case '.mat'

                    % Load .mat file
                    tic; fprintf('Loading h5 info from .mat file...\n');
                    mat_data = load(in_file);
                    fprintf(' Complete in %.1f seconds\n', toc);

                    % Error check .mat file had all necessary flds
                    found_flds = fieldnames(mat_data);
                    flds_not_found_idx = ~ismember(flds, found_flds);
                    if ~all(ismember(flds, found_flds))
                        msg = sprintf('  %s\n', flds{flds_not_found_idx});
                        error('mat file is missing the following fields: %s', msg)
                    end

                    % Assign data in mat_data to obj
                    for i = 1:length(found_flds)
                        obj.(found_flds{i}) = mat_data.(found_flds{i});
                    end

                otherwise
                    error('Unexpected file extension')
            end
        end

        function [N, entries] = counts(obj, stations, fld, varargin)
            % [N, entries] = counts(obj, stations, fld, varargin)
            %    Returns counts from  histcounts and makes
            %    plots showing histogram, probability and cumulative
            %    density functions for data in fld at stations
            %
            %    Additional optional parameters can be given as Name/Value
            %    pairs. See examples below.
            %
            %    'Edges' defines the edges of the bins.
            %    'FilterFields' and 'FilterRanges' together can be used to
            %    filter the data to a subset where the parameter in
            %    FilterFields matches the range in FilterRanges.
            %
            %
            % Example:
            %    % Create GCSAL object
            %    g = GCSAL.GCSAL('gcsal.h5.info.mat');
            %
            %    % Choose some stations. In this case stations within 2
            %    % degrees of the equator
            %    stations = g.station_search('LatLong', [-2 2 -180 180]);
            %
            %    % Histogram for all windspeeds
            %    [N, entries] = g.counts(stations, 'wspd');
            %
            %    % Define custom bin_edges
            %    bin_edges = 0:1:80;
            %
            %    % Do counts for windspeed filtered on geopotential
            %    % altitude between 20 and 30 km and with custom bin edges
            %    [N, entries] = g.counts(stations, 'wspd', 'FilterFields', {'gph'}, ...
            %    'FilterRanges', {[20 30]}, 'Edges', bin_edges);
            %
            %    % Now additionally filter on measurements taken in August
            %    % between 4 and 10 am
            %    [N, entries] = g.counts(stations, 'wspd', 'FilterFields', {'gph', 'month', 'hour'}, ...
            %     'FilterRanges', {[0 8], 8, [4 10]}, 'Edges', bin_edges);

            % Return on empty stations
            if isempty(stations)
                warning('Stations is empty, cannot count')
                N = []; entries = [];
                return
            end

            % Parse Name/Value pairs from input
            p = inputParser;
            addOptional(p, 'Edges', GCSAL.GCSAL.default_bin_edges(fld));
            addOptional(p, 'FilterFields', {});
            addOptional(p, 'FilterRanges', {});
            addOptional(p, 'Plot', []);
            addOptional(p, 'Verbose', []);
            parse(p, varargin{:});

            % Rename for convenience
            edges = p.Results.Edges;
            fltr_flds = p.Results.FilterFields;
            fltr_ranges = p.Results.FilterRanges;
            obj.plot_mode = (~ismember('Plot', p.UsingDefaults) && ...
                             (p.Results.Plot == true)) || ...
                             ismember('Plot', p.UsingDefaults);
            obj.quiet_mode = (~ismember('Verbose', p.UsingDefaults) && ...
                             (p.Results.Verbose == false));

            % Load data from stations and all fields
            entries = obj.query(stations, fld, fltr_flds, fltr_ranges);

            % Get counts
            [N] = GCSAL.GCSAL.histcounts(entries, fld, edges, obj.quiet_mode);
            [pdf, cdf] = GCSAL.GCSAL.counts2pdf(N, 2);

            % calculate bin centers from edges
            centers = GCSAL.GCSAL.get_bin_centers(edges);

            % Extract parameter definitions
            def = obj.find_def(fld);

            % Construct labels from definitions
            label = GCSAL.GCSAL.get_label(def);

            % Get string describing current filters in place for use in
            % title
            title_str = GCSAL.GCSAL.description_from_filters(fltr_flds, fltr_ranges);

            if (obj.plot_mode)
                % Make figure
                figure;

                % Plot histogram
                subplot(3,1,1)
                histogram(vertcat(entries.(fld)), edges)
                title(sprintf('Histogram\n%s', title_str))
                xlabel(label)
                ylabel('# of occurences')
                
                % Plot probability density function
                subplot(3,1,2)
                plot(centers, pdf, '-x')
                title(sprintf('Probability Density Function\n%s', title_str))
                
                xlabel(label)
                ylabel('Probability of occuring')
                
                % Plot cumulative density funciton
                subplot(3,1,3)
                plot(centers, cdf, '-x')
                title(sprintf('Cumulative Density Function\n%s', title_str))
                
                xlabel(label)
                ylabel('Probability of exceeding')
            end
            
            obj.plot_mode = true;
            obj.quiet_mode = false;
        end

        function [N, entries, stats] = counts2(obj, stations, x_fld, y_fld, varargin)
            % [N, entries] = counts2(obj, stations, x_fld, y_fld, varargin)
            %    Returns two dimensional counts from  histcounts2 and makes
            %    plots showing two dimensional probability and cumulative
            %    density functions comparing data in x_fld and y_fld of
            %    stations.
            %
            %    Additional optional parameters can be given as Name/Value
            %    pairs. See examples below.
            %
            %    'XEdges', 'YEdges' defines the edges of the bins.
            %    'FilterFields' and 'FilterRanges' together can be used to
            %    filter the data to a subset where the parameter in
            %    FilterFields matches the range in FilterRanges.
            %
            % Example:
            %    % Create GCSAL object
            %    g = GCSAL.GCSAL('gcsal.h5.info.mat');
            %
            %    % Choose some stations. In this case stations within 2
            %    % degrees of the equator
            %    stations = g.station_search('LatLong', [-2 2 -180 180]);
            %
            %    % Do counts between gph and wspd
            %    [N, entries] = g.counts2(stations, 'gph', 'wspd');
            %
            %    % Do counts between gph and pressure with custom bin
            %    % edges
            %    [N, entries] = g.counts2(stations, 'gph', 'press', ...
            %          'XEdges', 0:0.5:40, 'YEdges', 0:2000:100000);
            %
            %    % Do counts for data measured between 6 and 10 am in August
            %    [N, entries] = g.counts2(stations, 'gph', 'wspd', ...
            %           'FilterFields', {'hour', 'month'}, ...
            %           'FilterRanges', {[6 10], [8 8]});

          % Return on empty stations
            if isempty(stations)
                warning('Stations is empty, cannot count')
                N = []; entries = [];
                return
            end

            % Parse Name/Value pairs from input
            p = inputParser;
            addOptional(p, 'XEdges', GCSAL.GCSAL.default_bin_edges(x_fld));
            addOptional(p, 'YEdges', GCSAL.GCSAL.default_bin_edges(y_fld));
            addOptional(p, 'FilterFields', {});
            addOptional(p, 'FilterRanges', {});
            addOptional(p, 'Plot', []);
            addOptional(p, 'Verbose', []);
            parse(p, varargin{:});

            % Rename for convenience
            x_edges = p.Results.XEdges;
            y_edges = p.Results.YEdges;
            fltr_flds = p.Results.FilterFields;
            fltr_ranges = p.Results.FilterRanges;
            obj.plot_mode = (~ismember('Plot', p.UsingDefaults) && ...
                             (p.Results.Plot == true)) || ...
                             ismember('Plot', p.UsingDefaults);
            obj.quiet_mode = (~ismember('Verbose', p.UsingDefaults) && ...
                             (p.Results.Verbose == false));

            % Load data from stations and all fields
            entries = obj.query(stations, {x_fld, y_fld}, fltr_flds, fltr_ranges);

            % Get counts
            [N] = GCSAL.GCSAL.histcounts2(entries, x_fld, y_fld, ...
                x_edges, y_edges, obj.quiet_mode);
            [pdf, cdf] = GCSAL.GCSAL.counts2pdf(N, 2);
            pdf = pdf';
            cdf = cdf';

            % calculate bin centers from edges
            x_centers = GCSAL.GCSAL.get_bin_centers(x_edges);
            y_centers = GCSAL.GCSAL.get_bin_centers(y_edges);

            % Extract parameter definitions
            x_def = obj.find_def(x_fld);
            y_def = obj.find_def(y_fld);

            % Construct labels from definitions
            x_label = GCSAL.GCSAL.get_label(x_def);
            y_label = GCSAL.GCSAL.get_label(y_def);

            stats.x = x_centers;
            stats.y = y_centers;
            stats.cdf = cdf;
            stats.pdf = pdf;
            
            if (obj.plot_mode)
                % Make figure
                figure;
                
                % subplot for contour of cumulative density function
                subplot(3,1,1)
                [C, h] = contourf(x_centers, y_centers, cdf, ...
                    [0.05:0.05:0.95 0.99]);
                clabel(C,h, 0.1:0.2:0.9, 'LabelSpacing', 600, 'FontSize', 18);
                xlabel(x_label)
                ylabel(y_label)
                title('Percentile by Altitude')
                
                % subplot for surf of probability density function
                subplot(3,1,2)
                surf(x_centers, y_centers, pdf, 'LineStyle', 'None');
                view([0 0 1])
                xlabel(x_label)
                ylabel(y_label)
                title('Probability Density Function by Altitude')
                h_colorbar = colorbar;
                ylabel(h_colorbar, 'probability')
                
                % subplot for # of samples
                subplot(3,1,3)
                plot(x_centers, sum(N,2)/1000, '-x');
                xlabel(x_label)
                ylabel('Thousands of Counts')
                title('Sample size')
            end
            
            obj.plot_mode = true;
            obj.quiet_mode = false;
        end

        function [N, entries] = countsN(obj, stations, resolutions, varargin)
            % [N, entries] = countsN(obj, stations, resolutions, varargin)
            %    Returns N dimensional counts similar to histcounts but in
            %    N dimensions. Counts are performed based on the fields in
            %    the struct array resolutions. Each struct in resolutions
            %    must contain a variable name in the field 'fld' and
            %    optionally may have bin edges specified in the field
            %    'edges'.
            %
            %    Additional optional parameters can be given as Name/Value
            %    pairs. See examples below.
            %
            %    'FilterFields' and 'FilterRanges' together can be used to
            %    filter the data to a subset where the parameter in
            %    FilterFields matches the range in FilterRanges.
            %
            % Examples:
            %   % Get some stations
            %   stations = g.station_search('Countries', 'Brazil');
            %
            %   % Make 5-dimensional count matrix with default bin edges and no filtering
            %   resolutions = struct('fld', {'lat', 'lon', 'gph', 'month', 'wspd'});
            %   N = g.countsN(stations, resolutions);
            %
            %   % Add custom bin edges to gph field and limit data to only data between 6
            %   % and 10 am
            %   resolutions(3).edges = 0:1:80;
            %   N = g.countsN(stations, resolutions, 'FilterFields', {'hour'}, ...
            %                 'FilterRanges', [6 10]);


          % Return on empty stations
            if isempty(stations)
                warning('Stations is empty, cannot count')
                N = []; entries = [];
                return
            end

            % Parse Name/Value pairs from input
            p = inputParser;
            addOptional(p, 'FilterFields', {});
            addOptional(p, 'FilterRanges', {});
            addOptional(p, 'Plot', []);
            addOptional(p, 'Verbose', []);
            parse(p, varargin{:});

            % Rename for convenience
            fltr_flds = p.Results.FilterFields;
            fltr_ranges = p.Results.FilterRanges;
            obj.plot_mode = (~ismember('Plot', p.UsingDefaults) && ...
                             (p.Results.Plot == true)) || ...
                             ismember('Plot', p.UsingDefaults);
            obj.quiet_mode = (~ismember('Verbose', p.UsingDefaults) && ...
                             (p.Results.Verbose == false));

            % fill in edges
            for i = 1:length(resolutions)
                if ~isfield(resolutions(i), 'edges') || isempty(resolutions(i).edges)
                    resolutions(i).edges = GCSAL.GCSAL.default_bin_edges(resolutions(i).fld);
                end
            end

            flds = {resolutions.fld};
%             edges = {resolutions.edges};

            % Load data from stations and all fields
            entries = obj.query(stations, flds, fltr_flds, fltr_ranges);

            % Get counts
            [N] = GCSAL.GCSAL.histcountsN(entries, resolutions, obj.quiet_mode);

            obj.plot_mode = true;
            obj.quiet_mode = false;
        end

        function entries = query(obj, stations, params, fltr_flds, fltr_rngs)
            % entries = query(obj, station_ids, params, fltr_flds, fltr_rngs)
            %   Returns a struct array with each element containing the
            %   data for params for a station in stations.
            %
            %   Finds the data either by reading from the H5 file located
            %   at obj.h5_fname or by findind the data cached in
            %   obj.entries or obj.headers.
            %
            % INPUTS
            %     stations - A string, cell array of strings of struct array
            %                containing the station ids from which to load
            %                data
            %       params - params can be a string or cell array of strings.
            %                params can be either from the entry data or
            %                header data but at least one string in params
            %                must be from entry data
            %    fltr_flds - (optional) cell array of filtering parameter names
            %    fltr_rngs - (optional) cell array of filtering ranges
            %
            %   For a list of available by params see g.defs.header.params
            %   and g.defs.entries.params
            %
            %   If params or staiton_ids is empty, returns an empty struct
            %
            %   Examples:
            %      % Create GCSAL object
            %      g = GCSAL.GCSAL('gcsal.h5.info.mat');
            %
            %      % Get stations located in Botswana
            %      stations = g.station_search('Countries', {'Botswana'});
            %
            %      % Get all geopotential height and windspeed data as well
            %      % as hour, month, and year data
            %      entries1 = g.query(stations, {'gph', 'wspd', 'hour', 'month', 'year'});
            %
            %      % Plot distribution of hours and years for the data in entries1
            %      figure; histogram(vertcat(entries1.hour))
            %      figure; histogram(vertcat(entries1.year))
            %
            %      % Get gph and wspd data measured between 6 and 4 pm
            %      entries2 = g.query(stations, {'gph', 'wspd'}, 'hour', [6 16]);
            %
            %      % Plot distribution of hours for the data in entries2
            %      figure; histogram(vertcat(entries2.hour))
            %
            %      % Get data corresponding only to measuresments taken
            %      % in August between 4am and Noon and in the years 1990
            %      % to 1999
            %      entries3 = g.query(stations, {'gph', 'wspd'}, ...
            %                         {'month', 'hour', 'year'}, ...
            %                         {8, [4 12], [1990 1999]});
            %
            %      % Plot distribution of years for the data in entries3
            %      figure; histogram(vertcat(entries3.year))


            % Return empty struct array if no stations ids or params given
            if isempty(stations) || isempty(params)
                entries = struct([]);
                return
            end

            % Set filter parameters to empty cell arrays if not given
            if ~exist('fltr_flds', 'var')
                fltr_flds = {};
            end

            if ~exist('fltr_rngs', 'var')
                fltr_rngs = {};
            end

            % Ensure filter parameters are cell arrays
            fltr_flds = cellstr(fltr_flds);
            if ~iscell(fltr_rngs)
                fltr_rngs = {fltr_rngs};
            end

            % Error check varargin was input in pairs
            if length(fltr_flds) ~= length(fltr_rngs)
                error('Expected length of filter_flds and filter_ranges to match')
            end

            % Add all restriction range fields to params since they need to
            % be loaded as well
            params = unique([params fltr_flds]);

            % Get all header and entries parameters
            all_header_params = fieldnames(obj.defs.header.params);
            all_entries_params = fieldnames(obj.defs.entries.params);

            % Verify all params can be  found in either header or entries
            found = ismember(params, [all_header_params; all_entries_params]);
            if any(~found)
                msg = sprintf('  %s\n', params{~found});
                error('The following params were invalid: %s\n', msg)
            end

            % Determine whether each param is part of the header data or
            % entry data
            curr_header_params = intersect(params, all_header_params);
            curr_entries_params = intersect(params, all_entries_params);

            % Make sure at least one param is in entries group
            if isempty(curr_entries_params)
                error('At least one param must be an entry')
            end

            %Ensure stations is a character array
            station_ids = GCSAL.GCSAL.station_id_str(stations);

            % Load data from entries
            entries = obj.load_from_stations('entries', station_ids, curr_entries_params);

            % Add data from headers
            entries = obj.add_header_params_to_entries(entries, curr_header_params);

            % Filter data according to range limits
            if (~obj.quiet_mode)
                tic;
                fprintf('Applying filters... ');
            end
            for i = 1:length(fltr_rngs)
                entries = GCSAL.GCSAL.filter_data_by_range(entries, fltr_flds{i}, fltr_rngs{i});
            end
            if (~obj.quiet_mode)
                fprintf('Complete in %.1f seconds\n', toc);
            end

            % Clear entries if nargout is 0 so we don't use any RAM in the
            % base workspace for "ans"
            if nargout == 0
                entries = [];
            end
        end

        function stations_match = station_search(obj, varargin)
            % stations_match = station_search(obj, varargin)
            %   Search for stations by Latitude and Longitude, or by
            %   Country, or by name. Search criteria are given in
            %   Name/Value pairs with the key words 'LatLong, 'Countries',
            %   and "IDRegex". If multiple criteria are given they are
            %   combined with AND.
            %
            % Examples:
            %    % Create GCSAL object
            %    g = GCSAL.GCSAL('gcsal.h5.info.mat');
            %
            %   % Find stations within 25 degrees of the equator
            %   stations1 = g.station_search('LatLong', [-25 25 -180 180]);
            %
            %   % Find stations in Brazil or India and within 25 degrees of
            %   % the equator
            %   stations2 = g.station_search('Countries', {'Brazil', 'India'}, ...
            %                               'LatLong', [-25 25 -180 180]);
            %
            %   % Find stations with IDs beginnign with the letter A.
            %   % Note that in regex ^ means beginning of the line
            %   stations3 = g.station_search('IDRegex', '^/A');
            %
            %   % Find stations in Brazil or India AND within 25 degrees of
            %   % the equator AND with station IDs ending in 5
            %   % Note that in regex $ means end of the line
            %   stations4 = g.station_search('Countries', {'Brazil', 'India'}, ...
            %                               'LatLong', [-25 25 -180 180], ...
            %                               'IDRegex', '5$');

            % Parse varargin
            p = inputParser;
            addOptional(p, 'Countries', []);
            addOptional(p, 'IDRegex', []);
            addOptional(p, 'Lat', []);
            addOptional(p, 'LatLong', []);
            addOptional(p, 'Nearest', []);
            addOptional(p, 'Number', []);
            addOptional(p, 'Range', []);
            addOptional(p, 'Plot', []);
            addOptional(p, 'Verbose', []);
            parse(p, varargin{:});

            obj.plot_mode = (~ismember('Plot', p.UsingDefaults) && ...
                             (p.Results.Plot == true)) || ...
                             ismember('Plot', p.UsingDefaults);
            obj.quiet_mode = (~ismember('Verbose', p.UsingDefaults) && ...
                             (p.Results.Verbose == false));
                         
            % Plot the map of world with all stations marked
            figure; hold all;
            obj.plot_world_map();

            % Initialize station_ids_match to all station ids
            ids_match = GCSAL.GCSAL.station_id_str(obj.stations);
            Lmax = length(ids_match);

            % If Nearest specified:
            if ~ismember('Nearest', p.UsingDefaults)
                
                % Find stations in lat/long range
                num = uint16(p.Results.Number);
                if (isempty(num) || (num < 1))
                    num = 1;
                elseif (num > Lmax)
                    num = Lmax;
                end
                [stations_nearest, arclen] = ... 
                    obj.stations_near_latlong(p.Results.Nearest, num);
                
                % Report # stations found
                L = length(stations_nearest);
                if (~obj.quiet_mode)
                    fprintf('%d stations found near lat/long\n', L)
                end
                
                % Find intersect of stations found so far and currently
                % found stations
                curr_ids = GCSAL.GCSAL.station_id_str(stations_nearest);
                ids_match = intersect(ids_match, curr_ids, 'rows');
            end
            
            % If LatLong specified:
            if (~ismember('LatLong', p.UsingDefaults) || ...
                    ~ismember('Lat', p.UsingDefaults))

                if ~ismember('Lat', p.UsingDefaults)
                    % Find stations in lat/long range
                    range = single(p.Results.Range);
                    if (isempty(range) || (range < 1.2))
                        range = 1.2;
                    end
                    % Find stations in lat/long range
                    lat = p.Results.Lat(1);
                    lon = p.Results.Lat(2);
                    box = [(lat - range) (lat + range) lon lon];
                else
                    box = p.Results.LatLong;
                end
                
                % Find stations in lat/long range
                [stations_in_range, latbox, longbox] = ...
                    obj.stations_from_latlong(box);

                % Report # stations found
                L = length(stations_in_range);
                if (~obj.quiet_mode)
                    fprintf('%d stations found in lat/long range\n', L)
                end

                % Find intersect of stations found so far and currently
                % found stations
                curr_ids = GCSAL.GCSAL.station_id_str(stations_in_range);
                ids_match = intersect(ids_match, curr_ids, 'rows');

                % Highlight lat/long search box
                plot(longbox, latbox, 'b-', 'LineWidth', 2)
            end

            % If Countries specified:
            if ~ismember('Countries', p.UsingDefaults)

                % Find stations in countries
                [stations_in_countries, countries_match] = ...
                    obj.stations_from_countries(p.Results.Countries);

                % Report # stations found
                L = length(stations_in_countries);
                if (~obj.quiet_mode)
                    fprintf('%d stations found in countries\n', L)
                end

                % Find intersect of stations found so far and currently
                % found stations
                curr_ids = GCSAL.GCSAL.station_id_str(stations_in_countries);
                ids_match = intersect(ids_match, curr_ids, 'rows');

                % Highlight border of countries searched
                for i = 1:length(countries_match)
                    plot(countries_match(i).Lon, countries_match(i).Lat, ...
                        'b-', 'linewidth', 2)
                end
            end

            % IfIDRegex specified:
            if ~ismember('IDRegex', p.UsingDefaults)

                % Found stations matching IDRegex
                stations_from_regex = obj.stations_from_regex(p.Results.IDRegex);

                % Report # stations found
                L = length(stations_from_regex);
                if (~obj.quiet_mode)
                    fprintf('%d stations found matching search_str\n', L)
                end

                % Find intersect of stations found so far and currently
                % found stations
                curr_ids = GCSAL.GCSAL.station_id_str(stations_from_regex);
                ids_match = intersect(ids_match, curr_ids, 'rows');

                % Highlight stations matching IDRegex
                GCSAL.GCSAL.plot_stations(stations_from_regex, ...
                    'bo', 'MarkerSize', 6);

            end

            % Convert station ids to stations struct array
            stations_match = obj.find_stations(ids_match);
            
            % Report # stations found
            if (~obj.quiet_mode)
                fprintf('%d stations found combined\n', length(stations_match))
            end

            % Highlight stations found
            GCSAL.GCSAL.plot_stations(stations_match, 'r+');

            if ~ismember('Nearest', p.UsingDefaults)
                for i = 1:num
                    stations_match(i).arclen = arclen(i);
                end
            end
            
            obj.plot_mode = true;
            obj.quiet_mode = false;
        end

        function plot_world_map(obj, include_stations)
            % plot_world_map(obj, include_stations)
            %   Plots a map of the world based on the country borders in
            %   obj.countries. If include_stations is true then will also
            %   put a mark for each station in obj.stations. If
            %   include_stations is not given then it defaults to true.


            % Set default value
            if ~exist('include_stations', 'var')
                include_stations = true;
            end

            % Plot world map
            GCSAL.Map.world_map(obj.countries);

            % Plot stations
            if include_stations
                obj.plot_stations(obj.stations, 'k.');
            end
        end

        function country_matches = find_countries(obj, country_names)
            % country_matches = find_countries(obj, country_names)
            %   Returns a struct array corresponding to elements of
            %   obj.countries with a name matching any string in
            %   country_names. country_names can be a string or cell array
            %   of strings. Ignores case.

            % Ensure country_names is lower case
            country_names = lower(country_names);

            % Get all names in countries and ensure lower case
            all_countries = lower({obj.countries.name});

            % Find matches
            i_country_match = GCSAL.GCSAL.find_keys(all_countries, country_names);

            % index into countries
            country_matches = obj.countries(i_country_match);
        end

        function station_matches = find_stations(obj, station_ids)
            % station_matches = find_stations(obj, station_ids)
            %   Returns a struct array corresponding to elements of
            %   obj.stations with an id matching any string in
            %   station_ids. station_ids can be a string or cell array
            %   of strings.

            % Get all station ids in obj.stations
            all_station_ids = GCSAL.GCSAL.station_id_str(obj.stations);

            % Find matches
            i_station_match = GCSAL.GCSAL.find_keys(all_station_ids, station_ids);

            % Index into stations
            station_matches = obj.stations(i_station_match);

        end

        function [header_matches, i_header_match] = find_headers(obj, stations)
            % header_matches = find_headers(obj, ids)
            %   Returns a struct array corresponding to elements of
            %   obj.headers with an id matching any string in
            %   ids. ids can be a string or cell array
            %   of strings.

            % Get all header ids in obj.headers
            all_ids = GCSAL.GCSAL.station_id_str(obj.headers);

            station_ids = GCSAL.GCSAL.station_id_str(stations);
            % Find matches
            i_header_match = GCSAL.GCSAL.find_keys(all_ids, station_ids);

            % Index into headers
            header_matches = obj.headers(i_header_match);
        end

        function def = find_def(obj, varname)
            % def = find_def(obj, varname)
            %   Searches through all parameters in obj.defs and returns the
            %   struct whose name matches varname

            % Loop through all groups in obj.defs
            groups = fieldnames(obj.defs);
            for i = 1:length(groups)

                % Pull out the parameter names in param
                param_names = fieldnames(obj.defs.(groups{i}).params);

                % Use ismember to search for a match
                [~, idx ] = ismember(varname, param_names);

                % If a match is found return
                if idx
                    def = obj.defs.(groups{i}).params.(varname);
                    return
                end
            end

            % If we got here without hitting a return, no match was ever found
            error('Could not find defition for %s', varname)
        end

        function clear_entries_cache(obj)
            % clears the cached data in obj.entries. Do this if you are
            % running out of RAM

            obj.entries = struct();
        end

    end

    methods (Access = 'private')

        function entries = add_header_params_to_entries(obj, entries, params)
            % entries = add_header_params_to_entries(obj, entries, flds)
            %   For each entry in entries and each param in params, adds
            %   the data in header.(param) to entry.(param). The header is
            %   found based on matching station id. The data in
            %   header.(param) is expanded to match the length of data in
            %   entry based on a correspondence index.
            %
            %   If params is not given, then all params in header will be
            %   added.


            % If params is empty, then there is nothing to add, just return
            if isempty(params)
                return
            end


            % Check if params was an input
            if ~exist('params', 'var')
                % Params not given so add all params
                add_all_params = true;
            else
                % Params specific so do not add all params
                add_all_params = false;

                % Ensure params is a cell array
                params = cellstr(params);
            end

            % Find indices of headers that match station id with entries
            [~, i_header_match] = obj.find_headers(GCSAL.GCSAL.station_id_str(entries));

            % Convert logical array to indices
            i_header_match = find(i_header_match);

            % Add entry_idx to all headers that match entries
            obj.add_entry_idx_to_headers(i_header_match);  %#ok<FNDSB>

            % Determine whether to do waitbar
            L = length(entries);
            do_waitbar = (L > 1) && ~obj.quiet_mode;
            if do_waitbar
                h = waitbar(0, 'Adding header fields to entries');
            end
            
            if (~obj.quiet_mode)
                tic;
                fprintf('Adding header fields to entries... ');
            end

            % Loop through stations
            for i = 1:L

                % Get header for current station
                header = obj.find_headers(entries(i).id);

                % If add_all_params set params to all fields in header
                if add_all_params
                    params = fieldnames(header);
                end

                % Loop through params
                for j = 1:length(params)

                    % Get the current field name from params
                    fld = params{j};

                    % Try reading from cached
                    val = obj.read_from_cached_entries(header.id, fld);

                    % If not found in cache, read from header
                    if isempty(val)

                        % Get the data in header at fld
                        val = header.(fld);


                        if size(val, 1) == 1
                            % If val is a single row we need to duplicate it to the
                            % size of the entry data
                            val = repmat(val, length(header.entry_idx), 1);
                        else
                            % Apply entry_idx to val to expand val data to
                            % match length of entry data with correct
                            % correspondence
                            val = val(header.entry_idx);
                        end

                        % Set entry data to val
                        obj.cache_param(header.id, fld, val);

                    end

                    % add data to entries for return struct array
                    entries(i).(fld) = val;

                end

                % Update waitbar
                if do_waitbar && mod(i, ceil(L/50)) == 0
                    msg = sprintf('%d/%d: Adding header fields to entries for %s', i, L, header.id);
                    waitbar(i/L, h, msg);
                end
            end
            
            if (~obj.quiet_mode)
                fprintf('Complete in %.1f seconds\n', toc);
            end

            % Close waitbar
            if do_waitbar
                close(h);
            end

        end

        function [stations_nearest, arclen] = ...
                stations_near_latlong(obj, latlonposn, n)
            % [stations_nearest, arclen] = ...
            %      stations_near_latlong(obj, latlongrange)
            %   Returns an array of n-station structs for stations that
            %   are nearest the queried position.
            %   Additionally returns an array of arclength distances in meters
            %
            %   latlonposn must be a 2 element vector and is
            %   in degrees. Example:
            %          latlonposn = [9.999924 -84.205753]
            %
            
            % Return on empty input
            if isempty(latlonposn)
                stations_nearest = struct();
                return
            end
            
            % error check lat/long range
            if length(latlonposn) ~= 2
                error('Lat Lon Position must be a 2 element vector')
            end
            
            % Find stations in range as well as getting lat/long vectors
            % for plotting the search box
            [stations_nearest, arclen] = GCSAL.Map.find_nearest(...
                obj.stations, latlonposn(1), latlonposn(2), n);
            
        end
        
        function [stations_in_range, latbox, longbox] = ...
                stations_from_latlong(obj, latlongrange)
            % [stations_in_range, latbox, longbox] = ...
            %      stations_from_latlong(obj, latlongrange)
            %   Returns an array of station structs for stations that are
            %   located within the box defined by latlongrange.
            %
            %   Additionally returns longbox and latbox which can be used
            %   to plot the the searchbox that was used.
            %
            %   latlongrange must be a four element vector and is
            %   in degrees. Example:
            %          latlongrange = [-45 45 -180 180]
            %   would find all stations between -45 and 45 deg latitude
            %
            %   latlongrange does account for angle wrap around. Example:
            %          latlongrange = [45 -45 -180 180]
            %   would finda all stations with latitude above 45 deg or
            %   below -45.

            % Return on empty input
            if isempty(latlongrange)
                stations_in_range = struct();
                return
            end

            % error check lat/long range
            if length(latlongrange) ~= 4
                error('latlongrange must be a 4 element vector')
            end

            % Find stations in range as well as getting lat/long vectors
            % for plotting the search box
            [stations_in_range, latbox, longbox] = GCSAL.Map.find_in_lat_long_range(...
                obj.stations, latlongrange(1:2), latlongrange(3:4));

        end

        function [stations_in_countries, countries_match] = ...
                stations_from_countries(obj, country_names)
            % [stations_in_countries, countries_match] = stations_from_countries(obj, country_names)
            %   Returns an array of station structs for stations that
            %   are located within the countries listed in country_names.
            %   Additionally returns a struct array for countries that
            %   match country_names

            % Get struct array of countries from countries matching
            % country_names
            countries_match = obj.find_countries(country_names);

            % Get list of all stations in matching countries
            station_ids = vertcat(countries_match.stations);

            % In case where no stations or countries were found ensure
            % station_ids is an empyt string
            if isempty(station_ids); station_ids = ''; end

            % Convert station ids to stations struct array
            stations_in_countries = obj.find_stations(station_ids);

        end

        function station_matches = stations_from_regex(obj, search_str)
            % station_matches = find_stations_regex(obj, search_str)
            %   Returns an array of station structs for stations whose ids
            %   match the regex pattern in search_str

            % Get all station ids in obj.stations
            all_station_ids = GCSAL.GCSAL.station_id_str(obj.stations);

            % Convert to cell array
            all_station_ids = cellstr(all_station_ids);

            % Call regexp
            regex_out = regexp(all_station_ids, search_str);

            %  Use cellfun to find which elements in all_station_ids had a
            %  match
            i_station_match = ~cellfun(@isempty, regex_out);

            % Index into stations
            station_matches = obj.stations(i_station_match);

        end

        function headers = load_all_headers(obj)
            % Find all headers in the h5_info struct and load the data from
            % the h5 file for those headers

            % Find all station names in h5 info struct
            all_station_ids = {obj.h5_info.Groups.Name};

            % Remove / from beginning of station ids
            all_station_ids(:,1) = [];

            % Use empty params to indicate we want to load all parameters
            params = {};
            headers = obj.load_from_stations('header', all_station_ids, params);
        end




        function out = load_from_stations(obj, group, station_ids, params)
            % Load the parameters listed in params from the data in group from
            % the H5 file for all stations in station_ids.
            %
            % If params is empty then all parameters will be loaded

            % Set default params to empty cell which will revert to loading
            % all parameters
            if ~exist('params', 'var')
                params = {};
            end

            % Handle case where station_ids are empty
            if isempty(station_ids)
                out = [];
                return
            end

            % Ensure station_ids is a cell array
            station_ids = cellstr(station_ids);

            % Initialize counter
            count = 1;

            % Decide whether to do wait bar
            L = length(station_ids);
            do_waitbar = (L > 1) && ~obj.quiet_mode;

            % Open waitbar
            if do_waitbar
                h = waitbar(0, 'Loading data from stations');
            end
            
            if (~obj.quiet_mode)
                tic;
                fprintf('Loading data from stations... ');
            end

            % Loop through all station ids
            for i = 1:L

                % Attempt to read group data
                tmp = obj.load_group(group, station_ids{i}, params);

                % Assign data to out struct if tmp is not empty
                if ~isempty(tmp)
                    % Assign data
                    out(count) = tmp; %#ok<AGROW>

                    % Increment counter
                    count = count + 1;
                end

                % Update waitbar
                if do_waitbar && mod(i, ceil(L/50)) == 0
                    msg = sprintf('%d/%d: Loading data for %s/%s', i, L, station_ids{i}, group);
                    waitbar(i/L, h, msg);
                end
            end
            if (~obj.quiet_mode)
                fprintf('Complete in %.1f seconds\n', toc);
            end

            % Close wait bar
            if do_waitbar
                close(h)
            end

            % If counter never incremented, return empty struct array
            if count == 1
                out = struct([]);
            end
        end

        function out = load_group(obj, group, station_id, params)
            % Load the parameters listed in params from the data in group
            % and from the station in station_id from the H5 file
            %
            % If params is empty then all parameters are loaded


            % Extract parameter definitions for the current group
            param_defs = obj.defs.(group);

            % Set default params to all parameters
            if ~exist('params', 'var') || isempty(params)
                params = fieldnames(param_defs.params);
            end

            % Initialize output
            out = [];

            % Find info for the current station_id in the top level h5_info
            station_info = GCSAL.GCSAL.h5info_find_children(obj.h5_info, station_id);

            % If station_info is empty return with warning that station_id
            % was not found
            if isempty(station_info)
                fprintf('%s not found\n', station_id)
                return
            end

            % Find group in station_info
            group_info = GCSAL.GCSAL.h5info_find_children(station_info, group);

            % Throw error if group not found
            if isempty(group_info)
                error('Group not found: %s', group)
            end

            % For each parameter in params, load the param
            for i = 1:length(params)
                curr = param_defs.params.(params{i});
                out.(curr.varname) = obj.load_param(curr, group_info);
            end

            % Add id to struct so that the header data for this struct can
            % be easily found
            if ~isfield(out, 'id')
                out.id = station_id;
            end
        end


        function data = load_param(obj, param_def, group_info)
            % Load the data corresponding to param_def and group_info.
            % Apply data conversion and function_handle as
            % specified in param_def
            %
            % INPUTS
            %    param_def - struct containing varname, type, and
            %                function_handle for the parameter to be read
            %   group_info - Child struct from h5info call on h5_fname that
            %                points to the data of interest


            % Get station id from group info
            info_for_fileparts = group_info.Name;
            info_for_fileparts = strrep(info_for_fileparts, '/', filesep);
            [id, group] = fileparts(info_for_fileparts);
            id(1) = [];

            % Try reading data from cached entries
            data = obj.read_from_cached_entries(id, param_def.varname);

            % If data is not empty, it param was found in cached entries so
            % we can return
            if ~isempty(data); return; end

            % Find info for the param in group_info based on its varname
            param_info = GCSAL.GCSAL.h5info_find_children(group_info, param_def.varname);

            % If param not found, return empty
            if isempty(param_info)
                data = [];
                return
            end

            % load the parameter from the H5 file using param_info
            data = GCSAL.H5.load(obj.h5_fname, param_info);

            % Parameter is a char, convert uint8 to char
            if strcmp(param_def.type, 'char')
                data = char(data);
            end

            % If parameter was returned as a double from H5.load, but is
            % not defined as a double then convert to a single for
            % efficiency
            if isa(data, 'double') && ~strcmp(param_def.type, 'double')
                data = single(data);
            end

            % Apply function from parameter definition
            if ~isempty(param_def.function_handle)
                data = param_def.function_handle(data);
            end

            % Cache the data in entries
            % Since obj is a pointer (inherits handle class) we can cache
            % the data without returning it
            if ~strcmp(group, 'header')
                obj.cache_param(id, param_def.varname, data);
            end

        end


        function out = read_from_cached_entries(obj, id, param)
            % Read data from cached entries if it exists otherwise return
            % empty vector

            out = [];
            if isfield(obj.entries, id)
                if isfield(obj.entries.(id), param)
                    out = obj.entries.(id).(param);
                end
            end

        end

        function cache_param(obj, id, param, value)
            % Keep data in memory in entries struct for fast loading
            % Since obj is a pointer (inherits handle class) we can set the
            % obj.entries without returning it

            if obj.do_cache_entries
                obj.entries.(id).(param) = value;
            end
        end

        function add_entry_idx_to_headers(obj, i_headers)
            % add_entry_idx_to_headers(obj, i_headers)
            %   For the structs in obj.headers(i_headers), add the
            %   entry_idx field. This field is an indexing vector for
            %   the correspondence between header data and entry data.


            % Decide whether to do wait bar
            L = length(i_headers);
            do_waitbar = (L > 1) && ~obj.quiet_mode;

            % Open waitbar
            if do_waitbar
                h = waitbar(0, 'Calculating header to entry idx');
            end

            % Loop through indices in i_headres
            for i = 1:length(i_headers)

                % Extract the current header
                header = obj.headers(i_headers(i));

                % Check if entry_idx has already been added to this header
                if ~isfield(obj.headers, 'entry_idx') || isempty(header.entry_idx)

                    % Get the entry_idx that corresponds header to entry
                    % data
                    entry_idx = GCSAL.GCSAL.header_to_entry_idx(header);

                    % Convert entry_idx to the smallest possible type
                    obj.headers(i_headers(i)).entry_idx = GCSAL.IGRA.Param.convert_to_min_int(entry_idx);

                    % Update waitbar
                    if do_waitbar && mod(i, ceil(L/50)) == 0
                        msg = sprintf('%d/%d: Calculating header to entry idx for %s', i, L, header.id);
                        waitbar(i/L, h, msg);
                    end
                end
            end

            % Close wait bar
            if do_waitbar
                close(h)
            end
        end
    end

    methods (Static)

        function entries = filter_data_by_range(entries, range_fld, range)
            % entries = filter_data_by_range(entries, range_fld, range)
            %   Filter the data in entries to keep only instances where
            %   entries.(range_fld) is in range.
            %
            %   range can be a two element vector in which case data can be
            %   anywhere between range(1) and range(2) inclusive. Or range
            %   can be a scalar in which case data must be exactly equal to
            %   range.


            % Error check on length of range
            if length(range) ~= 1 && length(range) ~=2
                error('range should be length 1 or 2')
            end

            % Loop through stations
            L = length(entries);
            for i = 1:L

                % Extract the data in range_fld
                val = entries(i).(range_fld);

                % Get index vector for values that are in range
                if isscalar(range)
                    % If range is a scalar match exactly
                    idx = val == range;
                else
                    % If range is two element vector  match between
                    % range(1) and range(2) inclusive
                    idx = val >= range(1) & val <= range(2);
                end

                % Loop through each parameter in the current entry and apply index
                params = fieldnames(entries(i));
                for j = 1:length(params)

                    % Apply index
                    if size(entries(i).(params{j}), 1) ~= 1
                        entries(i).(params{j}) = entries(i).(params{j})(idx);
                    end
                end

            end

        end

        function [counts] = histcounts(entries, fld, edges, quiet_mode)
            % counts = histcounts(entries, edges, x_fld)
            %  Pulls the data located in fld for every element
            %  in entries and returns the counts in each bin constructed by
            %  bin edges

            if (quiet_mode)
                tic;
                fprintf('Counting %s...', fld );
            end
            
            % Force first and last bins to include -/+ inf
            edges(1) = -inf;
            edges(end) = inf;

            % Pre-allocate counts matrix with all zeros. There are one
            % fewer bins than edges on each side of bin grid
            Nrows = length(edges)-1;
            counts = zeros(1, Nrows);

            % Loop through each entry
            for i = 1:length(entries)

                % Extract data for x and y from entries
                x = entries(i).(fld);

                % Count x with bins defined by edges and add the
                % result to the existing counts
                counts = counts + histcounts(x, edges);

            end
            
            if (quiet_mode)
                fprintf('Complete in %.1f seconds\n', toc);
            end

            % The following is a simpler and more vectorized way to do the
            % same as above but surpisingly, testing proved that the above
            % is faster
%             counts = histcounts(vertcat(entries.(fld)), edges);

        end

        function [counts] = histcounts2(entries, x_fld, y_fld, x_edges, y_edges, quiet_mode)
            % counts = histcounts(entries, x_edges, y_edges, x_fld, y_fld)
            %  Pulls the data located in x_fld and y_fld for every element
            %  in entries and returns the counts in each bin constructed by
            %  bin edges defined in x_edges and y_edges

            if (~quiet_mode)
                tic;
                fprintf('Counting %s vs %$s...', x_fld, y_fld );
            end

            % Force first and last bins to include -/+ inf
            x_edges(1) = -inf;
            x_edges(end) = inf;
            y_edges(1) = -inf;
            y_edges(end) = inf;

            % Pre-allocate counts matrix with all zeros. There are one
            % fewer bins than edges on each side of bin grid
            Nrows = length(x_edges)-1;
            Ncols = length(y_edges)-1;
            counts = zeros(Nrows, Ncols);

            % Loop through each entry
            for i = 1:length(entries)

                % Extract data for x and y from entries
                x = entries(i).(x_fld);
                y = entries(i).(y_fld);

                % Count x and y in grid made of x/y edges and add the
                % result to the existing counts
                counts = counts + histcounts2(x, y, x_edges, y_edges);
            end

            if (~quiet_mode)
                fprintf('Complete in %.1f seconds\n', toc);
            end

            % The following is a simpler and more vectorized way to do the
            % same as above but surpisingly, testing proved that the above
            % is faster
%             x = vertcat(entries.(x_fld));
%             y = vertcat(entries.(y_fld));
%             counts = histcounts2(x, y, edges);

        end

        function [counts] = histcountsN(entries, resolutions, quiet_mode)
            % counts = histcounts(entries, edges, x_fld)
            %  Pulls the data located in fld for every element
            %  in entries and returns the counts in each bin constructed by
            %  bin edges

            if (quiet_mode)
                tic;
                msg = sprintf('  %s\n', resolutions.fld);
                fprintf('Counting... \n%s', msg)
            end
            
            % Pre-allocate counts matrix with all zeros. There are one
            % fewer bins than edges on each side of bin grid
            N = cell(size(resolutions));
            for j = 1:length(resolutions)
                N{j} = length(resolutions(j).edges) - 1;
            end
            counts = zeros(N{:}, 'uint32');


            % Determine whether to do waitbar
            L = length(entries);
            do_waitbar = (L > 1) && ~quiet_mode;
            if do_waitbar
                h = waitbar(0, 'Doing counts');
            end

            % Loop through each entry
            for i = 1:length(entries)

                % Initialize bins to cell array of proper size
                bins = cell(size(resolutions));

                % Loop through parameter resolutions
                for j = 1:length(resolutions)
                    fld = resolutions(j).fld;
                    edges = resolutions(j).edges;

                    % Force first and last bins to include -/+ inf
                    edges(1) = -inf;
                    edges(end) = inf;

                    % Extract data for x and y from entries
                    x = entries(i).(fld);

                    % Use discretize to determine bin index for every data
                    % point in x
                    bins{j} = discretize(x, edges);

                end

                % Add to counts for each idx made by bins
                try
                    idx = sub2ind(size(counts), bins{:});
                    idx(isnan(idx)) = [];
                    for j = 1:length(idx)
                        counts(idx(j)) = counts(idx(j)) + 1;
                    end
                catch e
                    fprintf(['Counts encountered an error with the following ' ...
                        'station so it was skipped: %s\n'], entries(i).id);
                    disp(e.identifier)
                    disp(e.message)
                    %                     keyboard
                end

                % Update waitbar
                if do_waitbar && mod(i, ceil(L/50)) == 0
                    msg = sprintf('%d/%d: Doing counts for %s', i, L, entries(i).id);
                    waitbar(i/L, h, msg);
                end

            end

            if (quiet_mode)
                fprintf('complete in %.1f seconds\n', toc)
            end

            % Close wait bar
            if do_waitbar
                close(h)
            end

            % Convert to smallest possible integer data type to save space
            counts = GCSAL.IGRA.Param.convert_to_min_int(counts);

        end

        function [pdf, cdf] = counts2pdf(counts, dim)
            % [pdf, cdf] = counts2percentile_pdf(counts, dim)
            %    For a matrix of counts, calculates the probability density
            %    function along dimension dim.
            %
            %    Conceptually if counts were a vector then
            %             pdf = counts / sum(counts)
            %
            %    This function does the calculation but in a vectorized way
            %    along dimension dim.
            %
            %    Additionally calculates the cumulative density function as
            %    cdf = cumsum(pdf, dim)

            % Get total counts in each row/column
            total_counts = sum(counts, dim);

            % Normalize counts to get probability density function
            pdf = bsxfun(@rdivide, counts, total_counts);

            % Accumulat pdf to get cdf
            cdf = cumsum(pdf, dim);

        end

        function bin_centers = get_bin_centers(bin_edges)
            % bin_centers = get_bin_centers(bin_edges)
            %   Returns bin_centers corresponding to midway value between
            %   each pair of edges in bin_edges

            left_edge = bin_edges(1:end-1);
            right_edge = bin_edges(2:end);

            bin_centers = (right_edge + left_edge) / 2;

        end

        function idx = find_keys(key_index, keys_to_find)
            % idx = find_keys(key_index, keys_to_find)
            %   Returns logical array idx indicating which elements in
            %   key_index match any string in keys_to_find.
            %
            %   key_index and keys_to_find can be strings, string matrices,
            %   or cell arrays of strings.
            %
            %   Warning is thrown if not all elements in keys_to_find are
            %   found in key_index

            % Handle case where keys_to_find is empty. This is
            % required because cellstr turns emptys strings into {''} which
            % is not an empty cell but rather has length 1
            if isempty(keys_to_find)
                keys_to_find = {};
            end

            % Ensure keys and key_array are cell arrays
            keys_to_find = cellstr(keys_to_find);
            key_index = cellstr(key_index);

            % Use ismember to get logical array for existence of each
            % element of key_index present in keys
            idx = ismember(key_index, keys_to_find);

            % Warning check that all keys_to_find were found
            keys_not_found_idx = ~ismember(keys_to_find, key_index(idx));

            if any(keys_not_found_idx)
                keys_not_found = keys_to_find(keys_not_found_idx);
                msg = sprintf('  %s\n', keys_not_found{:});
                warning('The following keys were not found: \n%s', msg)
            end
        end

        function info_matches = h5info_find_children(info, child_name)
            % info_matches = h5info_find_children(info, child_name)
            %   Returns a struct array corresponding to elements of
            %   info.Groups with a Name matching any string in
            %   child_name. child_name can be a string or cell array
            %   of strings. info should be part of the data structure
            %   returned by h5info.

            % Form search key by combining info.Name with subfolder
            search_key = GCSAL.H5.fullpath(info.Name, child_name);

            % Get all children names in info
            all_children_names = {info.Groups.Name};

            % Find matches
            i_children_match = GCSAL.GCSAL.find_keys( all_children_names, search_key);

            % Index in info.Groups
            info_matches = info.Groups(i_children_match);
        end

        function p = plot_stations(stations_to_plot, varargin)
            % p = plot_stations(stations_to_plot, varargin)
            %  Plots the lat/long coordinates of stations_to_plot. Any
            %  additional inputs to the plot function can be included in
            %  varargin. Returns a handle to the line object for the plot
            %  call.
            %
            %  stations_to_plot can either by a struct array with lat/long
            %  as fields or a list of station id strings

            p = plot([stations_to_plot.lon], [stations_to_plot.lat], varargin{:});

        end

        function edges = default_bin_edges(param_name)
            % edges = default_bin_edges(param_name)
            % Returns default values for bin edges given a parameter name

            switch param_name
                case 'gph'
                    edges = 0:1:30;
                case 'press'
                    edges = 0:1000:100000;
                case 'temp'
                    edges = -100:1:40;
                case 'rh'
                    edges = 0:1:100;
                case 'dpdp'
                    edges = 0:1:100;
                case 'wspd'
                    edges = 0:2:60;
                case 'wdir'
                    edges = 0:15:360;
                case 'month'
                    edges = 1:12;
                case 'day'
                    edges = 1:31;
                case 'hour'
                    edges = 0:24;
                case 'lat'
                    edges = -90:2:90;
                case 'lon'
                    edges = -180:4:180;
                otherwise
                    error('unrecognized fld: %s', param_name)
            end

        end

        function label = get_label(def)
            % label = get_label(def)
            %   Returns the from a definition struct. This label
            %   can be used for and x or y labels on a plot.

            label = def.description;
            if ~isempty(def.units)
                label = [label ' (' def.units ')'];
            end
        end

        function title_str = description_from_filters(fltr_flds, fltr_rngs)
            % Returns a string describing the filters in fltr_flds and
            % fltr_ranges that can be used on a plot title

            % Ensure fltr_flds and fltr_ranges are cell arrays
            fltr_flds = cellstr(fltr_flds);
            if ~iscell(fltr_rngs)
                fltr_rngs = {fltr_rngs};
            end

            % Initialize string as empty
            title_str = '';

            % Loop through each filter element
            for i = 1:length(fltr_flds)

                % Create a string describing the fitler applied based on
                % whether filter range was min/max or single value
                if length(fltr_rngs{i}) == 1
                    msg = sprintf('%s = %g, ', fltr_flds{i}, fltr_rngs{i});
                elseif length(fltr_rngs{i}) == 2
                    msg = sprintf('%s = [%g to %g], ', fltr_flds{i}, ...
                        fltr_rngs{i}(1), fltr_rngs{i}(2));
                else
                    error('fltr_ranges length expected to be 1 or 2')
                end

                % Append msg to the title_str
                title_str = [title_str msg];  %#ok<AGROW>
            end

            % If title_str is empty remove the new line character at the
            % end of the string
            if ~isempty(title_str)
                title_str(end-1:end) = [];
            end
        end

        function stations_out = stations_intersect(stations1, stations2)
            % out = stations_intersect(stations1, stations2)
            %   Returns the intersect of struct arrays stations1 and
            %   stations2 based on their id fields


            [~, idx] = GCSAL.GCSAL.struct_set_operation(...
                stations1, stations2, 'id', @intersect);
            stations_out = stations1(idx);

        end

        function stations_out = stations_union(stations1, stations2)
            % out = stations_union(stations1, stations2)
            %   Returns the union of struct arrays stations1 and
            %   stations2 based on their id fields

            [~, ia, ib] = GCSAL.GCSAL.struct_set_operation(...
                stations1, stations2, 'id', @union);
            stations_out = [stations1(ia) stations2(ib)];

        end

        function stations_out = stations_setxor(stations1, stations2)
            % out = stations_setxor(stations1, stations2)
            %   Returns the setxor of struct arrays stations1 and
            %   stations2 based on their id fields

            [~, ia, ib] = GCSAL.GCSAL.struct_set_operation(...
                stations1, stations2, 'id', @setxor);
            stations_out = [stations1(ia) stations2(ib)];

        end

    end


    methods (Static, Access = 'private')

        function [c, ia, ib] = struct_set_operation(struct1, struct2, fld, operation)
            % Set the set operation on the data at struct1.(fld) and
            % struct2.(fld). operation can be intersect, union, or setxor.
            % Returns the outputs of the operation

            str1 = vertcat(struct1.(fld));
            str2 = vertcat(struct2.(fld));
            [c, ia, ib] = operation(str1, str2, 'rows');

        end

        function stations = initialize_stations(headers)
            % stations = initialize_stations(obj)
            %  Creates the stations struct array from

            % Initialize stations struct with headers.id and NaN for lat
            % and lon
            stations = struct('id', {headers.id}, 'lat', NaN, 'lon', NaN);

            % Loop through all headers
            for i = 1:length(headers)

                % Get current header
                header = headers(i);

                % If lat and lon are sclars then use them
                if isscalar(header.lat) && isscalar(header.lon)
                    stations(i).lat = header.lat;
                    stations(i).lon = header.lon;
                else
                    % If lat and lon are not scalars then check how big the
                    % the biggest difference is
                    x = abs(max(header.lat) - min(header.lat));
                    y = abs(max(header.lon) - min(header.lon));
                    d = sqrt(x^2 + y^2);

                    % As long as the difference isn't too big, use the mode
                    if d < .2
                        stations(i).lat = mode(header.lat);
                        stations(i).lon = mode(header.lon);
                    else
                        % Otherwise we will the station lat/lon as NaN.
                        % Stations that begin with ZZ are expected to move
                        % around a lot, but any station besides ZZ, report
                        % tit.
                        if ~strcmp(header.id([1 2]), 'ZZ')
                            fprintf(['Location for station %s was not used ' ...
                                'because it moved around by ~%g deg\n'], header.id, d);
                        end
                    end
                end
            end
        end

        function idx = get_entry_idx_in_range(header, range, range_fld)

            % Get header values at range_fld
            val = header.(range_fld);

            % Find index in header where value is in range
            idx_header = val >= range(1) & val <= range(2);

            % Extract array of entry lengths for header values in range
            numlevs_in_range = header.numlevs(idx_header);

            % Next we want to calculate the index offset for the start of
            % each header in range. First we need the index to the start of
            % every entry which we get by first getting the cumulative sum
            % of all numlevs
            cumulative_numlevs_all = cumsum(header.numlevs);

            % Then to get an idx offset we just start at 1 and remove the
            % last element of the cumulative sum
            idx_offset_all = [1; cumulative_numlevs_all(1:end-1)];


            % Now getting the index offset to the in range headers is just
            % indexing into idx_offset_all
            idx_offset = idx_offset_all(idx_header);

            % Initialize idx to be size of total from numlevs
            idx = zeros(sum(numlevs_in_range), 1);

            % initialize start index to 1
            start = 1;
            for i = 1:length(numlevs_in_range)

                % finish index is start + current # of entries - 1
                finish = start + numlevs_in_range(i) - 1;

                % Assign idx based on idx_offset and current number of
                % entries
                idx(start:finish) = (1:numlevs_in_range(i)) + idx_offset(i);

                % For next loop start where we left off
                start = finish+1;

            end

        end

        function entry_idx = header_to_entry_idx(header)
            % entry_idx = header_to_entry_idx(header)
            %   Returns an indexing vector for the correspondence between
            %   header data and entry data.
            %
            %   Each element in header corresponds to many elements in the
            %   entry data for the same station. This function determines
            %   the indexing vector that allows you to expand the header
            %   data to the length of the entry data with the header data
            %   copied for each instance of entry data for which it
            %   corresponds.


            % Initialize idx to be size of total from numlevs
            entry_idx = zeros(sum(header.numlevs), 1);

            % initialize start index to 1
            start = 1;
            for i = 1:length(header.numlevs)

                % finish index is start + current # of entries - 1
                finish = start + header.numlevs(i) - 1;

                % Assign idx based on idx_offset and current number of
                % entries
                entry_idx(start:finish) = i;

                % For next loop start where we left off
                start = finish+1;
            end
        end

        function station_ids = station_id_str(stations)
            % station_ids = station_id_str(stations)
            %   Returns a character array whose rows correspond to the the
            %   id field from the struct array stations.
            %
            %   If stations is already a character array returns stations.

            % Handle empty stations
            if ~exist('stations', 'var') || isempty(stations)
                station_ids = '';
                return
            end

            % Handle case where station_ids is input as a struct array by
            % extracting station id character array
            if isstruct(stations)
                station_ids = vertcat(stations.id);
            elseif ischar(stations)
                station_ids = stations;
            else
                error('stations expeted to be struct or char')
            end
        end


    end
end
