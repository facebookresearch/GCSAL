classdef Param
    % Param(param_def, text_mat)
    %    Reads the IGRA formatted text in text_mat according to the
    %    format definition in param_def. IGRA text data uses fixed column
    %    widths so a single Param reads only data from the fixed
    %    columns specified in param_def corresponding with a single data
    %    parameter.
    %
    %    The param_def must supply the following struct fields:
    %       varname - name used to store the data in an H5 file. Must
    %                 resolve to a valid Matlab struct field.
    %          type - 'uint8', 'uint16', 'uint32', 'uint64', 'int8',
    %                 'int16', 'int32', 'int64', 'single', 'double', or
    %                 'char'
    %       col_idx - index vector corresponding to the columns of the
    %                 source file relevant to the current entry
    %      bad_vals - Cell array of strings listing the code used to
    %                 indicate missing or erroneous data in the IGRA file
    %
    % INPUTS
    %  param_def - a struct containing varname, type, col_idx,
    %              bad_vals
    %   text_mat - a matrix of characters converted to uint8. Each
    %              row corresponding to a row of text from an IGRA
    %              formatted text file
    %
    % PROPERTIES
    %   data - stored data representation of text.
    %   idx  - indexing vector for data
    %   len  - length of uncompressed data
    %   def  - format definition struct

    % Copyright (c) 2017-present, Facebook, Inc.
    % All rights reserved.
    %
    % This source code is licensed under the BSD-style license found in the
    % LICENSE file in the root directory of this source tree. An additional grant
    % of patent rights can be found in the PATENTS file in the same directory.


    properties
        data % stored data representation of text
        idx  % indexing vector for data
        len  % length of uncompressed data
        def  % format definition struct
        i_unique % indexing vector if data was compressed based on unique parameters
    end

    methods
        function obj = Param(param_def, text_mat)
            % Param Constructor
            %
            % INPUTS
            %  param_def - a struct containing varname, type, col_idx,
            %              bad_vals
            %   text_mat - a matrix of characters converted to uint8. Each
            %              row corresponding to a row of text from an IGRA
            %              formatted text file

            % Error check fields in format definition struct
            found_flds = fieldnames(param_def);
            expected_flds = {'varname', 'type', 'col_idx', 'bad_vals'};
            if ~all(ismember(expected_flds, found_flds))
                error('def struct must contain varname, type, col_idx, and bad_vals fields')
            end

            % Store format definition for the data entry
            obj.def = param_def;

            % Error check the text_mat
            if ~isa(text_mat, 'uint8')
                error('text_mat must be a matrix of characters converted to uint8')
            end

            % Read the string matrix text_mat
            obj = obj.read_columns(text_mat);

        end

        function obj = read_columns(obj, text_mat)

            % Extract the relative columns of text
            txt = text_mat(:, obj.def.col_idx);

            % Compress txt by removing rows that match bad_vals and
            % creating an idx vector for mapping remaining rows (unless the
            % compression actually increases the data size in which case
            % leave it as is)
            [txt, obj.idx, obj.len] = GCSAL.IGRA.Param.compress_txt(txt, obj.def.bad_vals);

            % Convert txt to data based on the type in the format
            % definition
            [obj.data, obj.i_unique] = GCSAL.IGRA.Param.txt2data(txt, obj.def.type);

        end

        function h5write(obj, filename, dataset_prefix)
            % Write the data to h5 file filename with h5 path
            % dataset_prefix

            % Make h5 path for this varname
            h5_path = GCSAL.H5.fullpath(dataset_prefix, obj.def.varname);

            % Write data, idx, and len
            obj.h5_write_param(filename, h5_path, 'data')
            obj.h5_write_param(filename, h5_path, 'idx')
            obj.h5_write_param(filename, h5_path, 'i_unique')

            % Don't need len if data is already full length
            if length(obj.data) ~= obj.len
                obj.h5_write_param(filename, h5_path, 'len')
            end

        end

        function h5_write_param(obj, filename, dataset_prefix, var_name)
            % Helper for calling H5.create_and_write with proper h5 path
            h5_path = GCSAL.H5.fullpath(dataset_prefix, var_name);
            GCSAL.H5.create_and_write(filename, h5_path, obj.(var_name))
        end

    end

    methods (Static)


        function str = pad_left(str, desired_length)
            % Prepends str with enough spaces to make a string of
            % desired_length

            str_length = length(str);
            if str_length > desired_length
                error('str: %s is already longer than str_length: %f', str, desired_length)
            end

            % Pad beginning of str with blanks
            prefix = blanks(desired_length - str_length);
            str = [prefix str];

        end

        function val = convert_to_min_int(val)
            % Create an idx variable with values 1:max_value and of the
            % smallest type that can store max_value

            max_val = max(val(:));
            if isempty(max_val)
                max_val = 1;
            end
            if max_val < 1 || rem(max_val, 1) ~= 0
                error('max_value: %f must be a whole number greater than 0')
            end

            if max_val < intmax('uint8')
                val = uint8(val);
            elseif max_val < intmax('uint16')
                val = uint16(val);
            elseif max_val < intmax('uint32')
                val = uint32(val);
            elseif max_val < intmax('uint64')
                val = uint64(val);
            else
                error(['max_val: %f exceeds the maximum value than ' ...
                    'can be stored as an integer'], max_val)
            end
        end

        function [data, i_unique] = txt2data(txt, type)
            % Convert txt to type

            % Return on empty txt
            if isempty(txt)
                data = []; i_unique = [];
                return
            end

            % Since txt may be very repetitive it is more efficient to find
            % unique rows before trying to convert
            [unique_txt, ~, i_unique] = unique(txt, 'rows');
            i_unique = GCSAL.IGRA.Param.convert_to_min_int(i_unique);

            % Switch on whether data is int, real, or char
            switch type

                % Integer types
                case {'int8', 'int16', 'int32', 'uint8', 'uint16', 'uint32'}
                    unique_val = GCSAL.IGRA.Param.str2int(unique_txt, type);


                    % Float types
                case {'single', 'double'}
                    unique_val = GCSAL.IGRA.Param.str2float(unique_txt, type);

                    % String types
                case 'char'
                    unique_val = unique_txt;

                    % Anything else
                otherwise
                    error('Unrecognized type: %s', type)
            end


            if size(unique_val, 1) == 1
                % If there is only one unique value, return that value with
                % i_unique empty
                i_unique = [];
                data = unique_val;
            else
                % Otherwise determine if we can save space by representing
                % data in unique form

                if GCSAL.IGRA.Param.compare_bytes_unique(i_unique, unique_val)
                    % keep data in unique form
                    data = unique_val;
                else
                    % Inverse the unique call and make data the full
                    % uncompressed data
                    data = unique_val(i_unique,:);

                    % Revert i_unique to empty vector
                    i_unique = [];
                end
            end
        end

        function do_compression = compare_bytes_unique(i_unique, unique_val)
            %  do_compression = compare_bytes_unique(i_unique, unique_val)
            %    Returns true if bytes needed to store i_unique and
            %    unique_val is less than bytes needed to store
            %    unique_val(i_unique,:).
            %
            %    This is generally true if i_unique is a smaller data type
            %    than unique_val and unique_val is significantly shorter
            %    than unique_val(i_unique,:)

                % Figure out bytes requred to store data with i_unique and
                % unique_val
                bytes_i_unique= whos('i_unique');
                bytes_i_unique = bytes_i_unique.bytes;

                bytes_unique = whos('unique_val');
                bytes_unique = bytes_unique.bytes;

                bytes_compressed = bytes_i_unique + bytes_unique;

                % Figure out bytes to store uncompressed data
                % Do this by multiplying bytes_unique by the
                % ratio of the size of data for compressed and uncompressed
                length_uncompressed = length(i_unique);
                length_compressed = size(unique_val, 1);
                bytes_uncompressed = bytes_unique*length_uncompressed/length_compressed;

                % Return comparison of bytes_compressed and uncompressed
                do_compression = bytes_compressed < bytes_uncompressed;

        end

        function int = str2int(str_mat, type)
            % Custom vectorized str2int function. Processes the character
            % matrix str_mat by columns, converting the characters in each
            % column to a number and then adding up the numbers from each
            % column. This is faster than str2double or str2num because it
            % relies on the characters in str_mat being well behaved.
            %
            % Additional the character matrinx str_mat in this case is
            % represented as uint8 for efficiency
            %
            % Assumptions:
            %   - The only characters in str_mat are ' -0123456789'
            %   - Every row of str_mat is a single number
            %   - Every row of str_mat is equal width with left padding
            %   - Every row is ordered blanks, negative sign, numerals from
            %     left to right

            % Convert char to uint8 if necessary
            if isa(str_mat, 'char')
                str_mat = uint8(str_mat);
            end

            % Get the size of str_mat for preallocation
            [rows, cols] = size(str_mat);

            % Pre-allocate int with zeros of the correct type
            int = zeros(rows,1, type);

            % Initialize place_factor with ones. This rep100resents the value
            % of the current column (like 1, 10, 100, 1000, etc.)
            place_factor = ones(rows,1, type);

            % Convert each character of str_mat to it's integer numeral
            % (blanks become 0 and negative sign becomes -1)
            numerals = GCSAL.IGRA.Param.char2numerals(str_mat);

            % Convert numerals to proper type
            numerals = cast(numerals, type);

            % Go column by column to sum each row vector of numerals into a
            % single value per row in a vectorized manner. Start with right
            % most column for the ones place
            for i = cols:-1:1

                % Extract the current column, going right to left
                curr_col = numerals(:,i);

                % Find any numerals in the current column that are
                % negative
                is_neg = curr_col == -1;

                % negate where isneg
                int(is_neg) = -1*int(is_neg);

                % Apply the place_factor multiplier to the curr_col where
                % not is_neg
                int(~is_neg) = int(~is_neg) + place_factor(~is_neg) .* curr_col(~is_neg);

                % Increment place_factor by x10 since we are in base 10
                place_factor = place_factor*10;
            end

            % error check on data type
            if any(int(:) >= intmax(type))
                error('Data overflow for type: %s', type)
            end
        end

        function int = char2numerals(char_mat)
            % converts matrix of characters char_mat to a matrix of
            % numerals. The characters of char_mat are represented as uint8
            % for efficiency
            %
            % The only acceptable characters in char_mat are ' -0123456789'

            % Convert char to uint8 if necessary
            if isa(char_mat, 'char')
                char_mat = uint8(char_mat);
            end

            % Create key/value map for converting uint8 characters to
            % numerals. blanks become 0 and negative signs become -1
            key = uint8(' -0123456789');
            val = [0 -1 0 1 2 3 4 5 6 7 8 9];

            % Match all characters in char_mat with key
            [test, key_idx] = ismember(char_mat, key);
            if ~all(test(:))
                error('Encountered an unrecognized character. The only recognized characters are -0123456789 and blank')
            end

            % Use the matching indices from ismember to index to the proper
            % values in val
            int = val(key_idx);

            % int defaults to a row vector if the char_mat is a column
            % vector, in this case transpose to size of in matches size of
            % char_mat
            if size(char_mat,2) == 1
                int = int';
            end

            % Error check on size
            if ~all(size(int) == size(char_mat))
                error('size error')
            end

        end

        function float = str2float(str_mat, type)
            % Uses str2num to convert the character matrix str_mat to a
            % number then applies the necessary type conversion

            float = str2num(char(str_mat)); %#ok<ST2NM>
            float = cast(float, type);
        end

        function int_array = bits2ints(bits)
            % Convert a logical array bits to an array of integers. A
            % logical array actually uses 1 byte (8 bits) to represent each
            % boolean. By converting the logical array to integers you can
            % reduce the size in memory by 8 times.

            % convert to uint32 and column vector
            bits = uint32(logical(bits(:)));

            % Want to convert to array of uint32 whereeach integer has 32
            % bits. We will reshape bits be an Nx32 matrix, but bits may
            % not be divisible by 32 so we will pad with a prefix of zeros
            % as necessary.
            L = 32;

            % First, keep track of how many bits we started with so the
            % first bit can be disambiguated upon decoding
            N_bits = uint32(length(bits));

            % Calculate how much padding is needed
            padding = L -rem(length(bits), L);
            if padding == L; padding = 0; end

            % Prepend padding
            bits = [zeros(padding, 1); bits];

            % Reshape
            bits = reshape(bits, L, [])';

            % convert each group of 32 bits to a uint32 by adding up each
            % column multiplied 2^i_col
            int_array = zeros(size(bits,1), 1, 'uint32');
            for i_col = 1:L
                int_array = int_array + bits(:,L+1-i_col).*2.^(i_col-1);
            end

            % Finally prefix N_bits at the beginning so that this array can
            % be decoded without ambiguity about the bits that were padded
            int_array = [N_bits; int_array];

        end

        function bits = ints2bits(int_array)
            % Convert an array of integers (assumed to be uint32 in this
            % implementation) to bits represented as a logical array.

            % First integer in int_array contains the # of bits stored
            N_bits = int_array(1);
            int_array(1) = [];

            % To convert an integer to bits we need to keep dividing by 2
            % and checking the remainder
            L = 32;
            int_array = double(int_array); % this ensures division by 2 works properly
            bits = zeros(32, length(int_array));
            for i = 1:L
                % Check the remainder when dividing by 2
                bits(L+1-i, :) = mod(int_array, 2);

                % Reduce by half
                int_array = floor(int_array/2);
            end

            % Ensure column vector
            bits = bits(:);

            % We may have some extra bits that were added by the padding
            % process caused by N_bits not being exactly divisible by 32.
            % This step removes any extra bits that were prefixed
            bits_to_remove = length(bits) - N_bits;
            bits(1:bits_to_remove) = [];

        end

        function idx = compress_idx(idx_good, original_length)
            % idx can be represented as either a list of the good indices,
            % a list of the bad indices, or a logical array of bits. This
            % function calculates which one will be most efficient

            % If idx_good is empty, then we do not need to continue
            if isempty(idx_good)
                idx = [];
                return
            end

            % If idx_good is the same length as original_length then no
            % compression occured and idx can be returned empty
            if length(idx_good) == original_length
                idx = [];
                return
            end

            % Crete a logical index version of idx
            idx_logical = zeros(original_length,1, 'uint8');
            idx_logical(idx_good) = 1;

            % Create the inverse idx
            idx_bad = find(~idx_logical);

            % Convert to the min int type to save space
            idx_bad = GCSAL.IGRA.Param.convert_to_min_int(idx_bad); %#ok<FNDSB>

            % Convert logical from bits to int to save space
            idx_logical = GCSAL.IGRA.Param.bits2ints(idx_logical);

            % Find which idx corresponds to the fewest bytes
            var_info(1) = whos('idx_good');
            var_info(2) = whos('idx_bad');
            var_info(3) = whos('idx_logical');
            [~, idx_type] = min([var_info.bytes]);

            % Set idx based on the idx_type found to be most efficient
            switch idx_type
                case 1
                    idx = idx_good;
                case 2
                    idx = idx_bad;
                case 3
                    idx = idx_logical;
            end

            % Prepend the idx_type to the idx array
            idx = [idx_type; idx];

        end

        function idx_out = uncompress_idx(idx_in, original_length)

            % First value in idx_in should be encoding of index type
            idx_type = idx_in(1);
            idx_in(1) = [];

            % Switch on index type
            switch idx_type

                case 1
                    % idx_in is already good_vals
                    idx_out = idx_in;

                case 2
                    % idx_in is bad_vals and needs to be inversed
                    logical_array = ones(original_length, 1);
                    logical_array(idx_in) = 0;
                    idx_out = find(logical_array);

                case 3
                    % idx_in is logical bits represented as integer aray
                    logical_array = GCSAL.IGRA.Param.ints2bits(idx_in);
                    idx_out = find(logical_array);

                otherwise
                    error('Unrecognized idx type')
            end

        end

    end

    methods (Static, Access = 'private')
        function data = unique_inverse(unique_val, i_unique)
            % Use i_unique to index unique_val and return the original
            % ordering of the data before unique was called. If, however,
            % there is only one unique element, then return just the unique
            % value

            if size(unique_val, 1) == 1
                data = unique_val;
            else
                data = unique_val(i_unique,:);
            end
        end

        function [txt, idx, original_length] = compress_txt(orig_txt, bad_vals)

            % Cache # of rows in orig_txt
            original_length = size(orig_txt, 1);
            original_length = GCSAL.IGRA.Param.convert_to_min_int(original_length);

            % Create compressed version of txt by removing rows that match
            % bad_vals and keeping track of the idx for the remaining rows
            [txt, idx] = GCSAL.IGRA.Param.remove_bad_vals(orig_txt, bad_vals);

            % Compress idx as efficiently as possible
            idx = GCSAL.IGRA.Param.compress_idx(idx, original_length);


        end

        function [txt, idx] = remove_bad_vals(txt, bad_vals)

            % Get size
            [N_rows, N_cols] = size(txt);

            % Initialize indexing to incldue all rows
            idx = (1:N_rows)';
            idx = GCSAL.IGRA.Param.convert_to_min_int(idx);

            % Removes rows of txt that match any of the strings in
            % bad_vals

            % Loop through the list of bad_vals in the format definition
            % Each bad_val should be a string
            for i = 1:length(bad_vals)

                % Ensure bad_val is the correct width by padding
                curr_bad_val = GCSAL.IGRA.Param.pad_left(bad_vals{i}, N_cols);

                % Convert char to uint8 to match format of txt matrix
                curr_bad_val = uint8(curr_bad_val);

                % Find compare curr_bad_val to each row of txt
                matching_characters = bsxfun(@eq, curr_bad_val, txt);

                % Find rows where all characters match
                matching_rows = all(matching_characters, 2);

                % Remove matching_rows from txt and idx
                txt(matching_rows, :) = [];
                idx(matching_rows) = [];
            end
        end
    end
end
