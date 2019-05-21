function [ out ] = load( filename, info )
% Copyright (c) Facebook, Inc. and its affiliates.
%
% This source code is licensed under the MIT license found in the
% LICENSE file in the root directory of this source tree.
%
% [ out ] = load( filename, info )
%   Load and uncompress data in h5 file filename at dataset find in info.
%   info is a specific Group found in the struct return by h5info()


% List of expected flds
flds = {'data', 'len', 'idx', 'i_unique'};

% Extract data set names, sort for stable comparison below
dsets = {info.Datasets.Name};

% Error check for unexpected fields
unexpected_flds = intersect(setxor(dsets, flds), dsets);
if ~isempty(unexpected_flds)
    msg = sprintf('  %s\n', unexpected_flds{:});
    error(['Unexpected datasets encountered:\n%s' ...
    'Expected datasets are data, len, and idx.'], msg) %#ok<SPERR>
end

% Loop through expected flds and read data if available
for i = 1:length(flds)
    if ismember(flds{i}, dsets)
        compressed_data.(flds{i}) = h5read(filename, GCSAL.H5.fullpath(info.Name, flds{i}));
    end
end

% Uncompress the data based on the datasets that were found
out = uncompress_data(compressed_data);

end


function out = uncompress_data(in)
% Returns uncompressed vector represented in uncompressed data
% struct in.

% First use i_unique to expand data if it was included
if isfield(in, 'i_unique')

    % Use i_unique to expand data to size of idx
    in.data = in.data(in.i_unique,:);

    % Remove i_unique field
    in = rmfield(in, 'i_unique');
end

% Extract fields from compressed data struct
flds = sort(fieldnames(in));

% Select decompression method based on fields found in uncompressed data
% struct
if isequal(flds, {'data'; 'idx'; 'len'})

    out = uncompress_data_idx_len(in.data, in.idx, in.len);

elseif isequal(flds, {'data'; 'len'})

    out = uncompress_data_len(in.data, in.len);

elseif isequal(flds, {'data'})
    % If only data is given, there is nothing to do
    out = in.data;
elseif isequal(flds, {'len'})
    % If only length is given return properly sized NaN vector
    out = NaN(in.len, 1);
else
    msg = sprintf('  %s\n', flds{:});
    error('Unexpected combination of datasets found:\n%s', msg)
end

end

function out = uncompress_data_idx_len(data, idx, len)
% Returns uncompressed vector from data, idx, and len compressed
% representation

    % Extract size of data input
    [N_data_rows, N_data_cols] = size(data);

    % Initialize out to properly sized NaN
    out = NaN(len, N_data_cols);

    % Uncompress idx
    idx = GCSAL.IGRA.Param.uncompress_idx(idx, len);

    % data should be the same size as idx unless there was only a single
    % data value, in which case it should replicated to match the length of
    % idx
    if N_data_rows == 1
        data = repmat(data, length(idx),1);
    end

    if size(data, 1) ~= length(idx)
        error('Numer of rows in data does not match length of idx')
    end

    % index data into out
    out(idx,:) = data;
end

function out = uncompress_data_len(data, len)
% Returns uncompressed vector from data, and len compressed
% representation

    if size(data,1) == 1
        % Only a single value for data given so replicate to match len
%         out = repmat(data, len, 1);
        out = data;
    else
        % Verify that length of data matches len
        if size(data,1) ~= len
            error('Number of rows in data does not match len')
        end

        % Nothing else to do data is already uncompressed
        out = data;
    end
end
