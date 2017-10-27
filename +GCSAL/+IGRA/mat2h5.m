function [ ] = mat2h5( mat_filename, h5_filename, station_id )
% Load mat_filename and write the contents in h5 format to h5_filename with
% station_id as h5 path root

% Copyright (c) 2017-present, Facebook, Inc.
% All rights reserved.
%
% This source code is licensed under the BSD-style license found in the
% LICENSE file in the root directory of this source tree. An additional grant
% of patent rights can be found in the PATENTS file in the same directory.


% load mat file
mat = load(mat_filename);

% Loop through the fields
flds = fieldnames(mat);
for i = 1:length(flds)

    % Construct h5 path
    h5path = GCSAL.H5.fullpath(station_id, flds{i});

    % h5write all fields
    h5write_all_params(mat.(flds{i}), h5_filename, h5path)
end

end

function h5write_all_params(data, h5filename, datasetname)
% Loop through all parameters in data and call h5write with h5filename and
% datasetname. data should be a struct of GCSAL.IGRA.Param objects

% Loop through fields
flds = fieldnames(data);
for i = 1:length(flds)
    % Extract Param object from data struct
    curr_param = data.(flds{i});

    % Call h5write method of Param object
    curr_param.h5write(h5filename, datasetname);
end
end
