function create_and_write(filename, datasetname, data)
% Copyright (c) Facebook, Inc. and its affiliates.
%
% This source code is licensed under the MIT license found in the
% LICENSE file in the root directory of this source tree.
%
% [] = create_and_write(filename, datasetname, data)
%   Uses h5create and h5write to write data to datasetname in filename
%   If data is empty, does nothing


if ~isempty(data)
    h5create(filename, datasetname, size(data), 'Datatype', class(data))
    h5write( filename, datasetname, data)
end

end
