function create_and_write(filename, datasetname, data)
% [] = create_and_write(filename, datasetname, data)
%   Uses h5create and h5write to write data to datasetname in filename
%   If data is empty, does nothing

% Copyright (c) 2017-present, Facebook, Inc.
% All rights reserved.
%
% This source code is licensed under the BSD-style license found in the
% LICENSE file in the root directory of this source tree. An additional grant
% of patent rights can be found in the PATENTS file in the same directory.


if ~isempty(data)
    h5create(filename, datasetname, size(data), 'Datatype', class(data))
    h5write( filename, datasetname, data)
end

end
