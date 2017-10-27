function [  ] = recursive_load( h5_file, info )
% Recursively load all data in h5_file that is a child of the Groups in
% info. This is for test purpose only

% Copyright (c) 2017-present, Facebook, Inc.
% All rights reserved.
%
% This source code is licensed under the BSD-style license found in the
% LICENSE file in the root directory of this source tree. An additional grant
% of patent rights can be found in the PATENTS file in the same directory.


if isempty(info.Datasets)
    for i = 1:length(info.Groups)
        GCSAL.H5.recursive_load( h5_file, info.Groups(i));
    end
else
    GCSAL.H5.load(h5_file, info);
end


end
