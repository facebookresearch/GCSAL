function [  ] = recursive_load( h5_file, info )
% Copyright (c) Facebook, Inc. and its affiliates.
%
% This source code is licensed under the MIT license found in the
% LICENSE file in the root directory of this source tree.
%
% Recursively load all data in h5_file that is a child of the Groups in
% info. This is for test purpose only


if isempty(info.Datasets)
    for i = 1:length(info.Groups)
        GCSAL.H5.recursive_load( h5_file, info.Groups(i));
    end
else
    GCSAL.H5.load(h5_file, info);
end


end
