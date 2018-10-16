function [  ] = recursive_load( h5_file, info )
% Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved
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
