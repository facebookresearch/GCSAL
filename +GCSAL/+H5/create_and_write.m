function create_and_write(filename, datasetname, data)
% Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved
%
% [] = create_and_write(filename, datasetname, data)
%   Uses h5create and h5write to write data to datasetname in filename
%   If data is empty, does nothing


if ~isempty(data)
    h5create(filename, datasetname, size(data), 'Datatype', class(data))
    h5write( filename, datasetname, data)
end

end
