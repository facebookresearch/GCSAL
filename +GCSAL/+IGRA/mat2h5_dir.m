function mat2h5_dir(in_dir, h5_filename, append, filespec)
% Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved
%
% [] = mat2h5_dir(in_dir, h5_filename, append, filespec)
%   Loads data in .mat files that match filespec and writes the data to
%   h5_filename. The flag append controls whether data will be appended to
%   h5_filename or if a new h5_filename will be created from scratch.
%
%   All inputs are optional and will revert to a default value if not provided.
%
% INPUTS
%          in_dir - Directory to look for .mat data files.
%                   Default: current working directory
%     h5_filename - File to write data to
%                   Default: './gcsal.h5'
%          append - Flag whether to append to existing h5 file or start new
%                   Default: true
%        filespec - filespec used to identify .mat files
%                   Default: '*-data.txt.mat'


% Set default values
if ~exist('in_dir', 'var')
    in_dir = pwd;
end

if ~exist('h5_filename', 'var')
    h5_filename = fullfile(in_dir, 'gcsal.h5');
end

if ~exist('append', 'var')
    append = true;
end

if ~exist('filespec', 'var')
    filespec = '*-data.txt.mat';
end

% Find all files matching filespec
filespec = fullfile(in_dir, filespec);
fileObj = dir(filespec);

% Sort files by size
[~, indices] = sort([fileObj.bytes], 'ascend');
fileObj = fileObj(indices);

% Calculate total size of all files in MB
N_files = length(fileObj);

% Get h5 file info
names = {};
if exist(h5_filename, 'file')
    if append
        info = h5info(h5_filename);
        names = {info.Groups.Name};
    else
        delete(h5_filename)
    end
end

fprintf('Processing %d files in %s\n', N_files, in_dir);

% Initialize counters
time_so_far = 0;
t1 = tic;

% Iterate through files and read data
for i = 1:N_files
    t2 = tic;
    curr = fileObj(i);
    curr_name = curr.name;
    fprintf('Reading %s %04d/%04d ', curr_name, i, N_files);

    % Construct station id from .mat filename
    station_id = ['/' curr_name(1:end-13)];

    % Check that station_id is not already in H5 file
    if ~ismember(station_id, names)
        % Load .mat file and write all params in header and entries to h5
        mat_filename = fullfile(in_dir, curr_name);
        GCSAL.IGRA.mat2h5( mat_filename, h5_filename, station_id )
    end

    curr_time = toc(t2);
    time_so_far = time_so_far + curr_time;
    prct_complete = i/N_files;
    total_time = time_so_far/prct_complete;
    avg_rate = time_so_far/i;
    fprintf('%3.0f%% Curr: %.2f sec, Avg: %.2f sec, %4.0f/%.0f sec\n', ...
        prct_complete*100, curr_time, avg_rate, time_so_far, total_time);
end
toc(t1)

end
