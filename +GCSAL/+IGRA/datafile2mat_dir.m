function datafile2mat_dir(in_dir, out_dir, overwrite_mat, filespec)
% Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved
%
% [] = datafile2mat_dir(in_dir, out_dir, use_mat_if_found, filespec)
%   Parses IGRA data files that match filespec and saves a .mat file
%   containing the data in out_dir. The flag overwrite_mat controls whether
%   IGRA files that already have a matching .mat file will be skipped or
%   overwritten.
%
%   All inputs are optional and will revert to a default value if not provided.
%
% INPUTS
%          in_dir - Directory to look for IGRA data files.
%                   Default: current working directory
%         out_dir - Directory to save output .mat files
%                   Default: current working directory
%   overwrite_mat - Flag whether to skip existing mat files or overwrite
%                   Default: true
%        filespec - filespec used to identify IGRA data files
%                   Default: '*-data.txt'


% Set default values
if ~exist('in_dir', 'var')
    in_dir = pwd;
end

if ~exist('out_dir', 'var')
    out_dir = pwd;
end

if ~exist('overwrite_mat', 'var')
    overwrite_mat = true;
end

if ~exist('filespec', 'var')
    filespec = '*-data.txt';
end

% Find all files ending in "-data.txt"
filespec = fullfile(in_dir, filespec);
fileObj = dir(filespec);

% % Sort files by size
% [~, indices] = sort([fileObj.bytes], 'ascend');
% indices = indices([1:300]);
% fileObj = fileObj(indices);

% Calculate total size of all files in MB
total_MB = sum([fileObj.bytes])/1e6;
N_files = length(fileObj);

fprintf('Reading %d files totalling %.1f MB in %s\n', N_files, total_MB, in_dir);

% Initialize counters
read_MB = 0;
time_so_far = 0;
t1 = tic;

% Iterate through files and read data
for i = 1:N_files
    t2 = tic;
    curr = fileObj(i);
    curr_name = curr.name;
    curr_MB = curr.bytes/1e6;
    fprintf('%04d: Reading %s, %5.1f MB', i, curr_name, curr_MB);

    GCSAL.IGRA.datafile2mat( fullfile(in_dir, curr_name), overwrite_mat, out_dir );

    read_MB = read_MB + curr_MB;
    curr_time = toc(t2);
    time_so_far = time_so_far + curr_time;
    prct_complete = read_MB/total_MB;
    total_time = time_so_far/prct_complete;
    curr_rate = curr_MB/curr_time;
    avg_rate = read_MB/time_so_far;
    fprintf(', %.0f/%.0f MB, %5.2f%% %.1f curr MB/s, %.1f avg MB/s, %.1f/%.1f seconds\n', ...
        read_MB, total_MB, prct_complete*100,  curr_rate, avg_rate, time_so_far, total_time);
end
toc(t1)
