function [ out ] = datafile2mat( filename, overwrite_mat, output_directory )
% [ out ] = datafile2mat( filename, overwrite_mat [opt], output_directory [opt] )
%  Parses IGRA data in filename and saves the data to a .mat file in
%  output_directory. Additionally returns the data in out.
%
%  If a .mat file with a matching name is found on the path and
%  overwrite_mat is set to false, this function returns an empty vector and
%  does not modify the existing .mat file or create a new one.
%
% INPUTS
%          filename - filename for IGRA data text file to be parsed
%     overwrite_mat - flag for whether to overwrite or ignore if a matching
%                     .mat file is found to already exist
%                     Default: true
%  output_directory - directory to save .mat file
%                     Default: current working directory

% Copyright (c) 2017-present, Facebook, Inc.
% All rights reserved.
%
% This source code is licensed under the BSD-style license found in the
% LICENSE file in the root directory of this source tree. An additional grant
% of patent rights can be found in the PATENTS file in the same directory.


% Set default values
if ~exist('overwrite_mat', 'var')
    overwrite_mat = true;
end

if ~exist('output_directory', 'var')
    output_directory = pwd;
end

% Format name of .mat file
[~, file_no_path, ext] = fileparts(filename);
mat_filename = fullfile(output_directory, [file_no_path ext '.mat']);

% If use_mat_if_found is true, look for mat file and load from there if
% found
if exist(mat_filename, 'file') && ~overwrite_mat
%     out = load(mat_filename);
    out = [];
    return
end

% Get formatting definitions
defs = GCSAL.IGRA.format_definitions( );

% open the file
[fid, msg] = fopen(filename);

if fid == (-1)
    error(message('MATLAB:fileread:cannotOpenFile', filename, msg));
end

% Throw error on bad file
if fid == -1; error('Could not find file: %s', filename); end

% Read text file as uint8. Working directly in uint8 is more efficient for
% operations on large datasets. Also we are safe to assume that the IGRA
% data files contain only UTF-8 characters so all characters can be
% represented with 1 byte instead of the 2 bytes of a char
%
% Additionally the entire text file is read in one line for speed. This is
% significantly faster than a while loop with fgetl().
try
    % read file
    orig_txt = fread(fid,'char=>uint8');
catch exception
    % close file
    fclose(fid);
	throw(exception);
end

% close file
fclose(fid);

% Lines of text associate with header information begin with #
i_hash = find(orig_txt == uint8('#'));

% Form indices into orig_txt for the location of header text characters.
% The start of each row of header_txt is given by i_hash and the width of
% each row of header text is fixed
header_idx = bsxfun(@plus, 0:defs.header.row_width, i_hash);
header_txt = orig_txt(header_idx);

% Find non-header text by simply removing the header text from the original
% text array
no_header_txt = orig_txt;
no_header_txt(header_idx) = [];

% Now reshape the non-header text so that each row is a line of text with
% fixed width.
if mod(length(no_header_txt), defs.entries.row_width+2) == 0
    entries_txt = reshape(no_header_txt, defs.entries.row_width+2, [])';
else
    error('File length not expected. Check for interrupted data')
end

% Parse header and entries text
out.header =  txt2params(defs.header, header_txt);
out.entries = txt2params(defs.entries, entries_txt);

% Save to mat file. Use -v6 for faster loading (v6 option does not compress
% data so you get larger files but faster loading)
save(mat_filename, '-v6', '-struct', 'out')

end


function out = txt2params(def, text_mat)
% For each Param in def create a Param object with the character array
% text_mat.

flds = fieldnames(def.params);
for i = 1:length(flds)
    out.(flds{i}) = GCSAL.IGRA.Param(def.params.(flds{i}), text_mat);
end

end
