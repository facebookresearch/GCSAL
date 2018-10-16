function out = fullpath(filepart1, varargin)
% Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved
%
% F = fullpath(filepart1, filepart2, ..., filepartN) builds a full
%     path specification F from the folders specified using / for the
%     fileseparator regardless of operating system. Input
%     arguments FOLDERNAME1, FOLDERNAME2, must be strings. The output
%     of fullfile is  equivalent to
%
%        F = [filepart1 / filepart2 / ... / filepartN]
%
%     except that care is taken to handle the cases when the folders begin
%     or end with a file separator.


% Error check for number of inputs
if length(varargin) < 1
    error('fullpath expects at least 2 inputs')
end

% Extract first input in varargin
filepart2 = varargin{1};

% Check that inputs are strings
if ~isa(filepart1, 'char') || ~isa(filepart2, 'char')
    error('fullpath expects string inputs')
end

% Remove / from end of filepart1 if it exists
if strcmp(filepart1(end), '/')
    filepart1(end) = [];
end

% Remove / from start of filepart2 if it exists
if strcmp(filepart2(1), '/')
    filepart2(1) = [];
end

% Concatenate fileparts with /
out = [filepart1 '/' filepart2];

% Recurse on remaining inputs
if length(varargin) > 1
    out = GCSAL.H5.fullpath(out, varargin{2:end});
end

end
