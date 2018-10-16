% Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved
%
% This example script takes the GCSAL process from downloading the source
% data from NOAA website through to creating a GCSAL Matlab object.
%
% The steps of this example script are only necessary if you want to
% download and re-build the GCSAL library from the original source data.
% Typically this is not necessary as you can use the provided
% .h5 file available on the website. Reasons you might want to do this:
%   1. Want to update the data with the latest measurements
%   2. Want to make a change to the way data is stored in the .h5 file


%% Set up paths and constants - Change the paths below as necessary
clear all; close all; clc;

%%%%%%%%%%%%%% CHANGE THESE %%%%%%%%%%%%%%%%%%%%%%%%
% Base directory. Text files should be in a directory called txt in this
% folder
base_dir = './raw_data/';

% Directory to save .h5 file
h5_dir = './h5_data/';

% Directory to code. The folder +GCSAL which contains this file should be
% in this directory
codebase_dir = './';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Directory containing IGRA .txt files
txt_dir = fullfile(base_dir, 'txt');

% filespec is used to find the source .txt files.
filespec = '*-data.txt'; % Use this to find process all .txt files
% filespec = 'BC*-data.txt'; % Use this for testing on a small subset files

% Directory to store .mat files
mat_dir = fullfile(base_dir, 'mat');

% Full path to .h5 file
h5_file = fullfile(h5_dir, 'gcsal.h5');

% Set up Matlab path
addpath(genpath(codebase_dir))

%% Step 1: Download IGRA data from NOAA

% Run "bash ./raw_data/download_igra.sh" on the command line from this directory
% This will take a while as you are downloading > 70 gb of data


%% Step 2: Convert IGRA txt file to .mat

% Set overwrite_mat true if you want to start from scratch and overwrite any
% existing .mat files that have already been made and exists on your path
overwrite_mat = true; %false;
GCSAL.IGRA.datafile2mat_dir(txt_dir, mat_dir, overwrite_mat, filespec);

%% Step 3: Convert .mat to .h5

% Set append_flag false if you want to start from scratch and clear the .h5
% file if it previously been made and exists on your path
append_flag = false; %true;
GCSAL.IGRA.mat2h5_dir(mat_dir, h5_file, append_flag, [filespec '.mat']);

%% Step 4: Load GCSAL object from h5 file

% The first time you create the object by pointing to the h5_file. This
% will create a .mat file which can be used after the first time
g = GCSAL.GCSAL(h5_file);


%%
% Now you can look in GCSAL.GCSAL_examples.m to learn about how to use the
% GCSAL object

% Get stations in Brazil
stations = g.station_search();

% Histogram for all windspeeds
t = tic;
[N, entries] = g.counts(stations, 'wspd');
toc(t)
