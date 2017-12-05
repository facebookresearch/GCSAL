%% Step 1 is to download the source .h5 file and set up your paths correctly

%clear;
close all; clc;

if (~exist('g','var') || ~isa(g,'GCSAL.GCSAL'))
    % The gcsal.h5 and gcsal.h5.info.mat files are available for download from the
    % website and should be placed in the h5_data directory.
    
    %%%%%%%%%%%%%% CHANGE THESE %%%%%%%%%%%%%%%%%%
    % Set this to wherever you put the gcsal.h5 file and gcsal.h5.info.mat
    % files downloaded from dropbox
    h5_dir = './h5_data/';
    
    % Directory to code. The folder +GCSAL which contains this file should be
    % in this directory
    codebase_dir = './';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Full path to .mat file with h5 info
    h5_file = fullfile(h5_dir, 'gcsal.h5');
    h5_mat_file = [h5_file '.info.mat'];
    
    % Set up Matlab path
    addpath(genpath(codebase_dir))
    
    %% Load GCSAL object from .mat file
    
    % This requires about 6 gb of RAM just to hold all of the header
    % information in memory
    
    % Normally you should load the GCSAL object from the .mat file but if it
    % doesn't exist on your path you can use the .h5 file. After using the .h5
    % file a .mat file will be created automatically for subsequent use
    if ~exist(h5_mat_file, 'file')
        g = GCSAL.GCSAL(h5_file);
    else
        g = GCSAL.GCSAL(h5_mat_file);
        g.h5_fname = h5_file;
    end
end

% Find 20 nearest stations to Capetown, South Africa

loc_lat = -33.938655;
loc_lon = 18.63863;
nsamples = 20;
stations = g.station_search('Nearest', [loc_lat loc_lon], 'Number', nsamples);

% Of the 20 nearest stations, grab the least number such that
% at each altitude between 18 and 25 km, there are at least 
% 300 wind samples per month

tooFew = true;
i = 0;

while tooFew
    i = i + 1;
    sttns = stations(1:i);
    lows = zeros(12,1);
    for k = 1:12
        [N, entries, stats] = g.counts2(sttns, 'gph', 'wspd', ...
            'FilterFields', {'month'}, 'FilterRanges', {k}, ...
            'Verbose', false, 'Plot', false);
        lows(k) = min(interp1(stats.x, sum(N,2), [18:25]));
    end
    if ((min(lows) >= 300) || (i == nsamples))
        tooFew = false;
    end
end

stations = stations(1:i);
fprintf('Need %d stations\n', length(stations));

figure;
plot([stations(:).arclen]/1e3,'-x');grid;
xlabel('Station Index');
ylabel('Station Distance [km]');
title('Station Distance from Capetown');

% Get wind speed/dir stats
[N1, entries1, stats1] = g.counts2(stations, 'gph', 'wspd', 'Verbose', false);
[N2, entries2, stats2] = g.counts2(stations, 'gph', 'wdir', 'Verbose', false);

% stats2.x is gph vector
% stats2.y is wdir vector

% At constant altitude
figure;plot(stats1.y,stats1.pdf(:,23));grid;
title('Yearly Wind Speed PDF At 22.5 km');
figure;plot(stats2.y,stats2.pdf(:,23));grid;
title('Yearly Wind Direction PDF At 22.5 km');

s = stats1.y;
p = stats1.pdf(:,23)';
wspds = datasample(s, 1e4, 'Weights', p, 'Replace', true);
figure;h1 = histogram(wspds,'BinMethod','sturges');grid;
title('Sampled Yearly Wind Speed PDF At 22.5 km');

s = stats2.y;
p = stats2.pdf(:,23)';
wdirs = datasample(s, 1e4, 'Weights', p, 'Replace', true);
figure;h2 = histogram(wdirs,'BinMethod','sturges');grid;
title('Sampled Yearly Wind Direction PDF At 22.5 km');

%% Reproduce distribution

% y1max = max(stats1.y);
% y1min = min(stats1.y);
% 
% y2max = max(stats2.y);
% y2min = min(stats2.y);
% 
% n = 1000;
% wind = zeros(n,1);
% dir = zeros(n,1);
% nsamples = 5;
% for j = 1:n
%     p = rand(nsamples,1);      % sample from uniform dist
%     y1 = p * (y1max - y1min) + y1min;   % Get wind values from sample
%     probs = interp1(stats1.y,stats1.pdf(:,23),y1);   % Probability of each wind
%     [~,idx] = sort(probs);  % Find max probability
%     wind(j) = y1(idx(end));
%     
%     p = rand(nsamples,1);      % sample from uniform dist
%     y2 = p * (y2max - y2min) + y2min;   % Get wind values from sample
%     probs = interp1(stats2.y,stats2.pdf(:,23),y2);   % Probability of each wind
%     [~,idx] = sort(probs);  % Find max probability
%     dir(j) = y2(idx(end));
% end
% figure;histogram(wind);figure;histogram(dir);
