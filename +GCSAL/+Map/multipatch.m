function patch_handles = multipatch( x, y, varargin )
% Copyright (c) Facebook, Inc. and its affiliates.
%
% This source code is licensed under the MIT license found in the
% LICENSE file in the root directory of this source tree.
%
% patch_handles = multipatch( x, y, varargin )
%   Like built in matlab function patch but allows for NaN values in x and
%   y to separate multiple patches. varargin allows for any additional
%   inputs to the patch function.
%
%   returns patch_handles, a vector of handles returned by each call to
%   patch()


if ~isequal(size(x), size(y))
    error('x and y must be same size')
end

% Find nan indices
xnan = isnan(x);
ynan = isnan(y);

% Find where either x or y is nan
anynan = find(xnan | ynan);

% anynan will be used for start/stop indices with any index in anynan being
% skipped over. To make the for loop smooth, add indices 0 and length+1 to
% anynan
anynan = [0 anynan length(x)+1];

% initialize patch_handles
patch_handles = [];

% Loop through anynan
for i = 2:length(anynan)

    % Choose idx between previous and next nan value
    idx = anynan(i-1)+1:anynan(i)-1;

    % Check that idx is not empty
    if ~isempty(idx)
        % Create a new patch with values at idx
        patch_handles(end+1) = patch(x(idx), y(idx), varargin{:}); %#ok<AGROW>
    end
end

end
