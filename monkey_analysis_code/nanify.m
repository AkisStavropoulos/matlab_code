function x = nanify(x,indx)

%% Set values to NaN
if nargin < 2
    indx = 1:numel(x);
end

x(indx) = nan;