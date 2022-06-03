function x = zerofy(x,indx)

%% Set values to Zero
if nargin < 2
    indx = 1:numel(x);
end

x(indx) = 0;