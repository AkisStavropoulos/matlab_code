function indx = get_inbetween_indices(startindx,stopindx)
% Get all indices between startindx and stopindx
% Input only indices

if length(startindx) > length(stopindx)
    startindx = startindx(1:end-1);
end
indx = [];
for i = 1:length(startindx)
    indx = [indx startindx(i):stopindx(i)];   
end