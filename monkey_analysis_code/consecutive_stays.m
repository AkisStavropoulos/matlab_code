function [numstays,staysleft] = consecutive_stays(x)

indx = find(diff(x));
numstays = [indx numel(x)] - [0 indx];
cnt = arrayfun(@(X) X-1:-1:0, numstays , 'un',0);
staysleft = cat(2,cnt{:});


