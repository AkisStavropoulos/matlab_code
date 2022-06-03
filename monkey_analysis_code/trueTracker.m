function [idxOut,startTrue,endTrue] = trueTracker(arrayIn,sr,timeWin)
% This function finds the indices of logical trues sequencies longer than the
% timeWindow, sr is the sampling rate

arrayIn = logical(arrayIn);

startTrue = [];endTrue = [];
for i = 1:numel(arrayIn)
    if i == 1
        if arrayIn(i)
            startTrue(end+1) = i;
        end
    else
        
        if arrayIn(i) && ~arrayIn(i-1)
            startTrue(end+1) = i;
        elseif ~arrayIn(i) && arrayIn(i-1)
            endTrue(end+1) = i-1;
        end 
    end
end

if numel(startTrue) > numel(endTrue)
endTrue(end+1) = numel(arrayIn);
end

cleanIdx = (endTrue-startTrue)./sr < timeWin;
endTrue(cleanIdx) = [];
startTrue(cleanIdx) = [];

idxOut = [];
for i = 1:numel(startTrue)
idxOut(end+1:end+numel(startTrue(i):endTrue(i))) = startTrue(i):endTrue(i);
end
idxOut = sort(idxOut );
end