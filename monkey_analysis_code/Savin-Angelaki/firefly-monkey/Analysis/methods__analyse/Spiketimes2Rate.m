function [nspk,timepoints] = Spiketimes2Rate(trials,timepoints,binwidth)

ntrls = length(trials);
timepoints = [timepoints(1)-binwidth timepoints timepoints(end)+binwidth];

multi_indx = any(arrayfun(@(x) any(size(x.tspk,2)>1), trials));
% force-column spikes times (needed for saccade case)
for i = 1:ntrls
    nsac(i) = size(trials(i).tspk,2); 
    trials_temp(i).tspk = trials(i).tspk(:);
end
% compute psth
nspk = hist(cell2mat({trials_temp.tspk}'),timepoints);
% nspk = histcounts(cell2mat({trials.tspk}'),timepoints); %Sina
if sum(arrayfun(@(x) isempty(x.tspk), trials_temp)) == numel(trials_temp); nspk(1:end) = nan; end
% throw away histogram edges
nspk = nspk(2:end-1); 
timepoints = timepoints(2:end-1);

% trial-average firing rates in units of spikes/s
if multi_indx == 1
nspk = nspk/(sum(nsac)*binwidth);    
else
nspk = nspk/(ntrls*binwidth);
end