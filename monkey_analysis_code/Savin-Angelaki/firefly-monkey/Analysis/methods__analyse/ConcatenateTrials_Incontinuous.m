function [xt,zt,yt,xt_pad,zt_pad,yt_pad] = ConcatenateTrials_Incontinuous(x,z,t_spk,t_events,timewindow)

% x, ts, and tspk are cell arrays of length N
% x{i}: time-series of stimulus in trial i
% z{i}: time of event in trial i
% tspk{i}: vector of spike times in trial i
% ts{i}: vector of time points in the ith trial
% timewindow: Nx2 array - the columns corresponds to start and end of analysis window
% e.g. to analyse all datapoints, timeindow(i,:) = [ts{i}(1) ts{i}(end)]

ntrls = length(t_events);
nevents = sum(cellfun(@(x) numel(x), t_events));
twin = mat2cell(timewindow,ones(1,ntrls));

%% concatenate data from different trials
% concatenate spikes
ind = find(cellfun(@(x) numel(x),t_events),1); if isempty(ind); ind=1; end
if ~(size(t_spk{ind},2) == length(t_events{ind})) % data as  spike times
    y = cellfun(@(x,y) hist(x,y),t_spk,t_events,'UniformOutput',false);
    if sum(cellfun(@(x) isempty(x),t_events))==numel(t_events); y = t_spk; end % fix for the case of no eye movements recorded
else % data already in spike counts
    y = t_spk;
end
t2 = cellfun(@(x) x(:),t_events,'UniformOutput',false);
y2 = cellfun(@(x) sum(x,1),y,'UniformOutput',false);
y2 = cellfun(@(x) x(:),y2,'UniformOutput',false); % transpose is to reshape to column vector
yt = cellfun(@(x,t,z) x((t>z(1) & t<z(2))),y2(:),t2(:),twin(:),'UniformOutput',false); % only keep events within t_targ and t_stop

% concatenate stimulus
xt = [];
if ~isempty(x)
    x2 = cellfun(@(x) x(:),x,'UniformOutput',false);
    xt = cellfun(@(x,t,z) x((t>z(1) & t<z(2))),x2(:),t2(:),twin(:),'UniformOutput',false); % only keep events within t_targ and t_stop
end

% concatenate events
zt = [];
if ~isempty(z)
    z2 = cellfun(@(x,y) [diff(y>x) ; 0],z(:),t2(:),'UniformOutput',false); % transpose is to reshape to column vector
    zt = cellfun(@(x,t,z) x(t>z(1) & t<z(2)),z2(:),t2(:),twin(:),'UniformOutput',false);
end

xt_pad = [];
zt_pad = [];
yt_pad = [];

if ~isempty(xt), xt = cell2mat_singleDouble(xt); end
if ~isempty(zt), zt = cell2mat_singleDouble(zt); end
yt = cell2mat_singleDouble(yt);