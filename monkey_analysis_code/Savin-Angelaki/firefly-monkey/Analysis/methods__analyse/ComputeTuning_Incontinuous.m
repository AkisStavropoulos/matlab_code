function tuningstats = ComputeTuning_Incontinuous(x,t_event,t_spk,timewindow,perievent_t,duration_zeropad,corr_lag,nbootstraps,tuning_prs,tuning_method,tuning_binrange)
% This function uses the average firing rate around the event, specified by perievent_t

ntrls = length(x);
if ntrls < nbootstraps % not enough trials
    tuningstats = [];
    return;
end
if nargin<10, tuning_binrange = []; end

%% concatenate data from different trials
[xt,~,yt,xt_pad,~,yt_pad] = ConcatenateTrials_Incontinuous(x,[],t_spk,t_event,timewindow);

%% calculate binwidth and nanify corss-correlation (not-applicable since we use mean firing rates)
ind = find(cellfun(@(x) numel(x),t_event),1); if isempty(ind); ind=1; end
if ~(size(t_spk{ind},2) == length(t_event{ind})) % data as  spike times
    temporal_binwidth = perievent_t;
    disp('WARNING: Number of events are not correct (ComputeTuning_Incontinuous.m). Please check!!');
else % any other data
    temporal_binwidth = perievent_t;
end
tuningstats.xcorr.val = nan;
tuningstats.xcorr.lag = nan;

%% compute tuning curves
if strcmp(tuning_method,'binning')
    nbins = tuning_prs.nbins1d_binning; % load predefined number of bins
    [tuningstats.tuning.stim,tuningstats.tuning.rate,tuningstats.tuning.pval] = NPregress_binning(xt,yt,temporal_binwidth,nbins,nbootstraps,tuning_binrange);
elseif strcmp(tuning_method,'k-nearest')
    k = arrayfun(tuning_prs.k_knn,numel(xt)); % compute k from predefined anonymous function
    nbins = tuning_prs.nbins1d_binning; % load predefined number of bins
    [tuningstats.tuning.stim,tuningstats.tuning.rate] = NPregress_knn(xt,yt,temporal_binwidth,k,nbins,nbootstraps);
elseif strcmp(tuning_method,'nadaraya-watson')
    kernel = tuning_prs.kernel_nw; % load kernel for smoothing
    bandwidth = tuning_prs.bandwidth_nw; % load kernel for smoothing
    [tuningstats.tuning.stim,tuningstats.tuning.rate] = NPregress_nw(xt,yt,temporal_binwidth,kernel,bandwidth,[],nbootstraps);
elseif strcmp(tuning_method,'local-linear')
    kernel = tuning_prs.kernel_locallinear; % load kernel for smoothing
    bandwidth = tuning_prs.bandwidth_locallinear; % load kernel for smoothing
    [tuningstats.tuning.stim,tuningstats.tuning.rate] = NPregress_locallinear(xt,yt,temporal_binwidth,kernel,bandwidth,[],nbootstraps);
end