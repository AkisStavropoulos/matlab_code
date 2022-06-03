clc
close all
clear all
load('/Volumes/WD Edo/firefly_analysis/LFP_band/DATASET/PPC+PFC+MST/m53s114.mat')
continuous = cell2mat({trials_behv.continuous}); % extract continuous channels
events = cell2mat({trials_behv.events});
for i = 1:length(continuous)
    sprintf('%d/%d',i,length(continuous))
    indx_stop = find(continuous(i).ts > events(i).t_stop, 1);
    indx_beg = find(continuous(i).ts > events(i).t_targ, 1);
    x0_monk(i) = continuous(i).xmp(indx_beg+1); y0_monk(i) = continuous(i).ymp(indx_beg+1);
    x_fly(i) = nanmedian(continuous(i).xfp(indx_beg:indx_stop)); 
    y_fly(i) = nanmedian(continuous(i).yfp(indx_beg:indx_stop));
    
    xx = -x0_monk(i) + continuous(i).xmp; % position relative to monkey start
    yy = -y0_monk(i) + continuous(i).ymp;
    
    % local robust smooth estimate position (to calculate velocity angle)
    xx2 = smooth(0.006.*(1:length(xx)),xx,0.1,'rloess');
    yy2 = smooth(0.006.*(1:length(yy)),yy,0.1,'rloess');
    dxx = xx2(2:end) - xx2(1:end-1);
    dyy = yy2(2:end) - yy2(1:end-1);
    
    % angle of velocity vector 
    alpha = (atan2(dyy,dxx)) / (pi*2) * 360;
    
    % angle of straightline path from position to target
    lin_alpha = (atan2(x_fly(i) - xx2(1:end-1), y_fly(i) - yy2(1:end-1))) / (pi*2) * 360;
    
    continuous(i).alpha = nan(size(continuous(i).xmp));
    continuous(i).lin_alpha = nan(size(continuous(i).xmp));
    
    continuous(i).alpha(indx_beg:indx_stop) = alpha(indx_beg:indx_stop);
    continuous(i).lin_alpha(indx_beg:indx_stop) = lin_alpha(indx_beg:indx_stop);
    
    continuous(i).alpha(continuous(i).v < 15) = nan;
    
    
end
