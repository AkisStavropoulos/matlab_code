function data = ReSample2FixedDt(data,dt,t_beg1)

%% Resample data to a fixed sampling rate
N = 6; % fraction of dt
ts = dt:dt:data.trial_time(end)+dt; ts = ts(:);
% data.trial_time = data.trial_time - (data.trial_time(1)-t_beg1);
% II = find(diff(data.trial_time)==0)+1;
% data.trial_time(II) = data.trial_time(II)+dt/N;

reset = 1; cnt = 0;
while reset
    if numel(unique(data.trial_time))~=numel(data.trial_time)        

        indx0 = [0 ; diff(data.trial_time)==0];
        upindx = find(diff(indx0) > 0);
        dnindx = find(diff(indx0) < 0);
        
        if length(upindx)>length(dnindx); dnindx(end+1) = numel(indx0); elseif length(upindx)<length(dnindx); dnindx(end) = []; end
        cons_samp = dnindx - upindx;
        cons_dt = data.trial_time(dnindx+1) - data.trial_time(upindx);
        steps = cons_dt./(cons_samp+1);
        
        for n = 1:numel(upindx)
        data.trial_time(upindx(n)+1:dnindx(n)) = data.trial_time(upindx(n)) + [steps(n):steps(n):steps(n)*cons_samp(n)];    
        end
        cnt=cnt+1;
        if cnt>500; error('Can''t fix duplicate frames. Check file for saving issues.'); end
    else
        reset = 0;
    end
    if any(diff(data.trial_time) < 0)
        error('Timestamps are going back in time!');
    end
end

%% Interpolate NaNs in eye data

nanx = isnan(data.Gx); t1 = 1:numel(data.Gx); data.Gx(nanx) = interp1(t1(~nanx), data.Gx(~nanx), t1(nanx), 'pchip');
nanx = isnan(data.Gy); t1 = 1:numel(data.Gy); data.Gy(nanx) = interp1(t1(~nanx), data.Gy(~nanx), t1(nanx), 'pchip');
nanx = isnan(data.Gz); t1 = 1:numel(data.Gz); data.Gz(nanx) = interp1(t1(~nanx), data.Gz(~nanx), t1(nanx), 'pchip');

nanx = isnan(data.RXctr); t1 = 1:numel(data.RXctr); data.RXctr(nanx) = interp1(t1(~nanx), data.RXctr(~nanx), t1(nanx), 'pchip');
nanx = isnan(data.RYctr); t1 = 1:numel(data.RYctr); data.RYctr(nanx) = interp1(t1(~nanx), data.RYctr(~nanx), t1(nanx), 'pchip');
nanx = isnan(data.RZctr); t1 = 1:numel(data.RZctr); data.RZctr(nanx) = interp1(t1(~nanx), data.RZctr(~nanx), t1(nanx), 'pchip');
nanx = isnan(data.LXctr); t1 = 1:numel(data.LXctr); data.LXctr(nanx) = interp1(t1(~nanx), data.LXctr(~nanx), t1(nanx), 'pchip');
nanx = isnan(data.LYctr); t1 = 1:numel(data.LYctr); data.LYctr(nanx) = interp1(t1(~nanx), data.LYctr(~nanx), t1(nanx), 'pchip');
nanx = isnan(data.LZctr); t1 = 1:numel(data.LZctr); data.LZctr(nanx) = interp1(t1(~nanx), data.LZctr(~nanx), t1(nanx), 'pchip');

nanx = isnan(data.RXnorm); t1 = 1:numel(data.RXnorm); data.RXnorm(nanx) = interp1(t1(~nanx), data.RXnorm(~nanx), t1(nanx), 'pchip');
nanx = isnan(data.RYnorm); t1 = 1:numel(data.RYnorm); data.RYnorm(nanx) = interp1(t1(~nanx), data.RYnorm(~nanx), t1(nanx), 'pchip');
nanx = isnan(data.RZnorm); t1 = 1:numel(data.RZnorm); data.RZnorm(nanx) = interp1(t1(~nanx), data.RZnorm(~nanx), t1(nanx), 'pchip');
nanx = isnan(data.LXnorm); t1 = 1:numel(data.LXnorm); data.LXnorm(nanx) = interp1(t1(~nanx), data.LXnorm(~nanx), t1(nanx), 'pchip');
nanx = isnan(data.LYnorm); t1 = 1:numel(data.LYnorm); data.LYnorm(nanx) = interp1(t1(~nanx), data.LYnorm(~nanx), t1(nanx), 'pchip');
nanx = isnan(data.LZnorm); t1 = 1:numel(data.LZnorm); data.LZnorm(nanx) = interp1(t1(~nanx), data.LZnorm(~nanx), t1(nanx), 'pchip');

nanx = isnan(data.conv_dist); t1 = 1:numel(data.conv_dist); data.conv_dist(nanx) = interp1(t1(~nanx), data.conv_dist(~nanx), t1(nanx), 'pchip');

%% Resample all variables
if isnan(data.trial_time(1)); data.trial_time(1) = 0; end

fieldnames = fields(data);
temp = [];
for n = 1:numel(fieldnames)
    if ~any(strcmp(fieldnames{n},{'trial_time','trial_num','phase','on_off','FFx','FFy','FFz','FFvel','rewarded','start_trial','stop_trial','end_trial'}))
        data.(fieldnames{n}) = interp1(data.trial_time,data.(fieldnames{n}),ts);
    elseif any(strcmp(fieldnames{n},{'trial_num','phase','on_off','FFx','FFy','FFz','FFvel'}))
        data.(fieldnames{n}) = interp1(data.trial_time,data.(fieldnames{n}),ts,'nearest');
    end
end

data.trial_time = ts; % interp1(data.trial_time,data.trial_time,ts);


