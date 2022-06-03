function data = ReSample2FixedDt_human(data,dt,t_beg1)

%% Resample data to a fixed sampling rate
N = 6; % fraction of dt
ts = dt:dt:data.trial_time(end)+dt; ts = ts(:);
% data.trial_time = data.trial_time - (data.trial_time(1)-t_beg1);
% II = find(diff(data.trial_time)==0)+1;
% data.trial_time(II) = data.trial_time(II)+dt/N;

reset = 1;
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

nanx = isnan(data.Gx0); t1 = 1:numel(data.Gx0); data.Gx0(nanx) = interp1(t1(~nanx), data.Gx0(~nanx), t1(nanx), 'pchip');
nanx = isnan(data.Gy0); t1 = 1:numel(data.Gy0); data.Gy0(nanx) = interp1(t1(~nanx), data.Gy0(~nanx), t1(nanx), 'pchip');
nanx = isnan(data.Gz0); t1 = 1:numel(data.Gz0); data.Gz0(nanx) = interp1(t1(~nanx), data.Gz0(~nanx), t1(nanx), 'pchip');

nanx = isnan(data.hitX); t1 = 1:numel(data.hitX); data.hitX(nanx) = interp1(t1(~nanx), data.hitX(~nanx), t1(nanx), 'pchip');
nanx = isnan(data.hitY); t1 = 1:numel(data.hitY); data.hitY(nanx) = interp1(t1(~nanx), data.hitY(~nanx), t1(nanx), 'pchip');
nanx = isnan(data.hitZ); t1 = 1:numel(data.hitZ); data.hitZ(nanx) = interp1(t1(~nanx), data.hitZ(~nanx), t1(nanx), 'pchip');

nanx = isnan(data.LPD); t1 = 1:numel(data.LPD); data.LPD(nanx) = interp1(t1(~nanx), data.LPD(~nanx), t1(nanx), 'pchip');
nanx = isnan(data.RPD); t1 = 1:numel(data.RPD); data.RPD(nanx) = interp1(t1(~nanx), data.RPD(~nanx), t1(nanx), 'pchip');
nanx = isnan(data.Lopen); t1 = 1:numel(data.Lopen); data.Lopen(nanx) = interp1(t1(~nanx), data.Lopen(~nanx), t1(nanx), 'pchip');
nanx = isnan(data.Ropen); t1 = 1:numel(data.Ropen); data.Ropen(nanx) = interp1(t1(~nanx), data.Ropen(~nanx), t1(nanx), 'pchip');

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


