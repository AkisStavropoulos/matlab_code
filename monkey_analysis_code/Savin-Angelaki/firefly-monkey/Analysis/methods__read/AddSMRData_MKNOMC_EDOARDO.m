function [trl,ch,t] = AddSMRData_MKNOMC_EDOARDO(data,cal,gains,prs)
% this version filters the data after dividing into trials
% data = ImportSMR('kl142.smr'); default_prs;
%% check channel headers (done)
plt = 0;
nch = length(data);
ch_title = cell(1,nch);
hdr = {data.hdr};
for i=1:nch
    if ~isempty(hdr{i})
        ch_title{i} = hdr{i}.title;
    else
        ch_title{i} = 'nan';
    end
end

%% channel titles (done)
chno.mrk = find(strcmp(ch_title,'marker')); 
% chno.kbd = find(strcmp(ch_title,'Keyboard')); 

chno.FFdraw = find(strcmp(ch_title,'FFDraw')); % trial start indicator 

chno.yle = find(strcmp(ch_title,'LDy')); 
chno.zle = find(strcmp(ch_title,'LDz'));
chno.yre = find(strcmp(ch_title,'RDy')); 
chno.zre = find(strcmp(ch_title,'RDz'));

chno.phi = find(strcmp(ch_title,'MonkeyYa')); % monkey yaw position
chno.xmp = find(strcmp(ch_title,'MonkeyX')); 
chno.ymp = find(strcmp(ch_title,'MonkeyY'));
chno.xfp = find(strcmp(ch_title,'FireflyX')); 
chno.yfp = find(strcmp(ch_title,'FireflyY'));
chno.v = find(strcmp(ch_title,'ForwardV')); 
chno.w = find(strcmp(ch_title,'AngularV'));
% chno.mtr = find(strcmp(ch_title,'MotorCom')); % yaw motor command
% chno.xac = find(strcmp(ch_title,'AccX')); % acceleration x-axis
% chno.yac = find(strcmp(ch_title,'AccY')); % acceleration y-axis
% 
% chno.vrol = find(strcmp(ch_title,'VelRoll')); % roll velocity 
% chno.vyaw = find(strcmp(ch_title,'VelYaw')); % yaw velocity 
% chno.vpit = find(strcmp(ch_title,'VelPitch')); % pitch velocity 


%% scale (done)
scaling.t = data(chno.mrk).hdr.tim.Scale*data(chno.mrk).hdr.tim.Units; % for markers
if isfield(chno, 'FFdraw'), scaling.FFdraw = data(chno.FFdraw).hdr.adc.Scale; offset.FFdraw = data(chno.FFdraw).hdr.adc.DC; end
try
scaling.xfp = data(chno.xfp).hdr.adc.Scale; offset.xfp = data(chno.xfp).hdr.adc.DC;
scaling.yfp = -data(chno.yfp).hdr.adc.Scale; offset.yfp = -data(chno.yfp).hdr.adc.DC;
catch
end
scaling.yle = data(chno.yle).hdr.adc.Scale; offset.yle = data(chno.yle).hdr.adc.DC;
scaling.yre = data(chno.yre).hdr.adc.Scale; offset.yre = data(chno.yre).hdr.adc.DC; 
scaling.zle = data(chno.zle).hdr.adc.Scale; offset.zle = data(chno.zle).hdr.adc.DC; 
scaling.zre = data(chno.zre).hdr.adc.Scale; offset.zre = data(chno.zre).hdr.adc.DC;

scaling.phi = data(chno.phi).hdr.adc.Scale; offset.phi = data(chno.phi).hdr.adc.DC;
scaling.xmp = data(chno.xmp).hdr.adc.Scale; offset.xmp = data(chno.xmp).hdr.adc.DC;
scaling.ymp = -data(chno.ymp).hdr.adc.Scale; offset.ymp = -data(chno.ymp).hdr.adc.DC; % !!!!!SIGN
scaling.v = data(chno.v).hdr.adc.Scale; offset.v = data(chno.v).hdr.adc.DC;
scaling.w = data(chno.w).hdr.adc.Scale; offset.w = data(chno.w).hdr.adc.DC;
% scaling.mtr = data(chno.mtr).hdr.adc.Scale; offset.mtr = data(chno.mtr).hdr.adc.DC;
% scaling.xac = data(chno.xac).hdr.adc.Scale; offset.xac = data(chno.xac).hdr.adc.DC;
% scaling.yac = data(chno.yac).hdr.adc.Scale; offset.yac = data(chno.yac).hdr.adc.DC;
% 
% scaling.vrol = data(chno.vrol).hdr.adc.Scale; offset.vrol = data(chno.vrol).hdr.adc.DC;
% scaling.vyaw = data(chno.vyaw).hdr.adc.Scale; offset.vyaw = data(chno.vyaw).hdr.adc.DC;
% scaling.vpit = data(chno.vpit).hdr.adc.Scale; offset.vpit = data(chno.vpit).hdr.adc.DC;
if isfield(chno,'microstim'), scaling.microstim = data(chno.microstim).hdr.adc.Scale; offset.microstim = data(chno.microstim).hdr.adc.DC; end
%% load relevant channels
chnames = fieldnames(chno); MAX_LENGTH = inf; dt = []; maxlength = [];
for i=1:length(chnames)
    if ~any(strcmp(chnames{i},'mrk'))
        try
        ch.(chnames{i}) = double(data(chno.(chnames{i})).imp.adc)*scaling.(chnames{i}) + offset.(chnames{i});
        dt = [dt prod(data(chno.(chnames{i})).hdr.adc.SampleInterval)];
        MAX_LENGTH = min(length(ch.(chnames{i})),MAX_LENGTH);
        maxlength = [maxlength length(ch.(chnames{i}))];
        catch
        end
    end
end

ch.('mrk') = double(data(chno.('mrk')).imp.tim)*scaling.t;

if length(unique(dt))==1
    dt = dt(1);
else
   error('channels must all have identical sampling rates');
end
% % fix length of data
% for i=1:length(chnames)
%     if ~any(strcmp(chnames{i},'mrk'))
%         ch.(chnames{i}) = ch.(chnames{i})(1:MAX_LENGTH); 
%     end
% end

%% define filter
sig = prs.filtwidth; %filter width
sz = prs.filtsize; %filter size
t2 = linspace(-sz/2, sz/2, sz);
h = exp(-t2.^2/(2*sig^2));
h = h/sum(h); % normalise filter to ensure area under the graph of the data is not altered
%% filter position and velocity channels
for i=1:length(chnames)
    if ~any(strcmp(chnames{i},{'mrk','yle','yre','zle','zre','microstim', 'FFdraw'}))
        try
        ch.(chnames{i}) = conv(ch.(chnames{i})(1:MAX_LENGTH),h,'same');
%         ch.(chnames{i}) = ch.(chnames{i})(sz/2+1:end);
        catch
            if ~any(strcmp(chnames{i},{'mrk','yle','yre','zle','zre','microstim', 'xfp', 'yfp', 'FFdraw'}))
                ch.(chnames{i}) = conv(ch.(chnames{i})(1:MAX_LENGTH),h,'same');
            end
        end
    end
end
ch.yle = ch.yle(1:MAX_LENGTH);
ch.yre = ch.yre(1:MAX_LENGTH);
ch.zle = ch.zle(1:MAX_LENGTH);
ch.zre = ch.zre(1:MAX_LENGTH);
try
ch.FFdraw = ch.FFdraw(1:MAX_LENGTH);
catch
end
% we are going to add xfp and yfp if they do not exist, to make the rest of
% the code run as originally intended
if ~isfield(ch, 'xfp')
    ch.xfp(size(ch.xmp, 1), 1) = 0;
    ch.yfp(size(ch.xmp, 1), 1) = 0;
end


ts = dt:dt:length(ch.(chnames{end-1}))*dt;
if any(strcmp(chnames,'microstim'))
    ts2 = dt_microstim:dt_microstim:length(ch.microstim)*dt_microstim; ch.microstim = interp1(ts2,ch.microstim,ts); 
    ch.microstim = ch.microstim(:);
end

%% import hand position
[ch.h1, ch.h2, isavailable] = ImportHandPosition(ch.v,ch.w,dt,prs);
if isavailable, chnames(end+1:end+2) = {'h1','h2'}; end % append two more channels

%% event markers (done)
markers = data(chno.mrk).imp.mrk(:,1); % from 'keyboard'
%% event times (done)
ts = dt:dt:dt*MAX_LENGTH;
ch.ts = ts'; % in seconds
t.events = double(data(chno.mrk).imp.tim)*scaling.t;
% 115 = s (start), 101 = e (end), 106 = j (juice), 117 = u (shutter close) , 118 = v (vestibular only) , 114 = r (end of block)
t.fstart = t.events(markers==1);
t.fstop = ts(end);
t.beg = t.events(markers ==2); % t.events(markers ==115); 
t.end = t.events(markers == 3); % t.events(markers == 101); % first marker is only start
t.beg = t.beg(1:length(t.end));
t.reward = t.events(markers==4); % t.reward = t.events(markers == 104);
% t.shutter = t.events(markers == 117);

% previous marking using FFdraw
jitter = prs.jitter_marker;
FFdraw = ch.FFdraw;
FFdraw(FFdraw>=.2)= 5;
FFdraw(FFdraw< .2)= 0;
pulse_max = max(FFdraw);
ch.FFdraw = FFdraw;

% make sure that last element of FFdraw is 0
% FFdraw(end) = 0;
   
% % sanity check
% figure;plot(ts,FFdraw);vline(t.beg);hold on;plot(ts,.01*ch.xmp);plot(ts,.01*ch.ymp);
% check that the FFdraw starts at 0 and not 5
count = 1;
while FFdraw(count) == 5
    FFdraw(count) = 0;
    count = count+1;
end
if count > 1
    disp('FFdraw starts at 5 for this block.')
end

% % use markers to choose events based on FFdraw (change markers)
% pulseOFF = find(diff(FFdraw) == -pulse_max);
% pulseON = find(diff(FFdraw) == pulse_max);
% 
% for i = 1:length(t.beg)       
% %         if i ~= length(t.beg)               % keep the last marker intact, comes from another marker(114)
%             indx_end = find(pulseON >= t.beg(i)/dt, 1);
%             if isempty(indx_end)
%                 t.beg(end) = [];
%                 t.end(end) = [];
%                 
%                 % keep track 
%                 trial_terminated = i-1;
%                 spike2terminated = ['Trials cut at ' num2str(trial_terminated) ', probably terminated Spike2 !'];
%                 disp(spike2terminated)
%             else
%                 t.end(i) = pulseON(indx_end)*dt;
%             end
% %         end
% end

%% refine t.beg to ensure it corresponds to target onset
% jitter = prs.jitter_marker;
% dPm__dt = [0 ; sqrt(diff(ch.ymp).^2 + diff(ch.xmp).^2)]; % derivative of monkey position
% [~,t_teleport] = findpeaks(dPm__dt,dt*(1:length(dPm__dt)),'MinPeakHeight',prs.minpeakprominence.monkpos); % detect peaks
% dPf__dt = [0 ; diff(ch.FFdraw)]; % derivative of firefly position
% [~,t_flyON] = findpeaks(dPf__dt,dt*(1:length(dPf__dt)),'MinPeakHeight',prs.minpeakprominence.flypos); % detect peaks
% t_teleport_trl = nan(length(t.beg),1); t_flyON_trl = nan(length(t.beg),1);
% for i=1:length(t.beg)
%     t_teleport_temp = t_teleport(t_teleport > (t.beg(i) - jitter) &  t_teleport < (t.beg(i) + jitter));
%     if ~isempty(t_teleport_temp), t_teleport_trl(i) = t_teleport_temp(end); end
%     t_flyON_temp = t_flyON(t_flyON > (t.beg(i) - jitter) &  (t_flyON < t.beg(i) + jitter));
%     if ~isempty(t_flyON_temp), t_flyON_trl(i) = t_flyON_temp(end); end
% end
% tflyON_minus_teleport = nanmedian(t_flyON_trl - t_teleport_trl);
% % set trial begin time equal to target onset
% t_beg_original = t.beg;
% for i=1:length(t.beg)
%     if ~isnan(t_flyON_trl(i)), t.beg(i) = t_flyON_trl(i);
%     elseif ~isnan(t_teleport_trl(i)), t.beg(i) = t_teleport_trl(i) + tflyON_minus_teleport;
%     end
% end
% t_beg_correction = t.beg - t_beg_original;

%% refine t.beg to ensure it corresponds to target onset
jitter = prs.jitter_marker;
dPm__dt = [0 ; sqrt(diff(ch.ymp).^2 + diff(ch.xmp).^2)]; % derivative of monkey position
[~,t_teleport] = findpeaks(dPm__dt,dt*(1:length(dPm__dt)),'MinPeakHeight',prs.minpeakprominence.monkpos); % detect peaks
if isfield(ch, 'FFdraw')
    t_flyON = dt*find(diff(ch.FFdraw)>3); % hope you enjoy the hardcode here :)
    %t_flyON(end+1) = t.beg(end); 
else
    dPf__dt = [0 ; sqrt(diff(ch.yfp).^2 + diff(ch.xfp).^2)]; % derivative of firefly position
    [~,t_flyON] = findpeaks(dPf__dt,dt*(1:length(dPf__dt)),'MinPeakHeight',prs.minpeakprominence.flypos); % detect peaks
end

t_teleport_trl = nan(length(t.beg),1); t_flyON_trl = nan(length(t.beg),1);
for i=1:length(t.beg)
    t_teleport_temp = t_teleport(t_teleport > (t.beg(i) - jitter) &  t_teleport < (t.beg(i) + jitter));
    if ~isempty(t_teleport_temp), t_teleport_trl(i) = t_teleport_temp(end); end
    t_flyON_temp = t_flyON(t_flyON > (t.beg(i) - jitter) &  (t_flyON < t.beg(i) + jitter));
    if ~isempty(t_flyON_temp), t_flyON_trl(i) = t_flyON_temp(end); end
end
tflyON_minus_teleport = nanmedian(t_flyON_trl - t_teleport_trl);
% set trial begin time equal to target onset
t_beg_original = t.beg;
for i=1:length(t.beg)
    if ~isnan(t_flyON_trl(i)), t.beg(i) = t_flyON_trl(i);
    elseif ~isnan(t_teleport_trl(i)), t.beg(i) = t_teleport_trl(i) + tflyON_minus_teleport;
    end
end
t_beg_correction = t.beg - t_beg_original;

%% keep this version of event indication for future use
% for i = 1:length(t.beg)
%     
%     % check  whether all markers are included in a pulse
%     indstop = find(pulseON < t.beg(i)/dt,1,'last');
%     indstart = find(pulseOFF > t.beg(i)/dt,1);
%     dif = pulseOFF(indstart) - pulseON(indstop);
%     % check whether marker is included in a pulse
%     if ~(dif>0)
%         fprintf(['marker number ' num2str(i) ' not within a pulse!!! Check AddSMRData.m\n'])
%     end
%     
%     %     while dif > 43
%     %         fprintf(['marker number ' num2str(i) ' not within a pulse!!! Check AddSMRData.m\n'])
%     %         % fix
%     %         indstartup = indstart + 1;
%     %         difup = pulseOFF(indstartup) - pulseON(indstop);
%     %
%     %         indstartdown = indstart - 1;
%     %         difdown = pulseOFF(indstartdown) - pulseON(indstop);
%     %
%     %         if difup <= 43
%     %             indstart = indstartup;
%     %         elseif difdown <= 43
%     %             indstart = indstartdown;
%     %         end
%     %
%     %         dif = pulseOFF(indstart) - pulseON(indstop);
%     %         if dif <= 43
%     %             fprintf(['marker number ' num2str(i) ' fixed.\n'])
%     %         else
%     %             fprintf(['marker number ' num2str(i) ' STILL NOT fixed.\n'])
%     %         end
%     %
%     %     end
%     
%     
%     % now change marker
%     indx_beg = find(pulseOFF >= t.beg(i)/dt, 1);
%     if isempty(indx_beg)
%         t.beg(i) = nan;
%         t.end(i) = nan;
%     else
%         t.beg(i) = pulseOFF(indx_beg)*dt;
%         
%         if i ~= length(t.beg)               % keep the last marker intact, comes from another marker(114)
%             indx_end = find(pulseON >= t.beg(i)/dt, 1);
%             t.end(i) = pulseON(indx_end)*dt;
%         end
%         
%     end
%     
%     
% end
% 
% % fix for wrong markers, added after the last additions (JSonly condition, moog neutralization after trial, etc..)
% t.beg(isnan(t.beg)) = [];
% t.end(isnan(t.end)) = [];
% % sanity check
% % figure;plot(ts,FFdraw);vline(t.beg);vline(t.end,'k');hold on;plot(ts,.01*ch.xmp);plot(ts,.01*ch.ymp);
%%
% t.events = sort([t.beg ; t.end ; t.shutter ; t.reward]);
% ch.FFdraw = FFdraw;
% 
% % now find events
% if find(markers == 104) % if FEEDBACK
%
%     indx_beg = find(diff(FFdraw) == -pulse_max); % choose the end of the pulse, since it comes after the reset of position
%     indx_end = find(diff(FFdraw) == pulse_max); % choose the start of the pulse, since it is the command to reset position
%     t.beg = ts(indx_beg(1:2:end));
%     t.end = ts(indx_end(1:2:end));
%     t.end(1) = [];
%     t.reward = ts(indx_beg(2:2:end));
%     t.beg = t.beg(1:length(t.end));
%     t.events = sort([ts(indx_beg) ts(indx_end)]);
% else
%     indx_beg = find(diff(FFdraw) == -pulse_max); % choose the end of the pulse, since it comes after the reset of position
%     indx_end = find(diff(FFdraw) == pulse_max); % choose the start of the pulse, since it is the command to reset position
%
%     t.events = sort([ts(indx_beg) ts(indx_end)]);
%     t.events(1) = [];     t.events(end) = []; % see FFdraw plot to check events, 1st event start, last event end
%
%     t.beg = t.events(1:2:end);
%     t.end = t.events(2:2:end);
%     t.beg = t.beg(1:length(t.end));
% end
% % indxon = find(round(diff(ch.FFdraw)) == pulse_max);
% % indxoff = find(round(diff(ch.FFdraw)) == -pulse_max);
% % pulse_duration  = indxoff - indxon;
%% Eye Calibration
%% Re-calibrate using EyeCal file, output offset and scale

if prs.recalibrate && ~isempty(gains) && ~isempty(cal)
    [G,NOREC] = ReCalibrateEyes(cal,gains,prs);
else
    G = []; NOREC = 0;
end

% determine whether there is eye recording
X = [ch.zle ch.zre ch.yle ch.yre];  
% Xstd = std(X(1:round(size(X,1)/3),:),1);
if all(prs.eyechannels == 0) % || NOREC == 1
    NOREC = 1;
    ch.zle(1:end) = 0;    ch.yle(1:end) = 0;
    ch.zre(1:end) = 0;    ch.yre(1:end) = 0;
    nanindx = zeros(length(ch.zle),1);
    blankindx = 1:length(ch.zle);
else
    NOREC = 0;
    %% Remove blinks and smooth
    % subtract offset to align blinks and detect shifts clearly, then replace
   
    if ~isempty(G)
    Xoffset = [G{2}(2) G{4}(2) G{1}(2) G{3}(2)];
    Xscale = [G{2}(1) G{4}(1) G{1}(1) G{3}(1)];
    X = X - Xoffset;
    end
    
    if 0;  figure;plot(X);grid on;ylim([-100 100]);  end
    
    [~,nanindx] = ReplaceWithNans(X, prs.blink_thresh_pos,prs.nanpadding); % patch nans of 60 ms around blinks
    X(nanindx,:) = nan;
    % interpolate
    X(end+1,:) = zeros(1,4);
    t1 = 1:numel(X(:,1));   X(nanindx,1) = interp1(t1(~nanindx), X(~nanindx,1), t1(nanindx), 'pchip');
    t1 = 1:numel(X(:,2));   X(nanindx,2) = interp1(t1(~nanindx), X(~nanindx,2), t1(nanindx), 'pchip');
    t1 = 1:numel(X(:,3));   X(nanindx,3) = interp1(t1(~nanindx), X(~nanindx,3), t1(nanindx), 'pchip');
    t1 = 1:numel(X(:,4));   X(nanindx,4) = interp1(t1(~nanindx), X(~nanindx,4), t1(nanindx), 'pchip');
    X = X(1:end-1,:);

    crazy_indx = find(abs(X) > 120);
    [I,J] = ind2sub(size(X),crazy_indx);
    X(I,J) = 0;

    % Remove one-sample spikes with median filter
    X = medfilt1(X,prs.medfiltwidth);
    
    % Savitzky-Golay filter
    order = 2;  frames = prs.sgfiltwidth;
    X = sgolayfilt(X,order,frames);
    
    
%     % low-pass filter out noise at 14 Hz
%     fc = 14; % Hz
%     fs = prs.fs_smr;
%     [b,a] = butter(6,fc/(fs/2));
%     X = filtfilt(b,a,X);
%     X(nanindx,:) = nan;

%     [~,nanindx1] = ReplaceWithNans(X-Xavg, prs.blink_thresh,prs.nanpadding,'normal');
%     [~,nanindx2] = ReplaceWithNans(X-Xavg, 10,prs.nanpadding,'derivative');
%     [~,nanindx3] = ReplaceWithNans([zeros(1,size(X,2)) ; diff(X-Xavg)], 10,prs.nanpadding,'derivative')
%     nanindx = nanindx1|nanindx2|nanindx3;
%     X(nanindx,:) = nan;
%     % remove shifts to re extract blinks
%     filtwidth = 60/dt;
%     Xavg = movmean(X,filtwidth,1,'omitnan');
%     [~,nanindx] = ReplaceWithNans(X-Xavg, prs.blink_thresh,prs.nanpadding,'normal');
%     X(nanindx,:) = nan;

    % re-scale if applicable
    if ~isempty(G)
        X = X./Xscale; % offset is already subtracted
        newscale = [G{2}(1) G{4}(1) G{1}(1) G{3}(1)];
        newoffset = [G{2}(2) G{4}(2) G{1}(2) G{3}(2)];
        X = X.*newscale + newoffset;
    end
    % create moving window of 1 min to detect shifts
    filtwidth = round(60/dt);
    Xavg = movmean(X,filtwidth,1,'omitnan');
    Xmovstd = movstd(X,filtwidth,1,'omitnan');
    if size(X,1) > 2*filtwidth % make sure block is long enough
        Xavgoffset = nanmean(Xavg(1:2*filtwidth,:),1);
    else
        Xavgoffset = nanmean(Xavg(1:end,:),1);
    end
    if 0;  figure;plot(Xavg);grid on;hline(Xavgoffset),vline(2*filtwidth);  end
        
    % subtract shifts
    X = X - Xavg;
    X = X + Xavgoffset;
    % interpolate
    nanperc = sum(isnan(X),1)/size(X,1);
    X(end+1,:) = zeros(1,4);
    t1 = 1:numel(X(:,1));   X(nanindx,1) = interp1(t1(~nanindx), X(~nanindx,1), t1(nanindx), 'pchip');
    t1 = 1:numel(X(:,2));   X(nanindx,2) = interp1(t1(~nanindx), X(~nanindx,2), t1(nanindx), 'pchip');
    t1 = 1:numel(X(:,3));   X(nanindx,3) = interp1(t1(~nanindx), X(~nanindx,3), t1(nanindx), 'pchip');
    t1 = 1:numel(X(:,4));   X(nanindx,4) = interp1(t1(~nanindx), X(~nanindx,4), t1(nanindx), 'pchip');
    X = X(1:end-1,:);
    crazy_indx = find(abs(X) > 120);
    [I,J] = ind2sub(size(X),crazy_indx);
    X(I,J) = 0;
    % Filter eye signal
%     sig = 5*prs.filtwidth; %filter width
%     sz = 10*prs.filtsize; %filter size
%     t2 = linspace(-sz/2, sz/2, sz);
%     h = exp(-t2.^2/(2*sig^2));
%     h = h/sum(h); % normalise filter to ensure area under the graph of the data is not altered
%     
%     X(1,:) = conv(X(1,:),h,'same'); X(2,:) = conv(X(2,:),h,'same');
%     X(3,:) = conv(X(3,:),h,'same'); X(4,:) = conv(X(4,:),h,'same');
    X([1:500 end-500:end],:) = 0; % avoid artifacts   
    
    ch.zle = X(:,1); ch.zre = X(:,2); ch.yle = X(:,3); ch.yre = X(:,4);
    if 0;  figure;plot([ch.zle ch.yle ch.zre ch.yre]);grid on;ylim([-100 100]);   end
    
end
blinks = nanindx;
%% Detect saccades
% take derivative of eye position = eye velocity
if all(prs.eyechannels ~= 0)
    dze = diff(0.5*(ch.zle + ch.zre));
    dye = diff(0.5*(ch.yle + ch.yre));
elseif prs.eyechannels(1) ~= 0
    dze = diff(ch.zle);
    dye = diff(ch.yle);
else
    dze = diff(ch.zre);
    dye = diff(ch.yre);
end
sig = 10*prs.filtwidth; %filter width
sz = 10*prs.filtsize; %filter size
t2 = linspace(-sz/2, sz/2, sz);
h = exp(-t2.^2/(2*sig^2));
h = h/sum(h); % normalise filter to ensure area under the graph of the data is not altered
v_eye_vel = dze/dt; 
h_eye_vel = dye/dt;
de = sqrt(dze.^2 + dye.^2); % speed of eye movement
de_smooth = conv(de,h,'same')/dt;

if 0
% % % % METHOD 1 % % % %

% apply threshold on eye speed
saccade_thresh = prs.saccade_thresh;
indx_thresh = de_smooth>saccade_thresh;
dindx_thresh = diff(indx_thresh);
t_saccade = find(dindx_thresh>0)*dt;

% remove duplicates by applying a saccade refractory period
min_isi = prs.min_intersaccade;
t_saccade(diff(t_saccade)<min_isi) = [];
t.saccade = t_saccade;

else
% % % % METHOD 2 - based on Larsson et al., 2013 % % % %
per_block = 1;

t.saccade = [];    t.sac_mag = [];    t.sac_vel = [];
if ~NOREC
    plt2 = 0;
    if per_block
        Tstart = t.beg(1);
        I = floor(t.beg(1)/dt);
        O = floor(t.end(end)/dt);
        [t_saccade,sac_mag,sac_vel] = DetectSaccades(ch.zle(I:O),ch.yle(I:O),ch.zre(I:O),ch.yre(I:O),prs,prs.dt,plt2);
        
        t_saccade = t_saccade + Tstart;        
        t.saccade = t_saccade(:);
        t.sac_mag = sac_mag(:);
        t.sac_vel = sac_vel(:);
        
    else
        
        for i = 1:length(t.beg)
            if i > 1; Tstart = t.end(i-1); else; Tstart = t.beg(i); end
            if i < length(t.beg); Tstop = t.beg(i+1); else; Tstop = t.end(i)+1; end
            I = floor(Tstart/dt);
            O = floor(Tstop/dt);
            [t_saccade,sac_mag,sac_vel] = DetectSaccades(ch.zle(I:O),ch.yle(I:O),ch.zre(I:O),ch.yre(I:O),prs,prs.dt,plt2);
            
            t_saccade = t_saccade + Tstart;
            validindx = find(t_saccade >= t.beg(i) & t_saccade <= t.end(i));
            t_saccade = t_saccade(validindx);
            sac_mag = sac_mag(validindx);
            sac_vel = sac_vel(validindx);
            
            t.saccade = [t.saccade ; t_saccade(:)];
            t.sac_mag = [t.sac_mag ; sac_mag(:)];
            t.sac_vel = [t.sac_vel ; sac_vel(:)];
        end
    end
end
end

if 0
figure;plot(t.sac_mag,t.sac_vel,'.');   xlabel('saccade magnitude');    ylabel('Peak Velocity'); xlim([0 60]); %ylim([0 120]);
XX = [ones(length(t.sac_mag),1) t.sac_mag(:)];
y = t.sac_vel(:);
c = regress(y,XX);
res_err = t.sac_vel(:) - XX*c;
hold on;  plot(XX(:,end),XX*c,'r'); title(['slope = ' num2str(c(end),3) ', \sigma_{res} = ' num2str(mean(abs(res_err)))]);
end
%% Use mahalinobis distance for saccade detection
% % apply threshold on eye speed
% thresh = prs.saccade_thresh; % deg/s
% indx_thresh = de_smooth>thresh;
% dindx_thresh = diff(indx_thresh);
% sac_on = find(dindx_thresh>0);
% sac_off = find(dindx_thresh<0); % use it to get the after saccade signal avg
% if sac_on(1) > sac_off(1) % fix event order
%     fixind = find(sac_on(1) < sac_off,1);
%     sac_off = sac_off(fixind:end);
% end
% if length(sac_off) < length(sac_on) % fix length
%     sac_on = sac_on(1:length(sac_off));
% end
% 
% % Mahalinobis distance for saccade detection (compare pre and post saccade eye position)
% nsamp = 275;    skipsamp_pre = 0;     skipsamp_post = 0;    sac_crit = [];
% if plt
% figure;plot(X);grid on;ylim([-100 100])
% vline(t.beg/dt,'k');vline(t.end/dt,'r')
% hold on;plot(de_smooth)
% % xlim([0 6*10^4]);
% ylim([-40 100])
% end
% for n = 2:length(sac_on)-2
%     % pre-saccade
%     LEpos_pre = [ch.zle(sac_on(n)-nsamp:sac_on(n)-skipsamp_pre) ch.yle(sac_on(n)-nsamp:sac_on(n)-skipsamp_pre)];
%     for zz = 1:size(LEpos_pre,2)
%         x = [1:size(LEpos_pre,1)];  X = [ones(length(x),1) x(:)];
%         y = LEpos_pre(:,zz);
%         [c]=regress(y,X);
%         LEpos_pre(:,zz) = LEpos_pre(:,zz) - X*c + LEpos_pre(end,zz);
%     end
%     REpos_pre = [ch.zre(sac_on(n)-nsamp:sac_on(n)-skipsamp_pre) ch.yre(sac_on(n)-nsamp:sac_on(n)-skipsamp_pre)];
%     for zz = 1:size(REpos_pre,2)
%         x = [1:size(REpos_pre,1)];  X = [ones(length(x),1) x(:)];
%         y = REpos_pre(:,zz);
%         [c]=regress(y,X);
%         REpos_pre(:,zz) = REpos_pre(:,zz) - X*c + REpos_pre(end,zz);
%     end
%     BEpos_pre = [0.5*(LEpos_pre(:,1) + REpos_pre(:,1))  0.5*(LEpos_pre(:,2) + REpos_pre(:,2))];
%     % sanity check
%     if 0
%         figure;plot(REpos_pre(:,zz));hold on;hold on;plot(X*c,'r');plot(REpos_pre(:,zz) - X*c + REpos_pre(end,zz),'k')
%     end
%     %
%     % post-saccade
%     LEpos_post = [ch.zle(sac_off(n)+skipsamp_post:sac_off(n)+nsamp) ch.yle(sac_off(n)+skipsamp_post:sac_off(n)+nsamp)];
%     for zz = 1:size(LEpos_post,2)
%         x = [1:size(LEpos_post,1)];  X = [ones(length(x),1) x(:)];
%         y = LEpos_post(:,zz);
%         [c]=regress(y,X);
%         LEpos_post(:,zz) = LEpos_post(:,zz) - X*c + LEpos_post(end,zz);
%     end
%     REpos_post = [ch.zre(sac_off(n)+skipsamp_post:sac_off(n)+nsamp) ch.yre(sac_off(n)+skipsamp_post:sac_off(n)+nsamp)];
%     for zz = 1:size(REpos_post,2)
%         x = [1:size(REpos_post,1)];  X = [ones(length(x),1) x(:)];
%         y = REpos_post(:,zz);
%         [c]=regress(y,X);
%         REpos_post(:,zz) = REpos_post(:,zz) - X*c + REpos_post(end,zz);
%     end
%     BEpos_post = [0.5*(LEpos_post(:,1) + REpos_post(:,1))  0.5*(LEpos_post(:,2) + REpos_post(:,2))];
%     % saccade magnitude
%     LEsac_mag = sqrt((nanmean(LEpos_pre(:,1))-nanmean(LEpos_post(:,1))).^2 + (nanmean(LEpos_pre(:,2))-nanmean(LEpos_post(:,2))).^2);
%     REsac_mag = sqrt((nanmean(REpos_pre(:,1))-nanmean(REpos_post(:,1))).^2 + (nanmean(REpos_pre(:,2))-nanmean(REpos_post(:,2))).^2);
%     BEsac_mag = sqrt((nanmean(BEpos_pre(:,1))-nanmean(BEpos_post(:,1))).^2 + (nanmean(BEpos_pre(:,2))-nanmean(BEpos_post(:,2))).^2);
%     sac_mag(n) = BEsac_mag;
%     % mahalanobis distance
%     LEmahal = sqrt(mahal(nanmean(LEpos_post),LEpos_pre));
%     REmahal = sqrt(mahal(nanmean(REpos_post),REpos_pre));
%     BEmahal = nanmean(sqrt(mahal(BEpos_post,BEpos_pre)));
% %     BEmahal = sqrt(mahal(nanmean(.5*(REpos_post+LEpos_post)),.5*(REpos_pre+LEpos_pre)));
%     if plt; vline(sac_on(n),'b');xlim([sac_on(n)-3*10^4 sac_on(n)+3*10^4]);disp(['mahal = ' num2str(BEmahal) ', mag = ' num2str(BEsac_mag)]);  end
%     
% %     if BEmahal > 20 || (LEmahal > 28 || REmahal > 28)
%     if (BEsac_mag > 3 && (BEmahal > 18 || (LEmahal > 28 || REmahal > 28)))...
%             || BEsac_mag > 6
%         sac_crit = [sac_crit sac_on(n)];
%         if plt; vline(sac_crit(end),'g');   end
% %         figure;plot([ch.zle(sac_on(n)-nsamp:sac_on(n)+nsamp) ch.yle(sac_on(n)-nsamp:sac_on(n)+nsamp)])
% %         figure;plot(LEpos_pre)
% %         figure;plot(LEpos_post)
% %         vline(sac_on(n),'y');vline(sac_off(n),'m')
%     end
% end
% t_saccade = sac_crit(:)*dt;

%% Interpolate nans (blinks)
if ~NOREC
nanx = isnan(ch.zle); t1 = 1:numel(ch.zle); ch.zle(nanx) = interp1(t1(~nanx), ch.zle(~nanx), t1(nanx), 'pchip');
nanx = isnan(ch.zre); t1 = 1:numel(ch.zre); ch.zre(nanx) = interp1(t1(~nanx), ch.zre(~nanx), t1(nanx), 'pchip');
nanx = isnan(ch.yle); t1 = 1:numel(ch.yle); ch.yle(nanx) = interp1(t1(~nanx), ch.yle(~nanx), t1(nanx), 'pchip');
nanx = isnan(ch.yre); t1 = 1:numel(ch.yre); ch.yre(nanx) = interp1(t1(~nanx), ch.yre(~nanx), t1(nanx), 'pchip');

% Re-NaN the blinks to not confuse them with eye movements in analysis
% ch.zle(nanindx,:) = nan;    ch.yle(nanindx,:) = nan;
% ch.zre(nanindx,:) = nan;    ch.yre(nanindx,:) = nan;
end
t.eyeblinks = blinks;
%% detect time points of fixation onsets
v_eye_vel = dze/dt; 
h_eye_vel = dye/dt;
de = sqrt(dze.^2 + dye.^2); % speed of eye movement
de_smooth = conv(de,h,'same')/dt;

fixateduration = prs.fixateduration; fixate_thresh = prctile(de_smooth,90); % set thresh to 90th prctile
fixateduration_samples = round(fixateduration/dt);
fixateindx = false(1,numel(ts) - round(2*fixateduration/dt));
for i=1:(numel(ts) - round(2*fixateduration/dt))
    if mean(de_smooth(i:i+fixateduration_samples)) < fixate_thresh && max(de_smooth(i:i+fixateduration_samples)) < 1.5*fixate_thresh, fixateindx(i) = true; end
end
fixation_switch = diff(fixateindx);
t.fix = ts(fixation_switch>0);

%% replace the broken eye coil (if any) with NaNs
% if var(ch.zle) < 10 || var(ch.zle) > 1000
%     ch.zle(:) = nan;
%     ch.yle(:) = nan;
% end
% if var(ch.zre) < 10 || var(ch.zre) > 1000
%     ch.zre(:) = nan;
%     ch.yre(:) = nan;
% end

%% detect start-of-movement and end-of-movement times for each trial
v_thresh = prs.v_thresh; w_thresh = prs.w_thresh;
v_time2thresh = prs.v_time2thresh;
v = ch.v; w = ch.w;
for j=1:length(t.end)
   % start-of-movement
   if j==1, t.move(j) = t.beg(j); % first trial is special because there is no pre-trial period
   else
       indx = find(v(ts>t.end(j-1) & ts<t.end(j)) > v_thresh,1); % first upward threshold-crossing
       if ~isempty(indx), t.move(j) = t.end(j-1) + indx*dt;
       else, t.move(j) = t.beg(j); end % if monkey never moved, set movement onset = target onset
   end
   % end-of-movement
   indx = find(abs(v(ts>t.move(j) & ts<t.end(j))) < v_thresh & abs(w(ts>t.move(j) & ts<t.end(j))) < w_thresh,1); % first downward threshold-crossing
   if ~isempty(indx), t.stop(j) = t.move(j) + indx*dt;
   else, t.stop(j) = t.end(j); end % if monkey never stopped, set movement end = trial end
   % if monkey stopped prematurely, set movement end = trial end
   if (t.stop(j)<t.beg(j) || (t.stop(j)-t.move(j))<prs.mintrialduration)
       % second attempt to locate t_stop (added 12-04-2019)
       indx = find(abs(v(ts>t.beg(j) & ts<t.end(j))) > v_thresh | abs(w(ts>t.beg(j) & ts<t.end(j))) > w_thresh,1,'last');
       if ~isempty(indx), t.stop(j) = t.beg(j) + indx*dt;
       else, t.stop(j) = t.end(j); end
   end
end
t.move = t.move(:);
t.stop = t.stop(:);
%% extract trials and downsample for storage
dt_original = dt;
dt = dt*prs.factor_downsample;
for j=1:length(t.end)
    % define pretrial period
    pretrial = max(t.beg(j) - t.move(j),0) + prs.pretrial; % extract everything from "movement onset - pretrial" or "target onset - pretrial" - whichever is first
    posttrial = prs.posttrial; % extract everything until "t_end + posttrial"
    for i=1:length(chnames)
        if ~any(strcmp(chnames{i},'mrk'))
            trl(j).continuous.(chnames{i}) = ch.(chnames{i})(ts>t.beg(j)-pretrial & ts<t.end(j)+posttrial);
            if any(strcmp(chnames{i},{'xfp','xmp','yfp','ymp'})) % set position values prior to target onset to nan
                trl(j).continuous.(chnames{i})(1:floor(pretrial/dt_original)) = nan;
            end
            if ~strcmp(chnames{i},'microstim')
                trl(j).continuous.(chnames{i}) = downsample(trl(j).continuous.(chnames{i}),prs.factor_downsample);
            else, trl(j).continuous.ts_microstim = (dt_original:dt_original:length(trl(j).continuous.microstim)*dt_original)' - pretrial;
            end
        end
    end
    trl(j).continuous.ts = (dt:dt:length(trl(j).continuous.(chnames{2}))*dt)' + ...
        ((t.beg(j)-pretrial)<0)*(t.beg(j)-pretrial) + ((t.beg(j)-pretrial)>0)*(-pretrial); % because not enough pretrial before 1st trial  
    trl(j).continuous.firefly = trl(j).continuous.ts>=0 & trl(j).continuous.ts<(0+prs.fly_ONduration);
    trl(j).events.t_beg = t.beg(j);
    trl(j).events.t_end = t.end(j);
    trl(j).events.t_move = t.move(j);
    trl(j).events.t_stop = t.stop(j);
    % saccade time
    trl(j).events.t_sac = t.saccade(t.saccade>(t.beg(j)-pretrial) & t.saccade<t.end(j));
    % fixation start times
    trl(j).events.t_fix = t.fix(t.fix>(t.beg(j)-3*pretrial) & t.fix<(t.end(j)+3*posttrial)); % wider search-range because fixation might be outside trial
    t.fix(t.fix>(t.beg(j)-3*pretrial) & t.fix<(t.end(j)+3*posttrial)) = []; % remove from list
    % reward time
    if any(t.reward>t.beg(j) & t.reward<t.end(j))
        trl(j).logical.reward = true;
        trl(j).events.t_rew = t.reward(t.reward>t.beg(j) & t.reward<t.end(j));
    else
        trl(j).logical.reward = false;
        trl(j).events.t_rew = nan;
    end
end

%% timestamps referenced relative to exp_beg
exp_beg = t.events(find(markers==1,1,'first'));
exp_end = t.events(find(markers==3,1,'last'));

for i=1:length(trl)
    trl(i).events.t_beg = trl(i).events.t_beg - exp_beg;
    trl(i).events.t_rew = trl(i).events.t_rew - exp_beg - trl(i).events.t_beg; % who cares about absolute time?!
    trl(i).events.t_end = trl(i).events.t_end - exp_beg - trl(i).events.t_beg;    
    trl(i).events.t_sac = trl(i).events.t_sac - exp_beg - trl(i).events.t_beg;
    trl(i).events.t_fix = trl(i).events.t_fix - exp_beg - trl(i).events.t_beg;
    trl(i).events.t_move = trl(i).events.t_move - exp_beg - trl(i).events.t_beg;
    trl(i).events.t_stop = trl(i).events.t_stop - exp_beg - trl(i).events.t_beg;
    trl(i).events.t_targ = 0;
    trl(i).events.t_beg_correction = t_beg_correction(i);
    trl(i).events.t_flyON = t_flyON_trl(i) - exp_beg;
    trl(i).events.t_flyON_minus_teleport = t_flyON_trl(i) - t_teleport_trl(i);
end

%% detect trials where targets appeared again during the trial
% for i=1:length(trl)
%     timeindx = trl(i).continuous.ts<trl(i).events.t_end;
%     dPf__dt = [0 ; sqrt(diff(trl(i).continuous.xfp(timeindx)).^2 + diff(trl(i).continuous.yfp(timeindx)).^2)];
%     if findpeaks(dPf__dt,dt*(1:length(dPf__dt)),'MinPeakHeight',prs.minpeakprominence.flypos)>0 % detect peaks
%         trl(i).logical.spurioustarg = true;
%     else
%         trl(i).logical.spurioustarg = false;
%     end
% end

%% downsample continuous data
for i=1:length(chnames)
    if ~any(strcmp(chnames{i},'mrk'))
        try
        temp = ch.(chnames{i})(ts>exp_beg & ts<exp_end);
        temp(isnan(temp)) = 0; 
        ch.(chnames{i}) = temp; 
        ch.(chnames{i}) = decimate(ch.(chnames{i}),prs.factor_downsample);
        catch
        end
    end
end
ts = ts(ts>exp_beg & ts<exp_end) - exp_beg;
ch.ts = decimate(ts,prs.factor_downsample); ch.ts = ch.ts(:);
ch.ntrls = length(trl);
%% extract trials (and downsample for storage)
% % dt = dt*prs.factor_downsample;
% for j=1:length(t.end)
%     % define pretrial period
% %     pretrial = prs.pretrial; % extract everything from "movement onset - pretrial" or "target onset - pretrial" - whichever is first
% %     posttrial = prs.posttrial; % extract everything until "t_end + posttrial"
% 
% % fix start of trial
% if sum(ts>t.beg(j) & ts<t.end(j)) <= 833; maxlength = sum(ts>t.beg(j) & ts<t.end(j)); else; maxlength = 833; end
% reset_detect = [0 ; diff(ch.ymp(ts>t.beg(j) & ts<t.end(j)))/dt];
% reset_detect = find(reset_detect(1:maxlength) >= 8000,1,'last');
% 
% if ~isempty(reset_detect)
%     t.beg(j) = t.beg(j) + ts(reset_detect);
% end
% 
%     for i=1:length(chnames)
%         if ~any(strcmp(chnames{i},'mrk'))
%             trl(j).continuous.(chnames{i}) = ch.(chnames{i})(ts>t.beg(j) & ts<t.end(j));
% %             trl(j).continuous.(chnames{i}) = downsample(trl(j).continuous.(chnames{i}),prs.factor_downsample);
%         end
%     end
% %     trl(j).continuous.ts = (dt:dt:length(trl(j).continuous.(chnames{2}))*dt)';
% %     trl(j).continuous.firefly = trl(j).continuous.ts>=0 & trl(j).continuous.ts<(0+prs.fly_ONduration);
%     trl(j).continuous.ts = (dt:dt:t.end(j)-t.beg(j))';
%     trl(j).continuous.firefly = logical(trl(j).continuous.FFdraw);
%     trl(j).events.t_beg = t.beg(j);
% %     indx = find(t.shutter >= t.beg(j) & t.shutter <= t.end(j));
% %     trl(j).events.t_shutter = t.shutter(indx);
%     trl(j).events.t_end = t.end(j);
% %     trl(j).events.t_move = t.move(j);
% %     trl(j).events.t_stop = t.stop(j);
%     % saccade time
%     sacindx = (t.saccade>t.beg(j)) & (t.saccade<t.end(j));
%     trl(j).events.t_sac = t.saccade(sacindx);
%     trl(j).events.sac_mag = t.sac_mag(sacindx);
%     trl(j).events.sac_vel = t.sac_vel(sacindx);
% %     trl(j).events.fix = t.fix(t.fix>(t.beg(j)) & t.fix<t.end(j));
% 
%     trl(j).continuous.blink = t.eyeblinks(ts > (t.beg(j)) & ts < t.end(j));
%     if NOREC
%         trl(j).events.NOREC =  blankindx(ts > (t.beg(j)) & ts < t.end(j));
%         if ~isempty(trl(j).events.NOREC) % trials that are fine before recording got screwed
%             trl(j).events.NOREC = 1;
%         else
%             trl(j).events.NOREC = 0;
%         end
%     else
%         trl(j).events.NOREC = 0;
%     end
%     % reward time
% %     if any(t.reward>t.beg(j) & t.reward<t.end(j))
% %         trl(j).logical.reward = true;
% %         trl(j).events.t_rew = t.reward(t.reward>t.beg(j) & t.reward<t.end(j));
% %     else
% %         trl(j).logical.reward = false;
% %         trl(j).events.t_rew = nan;
% %     end
%      % ptb time
% %     if any(t.ptb>t.beg(j) & t.ptb<t.end(j))
% %         trl(j).logical.ptb = true;
% %         trl(j).events.t_ptb = t.ptb(t.ptb>t.beg(j) & t.ptb<t.end(j));
% %     else
% %         trl(j).logical.ptb = false;
% %         trl(j).events.t_ptb = nan;
% %     end
% end
%% remove trials that include moog crash (detect as jump in position)
% if exist('spike2terminated')
%     moog_crash = [];
%     for i = 1:length(trl) 
%             Y = trl(i).continuous.ymp;
%             X = trl(i).continuous.xmp;
%             scany = diff(Y(1:15:end));
%             scanx = diff(X(1:15:end));
%             [maxy,indy] = max(scany);
%             [maxx,indx] = max(scanx);
%             if [maxx maxy] > 60
%                 moog_crash = [moog_crash i];
%             end
%     end
%     if ~isempty(moog_crash)
% %         trl(moog_crash) = [];
%         disp (['Trial ' num2str(moog_crash) ' MUST be removed because of possible moog crash !'])
%     else
%     disp('Moog Crash wasn''t found, make sure about reason of termination of Spike2')    
%     end
% else
%     moog_crash = [];
% end
%% define gaussian filter (done)
% % use Median instead of Gaussian filter? No
% sig = prs.filtwidth; %filter width
% sz = prs.filtsize; %filter size
% t2 = linspace(-sz/2, sz/2, sz);
% h = exp(-t2.^2/(2*sig^2));
% h = h/sum(h); % normalise filter to ensure area under the graph of the data is not altered

%% filter position, velocity and sensor channels (done)
% % define Butterworth filter
% SR = 1/dt;
% [b,a] = butter(2,30/(SR/2),'low'); % for the sensors
% screwed = []; % keep track of screwed trials
% for i=1:length(trl)
%     %     if i == 13 % use this to fix temporarily wrong trials, if any
%     %         keyboard;
%     %     end
%     for j = 1:length(chnames)
%         if ~any(strcmp(chnames{j},{'FFdraw','ts','mrk','firefly','yle','yre','zle','zre'}))
%             if isempty(trl(i).continuous.(chnames{j}))
%                 screwed = [screwed i];
%                 continue;
%             end
%             if any(strcmp(chnames{j},{'xac','yac','vrol','vyaw','vpit'}))
%                 trl(i).continuous.(chnames{j}) = filtfilt(b,a,trl(i).continuous.(chnames{j})); % butterowrth filter for the sensors
%             else
%                 %           trl(i).continuous.(chnames{j}) = medfilt1(trl(i).continuous.(chnames{j}),5); % median filter
%                 trl(i).continuous.(chnames{j}) = [trl(i).continuous.(chnames{j})(1)*ones(1,100)';...
%                     trl(i).continuous.(chnames{j});...
%                     trl(i).continuous.(chnames{j})(end)*ones(1,100)']; % adding first and last value in edges
%                 trl(i).continuous.(chnames{j}) = conv(trl(i).continuous.(chnames{j})(1:end),h,'same'); % gaussian filter for data
%                 %           ch.(chnames{i}) = ch.(chnames{i})(sz/2+1:end);
%                 trl(i).continuous.(chnames{j}) = trl(i).continuous.(chnames{j})(1+100:end-100);
%             end
%         end
%     end
% end
% if screwed
%     screwed = unique(screwed);
%     disp(['There are screwed trials No. ' num2str(screwed) '!!!!! Check AddSMRData'])
%     keyboard;
% end
%% Detect saccades (per trial)
%% detect saccade times - based on Larsson et al., 2013
if 0
for i = 1:length(trl)
    if ~NOREC
        plt2 = 0;
        
        ZLE = trl(i).continuous.zle;
        YLE = trl(i).continuous.yle;
        ZRE = trl(i).continuous.zre;
        YRE = trl(i).continuous.yre;
                
        [t_saccade,sac_mag,sac_vel] = DetectSaccades(ZLE,YLE,ZRE,YRE,prs,prs.dt,plt2);
        trl(i).events.t_sac = t_saccade(:);
        trl(i).events.sac_mag = sac_mag(:);
        trl(i).events.sac_vel = sac_vel(:);
    else
        trl(i).events.t_sac = [];
        trl(i).events.sac_mag = [];
        trl(i).events.sac_vel = [];
    end
end
end
%% set position values prior to target onset to nan
% for j=1:length(trl)
%     for i=1:length(chnames)
%         if any(strcmp(chnames{i},{'xmp','ymp'}))
%             trl(j).continuous.(chnames{i})(trl(j).continuous.ts<0) = nan; % target onset happens exactly at t_beg ?????
%         end
%     end
% end
%% debugging plots

% figure;plot(ch.ts,ch.v)
% hold on;plot(ch.ts,ch.FFdraw*100)
% vline(t.beg)
% vline(t.end,'k')

t.ntrls = length(trl);