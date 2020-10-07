function [trl,moog_crash,ch,t] = AddSMRData(data,cal,gains,prs)
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
chno.mrk = find(strcmp(ch_title,'Keyboard')); % keyboard as the marker
chno.tsi = find(strcmp(ch_title,'TrialSta')); % trial start indicator 

chno.yle = find(strcmp(ch_title,'LDy')); 
chno.zle = find(strcmp(ch_title,'LDz'));
chno.yre = find(strcmp(ch_title,'RDy')); 
chno.zre = find(strcmp(ch_title,'RDz'));

chno.phi = find(strcmp(ch_title,'MonkeyYa')); % monkey yaw position
chno.xmp = find(strcmp(ch_title,'MonkeyX')); 
chno.ymp = find(strcmp(ch_title,'MonkeyY'));
chno.v = find(strcmp(ch_title,'ForwardV')); 
chno.w = find(strcmp(ch_title,'AngularV'));
chno.mtr = find(strcmp(ch_title,'MotorCom')); % yaw motor command
chno.xac = find(strcmp(ch_title,'AccX')); % acceleration x-axis
chno.yac = find(strcmp(ch_title,'AccY')); % acceleration y-axis

chno.vrol = find(strcmp(ch_title,'VelRoll')); % roll velocity 
chno.vyaw = find(strcmp(ch_title,'VelYaw')); % yaw velocity 
chno.vpit = find(strcmp(ch_title,'VelPitch')); % pitch velocity 


%% scale (done)
scaling.t = data(chno.mrk).hdr.tim.Scale*data(chno.mrk).hdr.tim.Units; % for markers

scaling.tsi = data(chno.tsi).hdr.adc.Scale; offset.tsi = data(chno.tsi).hdr.adc.DC; % for pulse

scaling.yle = data(chno.yle).hdr.adc.Scale; offset.yle = data(chno.yle).hdr.adc.DC;
scaling.yre = data(chno.yre).hdr.adc.Scale; offset.yre = data(chno.yre).hdr.adc.DC; 
scaling.zle = data(chno.zle).hdr.adc.Scale; offset.zle = data(chno.zle).hdr.adc.DC; 
scaling.zre = data(chno.zre).hdr.adc.Scale; offset.zre = data(chno.zre).hdr.adc.DC;

scaling.phi = data(chno.phi).hdr.adc.Scale; offset.phi = data(chno.phi).hdr.adc.DC;
scaling.xmp = data(chno.xmp).hdr.adc.Scale; offset.xmp = data(chno.xmp).hdr.adc.DC;
scaling.ymp = -data(chno.ymp).hdr.adc.Scale; offset.ymp = -data(chno.ymp).hdr.adc.DC; % !!!!!SIGN
scaling.v = data(chno.v).hdr.adc.Scale; offset.v = data(chno.v).hdr.adc.DC;
scaling.w = data(chno.w).hdr.adc.Scale; offset.w = data(chno.w).hdr.adc.DC;
scaling.mtr = data(chno.mtr).hdr.adc.Scale; offset.mtr = data(chno.mtr).hdr.adc.DC;
scaling.xac = data(chno.xac).hdr.adc.Scale; offset.xac = data(chno.xac).hdr.adc.DC;
scaling.yac = data(chno.yac).hdr.adc.Scale; offset.yac = data(chno.yac).hdr.adc.DC;

scaling.vrol = data(chno.vrol).hdr.adc.Scale; offset.vrol = data(chno.vrol).hdr.adc.DC;
scaling.vyaw = data(chno.vyaw).hdr.adc.Scale; offset.vyaw = data(chno.vyaw).hdr.adc.DC;
scaling.vpit = data(chno.vpit).hdr.adc.Scale; offset.vpit = data(chno.vpit).hdr.adc.DC;

%% load relevant channels (done)
chnames = fieldnames(chno); MAX_LENGTH = inf; dt = [];
for i=1:length(chnames)
    if ~any(strcmp(chnames{i},'mrk'))
        ch.(chnames{i}) = double(data(chno.(chnames{i})).imp.adc)*scaling.(chnames{i}) + offset.(chnames{i});
        dt = [dt prod(data(chno.(chnames{i})).hdr.adc.SampleInterval)];
        MAX_LENGTH = min(length(ch.(chnames{i})),MAX_LENGTH);
    end
end

ch.('mrk') = double(data(chno.('mrk')).imp.tim)*scaling.t;

if length(unique(dt))==1
    dt = dt(1);
else
   error('channels must all have identical sampling rates');
end
%% fix length of data
for i=1:length(chnames)
    if ~any(strcmp(chnames{i},'mrk'))
        ch.(chnames{i}) = ch.(chnames{i})(1:MAX_LENGTH); % butterowrth filter for the sensors
    end
end
%% event markers (done)
markers = data(chno.mrk).imp.mrk(:,1); % from 'keyboard'
%% event times (done)
ts = dt:dt:dt*MAX_LENGTH;
ch.ts = ts'; % in seconds
t.events = double(data(chno.mrk).imp.tim)*scaling.t;
% 115 = s (start), 117 = u (shutter close) , 118 = v (vestibular only) , 114 = r (end of block)
t.beg = t.events(markers ==115 | markers==118); 
t.end = [t.beg ; t.events(markers == 114)]; t.end(1) = []; % first marker is only start
t.beg = t.beg(1:length(t.end));
t.shutter = t.events(markers == 117);
t.feedback = t.events(markers == 104);

% % previous marking using TSI
tsi = ch.tsi;
tsi(tsi>=.2)= 5;
tsi(tsi< .2)= 0;
pulse_max = max(tsi);

% make sure that last element of tsi is 0
% tsi(end) = 0;
   
% % sanity check
% figure;plot(ts,tsi);vline(t.beg);hold on;plot(ts,.01*ch.xmp);plot(ts,.01*ch.ymp);
% check that the TSI starts at 0 and not 5
count = 1;
while tsi(count) == 5
    tsi(count) = 0;
    count = count+1;
end
if count > 1
    disp('TSI starts at 5 for this block.')
end
% use markers to choose events based on tsi (change markers)
pulseOFF = find(diff(tsi) == -pulse_max);
pulseON = find(diff(tsi) == pulse_max);

for i = 1:length(t.beg)       
%         if i ~= length(t.beg)               % keep the last marker intact, comes from another marker(114)
            indx_end = find(pulseON >= t.beg(i)/dt, 1);
            if isempty(indx_end)
                t.beg(end) = [];
                t.end(end) = [];
                
                % keep track 
                trial_terminated = i-1;
                spike2terminated = ['Trials cut at ' num2str(trial_terminated) ', probably terminated Spike2 !'];
                disp(spike2terminated)
            else
                t.end(i) = pulseON(indx_end)*dt;
            end
%         end
end

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
% % figure;plot(ts,tsi);vline(t.beg);vline(t.end,'k');hold on;plot(ts,.01*ch.xmp);plot(ts,.01*ch.ymp);
%%
t.events = sort([t.beg ; t.end ; t.shutter ; t.feedback]);
ch.tsi = tsi;

% % now find events
% if find(markers == 104) % if FEEDBACK
%
%     indx_beg = find(diff(tsi) == -pulse_max); % choose the end of the pulse, since it comes after the reset of position
%     indx_end = find(diff(tsi) == pulse_max); % choose the start of the pulse, since it is the command to reset position
%     t.beg = ts(indx_beg(1:2:end));
%     t.end = ts(indx_end(1:2:end));
%     t.end(1) = [];
%     t.feedback = ts(indx_beg(2:2:end));
%     t.beg = t.beg(1:length(t.end));
%     t.events = sort([ts(indx_beg) ts(indx_end)]);
% else
%     indx_beg = find(diff(tsi) == -pulse_max); % choose the end of the pulse, since it comes after the reset of position
%     indx_end = find(diff(tsi) == pulse_max); % choose the start of the pulse, since it is the command to reset position
%
%     t.events = sort([ts(indx_beg) ts(indx_end)]);
%     t.events(1) = [];     t.events(end) = []; % see tsi plot to check events, 1st event start, last event end
%
%     t.beg = t.events(1:2:end);
%     t.end = t.events(2:2:end);
%     t.beg = t.beg(1:length(t.end));
% end
% % indxon = find(round(diff(ch.tsi)) == pulse_max);
% % indxoff = find(round(diff(ch.tsi)) == -pulse_max);
% % pulse_duration  = indxoff - indxon;
%% Eye Calibration
%% Re-calibrate using EyeCal file, output offset and scale
% extract scale and offset from data (watch the order: scale 1st, offset 2nd)
g{1} = [gains.yle_scale gains.yle_offset];  g{2} = [gains.zle_scale gains.zle_offset];
g{3} = [gains.yre_scale gains.yre_offset];  g{4} = [gains.zre_scale gains.zre_offset];
IOD = gains.eye_dist;
max_MSE_cal_fit = prs.max_MSE_cal_fit;
% check channel headers
if ~isempty(cal)
    nch_cal = length(cal);
    ch_title_cal = cell(1,nch_cal);
    hdr_cal = {cal.hdr};
    for i=1:nch_cal
        if ~isempty(hdr_cal{i})
            ch_title_cal{i} = hdr_cal{i}.title;
        else
            ch_title_cal{i} = 'nan';
        end
    end
    % channel titles
    chno_cal.yle = find(strcmp(ch_title_cal,'LDy'));
    chno_cal.zle = find(strcmp(ch_title_cal,'LDz'));
    chno_cal.yre = find(strcmp(ch_title_cal,'RDy'));
    chno_cal.zre = find(strcmp(ch_title_cal,'RDz'));
    chno_cal.tar = find(strcmp(ch_title_cal,'TextMark')); % nothing
    chno_cal.key = find(strcmp(ch_title_cal,'Keyboard')); % nothing
    chno_cal.mrk = find(strcmp(ch_title_cal,'marker')); % target positions
    
    % scale and offset
    scaling_cal.t = cal(chno_cal.mrk).hdr.tim.Scale*cal(chno_cal.mrk).hdr.tim.Units; % for markers
        
    scaling_cal.yle = cal(chno_cal.yle).hdr.adc.Scale; offset_cal.yle = cal(chno_cal.yle).hdr.adc.DC;
    scaling_cal.yre = cal(chno_cal.yre).hdr.adc.Scale; offset_cal.yre = cal(chno_cal.yre).hdr.adc.DC;
    scaling_cal.zle = cal(chno_cal.zle).hdr.adc.Scale; offset_cal.zle = cal(chno_cal.zle).hdr.adc.DC;
    scaling_cal.zre = cal(chno_cal.zre).hdr.adc.Scale; offset_cal.zre = cal(chno_cal.zre).hdr.adc.DC;
    
    % load channels
    chnames_cal = fieldnames(chno_cal); dt = [];
    for i=1:length(chnames_cal)
        if ~any(strcmp(chnames_cal{i},{'key','mrk','tar'}))
            ch_cal.(chnames_cal{i}) = double(cal(chno_cal.(chnames_cal{i})).imp.adc)*scaling_cal.(chnames_cal{i}) + offset_cal.(chnames_cal{i});
            dt = [dt prod(cal(chno_cal.(chnames_cal{i})).hdr.adc.SampleInterval)];
            MAX_LENGTH = min(length(ch_cal.(chnames_cal{i})),MAX_LENGTH);
        end
    end
    
    ch_cal.('tim') = double(cal(chno_cal.('mrk')).imp.tim(:,1))*scaling_cal.t;
    ch_cal.('mrk') = double(cal(chno_cal.('mrk')).imp.mrk(:,1));
    
    % fix length of data
    for i=1:length(chnames_cal)
        if ~any(strcmp(chnames_cal{i},{'key','mrk','tar'}))
            ch_cal.(chnames_cal{i}) = ch_cal.(chnames_cal{i})(1:MAX_LENGTH); 
        end
    end
    
    if length(unique(dt(~isnan(dt))))==1
        dt = dt(1);
    else
        error('channels must all have identical sampling rates');
    end
    ts_cal = dt:dt:dt*length(ch_cal.yle);
    ch_cal.ts = ts_cal';
    
    % save target positions
    [ZLEcal,YLEcal,ZREcal,YREcal,indR] = EyeCalPos(ch_cal,gains);
    ch_cal.zle_cal = ZLEcal;    ch_cal.yle_cal = YLEcal;
    ch_cal.zre_cal = ZREcal;    ch_cal.yre_cal = YREcal;
        
    % Remove blinks and smooth
    X = [ch_cal.zle ch_cal.zre ch_cal.yle ch_cal.yre];
    X = ReplaceWithNans(X, 3, 500,'normal');
    ch_cal.zle = X(:,1); ch_cal.zre = X(:,2); ch_cal.yle = X(:,3); ch_cal.yre = X(:,4);
    
    % interpolate
    nanx = isnan(ch_cal.yle); t1 = 1:numel(ch_cal.yle); ch_cal.yle(nanx) = interp1(t1(~nanx), ch_cal.yle(~nanx), t1(nanx), 'pchip');
    nanx = isnan(ch_cal.zle); t1 = 1:numel(ch_cal.zle); ch_cal.zle(nanx) = interp1(t1(~nanx), ch_cal.zle(~nanx), t1(nanx), 'pchip');
    nanx = isnan(ch_cal.yre); t1 = 1:numel(ch_cal.yre); ch_cal.yre(nanx) = interp1(t1(~nanx), ch_cal.yre(~nanx), t1(nanx), 'pchip');
    nanx = isnan(ch_cal.zre); t1 = 1:numel(ch_cal.zre); ch_cal.zre(nanx) = interp1(t1(~nanx), ch_cal.zre(~nanx), t1(nanx), 'pchip');

    sig = 12*prs.filtwidth; %filter width
    sz = 30*prs.filtsize; %filter size
    t2 = linspace(-sz/2, sz/2, sz);
    h = exp(-t2.^2/(2*sig^2));
    h = h/sum(h); % normalise filter to ensure area under the graph of the data is not altered
    
    ch_cal.zle = conv(ch_cal.zle,h,'same'); ch_cal.zre = conv(ch_cal.zre,h,'same');
    ch_cal.yle = conv(ch_cal.yle,h,'same'); ch_cal.yre = conv(ch_cal.yre,h,'same');
        
    % fit Gain and Offset
    indx = [];
    nanind = find(~isnan(YLEcal),1);
    nt = length(YLEcal);
    
    indx{1} = 1:length(YLEcal); EYE{1} = ch_cal.yle(indx{1});   CAL{1} = YLEcal(indx{1}); CAL{1}(CAL{1}== 0) = nan;
    indx{2} = 1:length(ZLEcal); EYE{2} = ch_cal.zle(indx{2});   CAL{2} = ZLEcal(indx{2}); CAL{2}(CAL{2}== 0) = nan;
    indx{3} = 1:length(YREcal); EYE{3} = ch_cal.yre(indx{3});   CAL{3} = YREcal(indx{3}); CAL{3}(CAL{3}== 0) = nan;
    indx{4} = 1:length(ZREcal); EYE{4} = ch_cal.zre(indx{4});   CAL{4} = ZREcal(indx{4}); CAL{3}(CAL{4}== 0) = nan;
    
    EYE{1}(end-1000:end) = nan;  EYE{2}(end-1000:end) = nan;  EYE{3}(end-1000:end) = nan;  EYE{4}(end-1000:end) = nan;
    
    % set rotation matrix
    
%     R = @(phi)[cosd(phi) -sind(phi)...
%         ; sind(phi) cosd(phi)];
%     
%     A = @(a)[a(1)  0  ;
%         0   a(2)];
%     
%     B = @(b)[b(1) b(2)];
%     
%     E{1} = [EYE{2}(:) EYE{1}(:)];
%     E{2} = [EYE{4}(:) EYE{3}(:)];
%     C{1} = [CAL{2}(:) CAL{1}(:)];
%     C{2} = [CAL{4}(:) CAL{3}(:)];
%     
%     fmin = @(p) (E{n}*A(p(2:3)) + B(p(4:5)));
%     for n = 1:length(E)
%         [params{n},fval{n}] = fminsearch(@(p) sum(sum((C{n} - fmin(p(:))).^2,2)),[0 0 0 0 0]);
%     end
%     
%     % sanity check
%     figure;
%     n = 1;
%     tempL = fmin(params{1});
%     yL = tempL(:,2);
%     zL = tempL(:,1);
%     subplot(1,2,1);plot(CAL{1},CAL{2},'o');hold on;
%     plot(yL,zL);title('Left Eye');
%     n = 2;
%     tempR = fmin(params{2});
%     yR = tempR(:,2);
%     zR = tempR(:,1);
%     subplot(1,2,2);plot(CAL{3},CAL{4},'o');hold on;
%     plot(yR,zR);title('Right Eye')

    % fit offset and scale
    for n = 1:length(indx)
        [c{n},fval{n}] = fminsearch(@(p) nansum((CAL{n}(:) - p(1)*EYE{n}(:) - p(2)).^2),[0 0]);%         [a{n}] = regress(CAL{n},[ones(length(CAL{n}),1) EYE{n}]);
    end
    
    plt=0;

    % use scale and offset from block data
    for n = 1:length(indx)
        if plt
        figure;
        % fit
        subplot(1,2,1);plot(CAL{n});hold on;plot(EYE{n}*c{n}(1) + c{n}(2));
        legend({'target position','eye position'});axis([nanind nt -30 30]);
        title([chnames_cal{n} ' - G_f_i_t = ' num2str(c{n}(1)) ', C_f_i_t = ' num2str(c{n}(2))])
        % data
        subplot(1,2,2);plot(CAL{n});hold on;plot(EYE{n}*g{n}(1) + g{n}(2));
        legend({'target position','eye position'});axis([nanind nt -30 30]);
        title([chnames_cal{n} ' - G_d_a_t_a = ' num2str(g{n}(1)) ', C_d_a_t_a = ' num2str(g{n}(2))])
        end
    end
    % plot target cross
    plt = 0;
    if plt
        YLEcal(isnan(YLEcal)) = 0;        ZLEcal(isnan(ZLEcal)) = 0;
        YREcal(isnan(YREcal)) = 0;        ZREcal(isnan(ZREcal)) = 0;
    figure;
    % fit
    subplot(2,2,1);plot(YLEcal,ZLEcal,'o');hold on;
    plot(EYE{1}*c{1}(1) + c{1}(2),EYE{2}*c{2}(1) + c{2}(2));title('Left Eye - fit')
    subplot(2,2,2);plot(YLEcal,ZLEcal,'o');hold on;
    plot(EYE{3}*c{3}(1) + c{3}(2),EYE{4}*c{4}(1) + c{4}(2));title('Right Eye - fit')
    % data
    subplot(2,2,3);plot(YREcal,ZREcal,'o');hold on;
    plot(EYE{1}*g{1}(1) + g{1}(2),EYE{2}*g{2}(1) + g{2}(2));title('Left Eye - data')
    subplot(2,2,4);plot(YREcal,ZREcal,'o');hold on;
    plot(EYE{3}*g{3}(1) + g{3}(2),EYE{4}*g{4}(1) + g{4}(2));title('Right Eye - data')
    end
    % choose calibration (manual or fitted)
    err_fit = [];   err_man = [];
    for n = 1:length(indx)
        % fit
        err_fit(n) = nanmean((CAL{n}(:) -(EYE{n}(:)*c{n}(1) + c{n}(2))).^2);
        % manual
        err_man(n) = nanmean((CAL{n}(:) -(EYE{n}(:)*g{n}(1) + g{n}(2))).^2);
    end
    
    for n = 1:length(err_fit)
        
        if (err_fit(n) == 0) || sum(isnan(EYE{n}))/length(EYE{n}) > .80 % problematic calibration (Baptiste)
            G(n) = g(n);
        else % normal calibrations
            if err_fit(n) < err_man(n)
                G(n) = c(n);
            else
                G(n) = g(n);
            end
        end
    end
    if any(err_fit > max_MSE_cal_fit)
        disp('EYE RECORDING THROWN AWAY BECAUSE OF BAD CALIBRATION, err = '); disp(num2str(err_fit(:)));
        suptitle('THROWN AWAY')
        NOREC = 1;
    else
        NOREC = 0;
    end
else
    G = [];
    NOREC = 1;
end
%% Eye Data
% determine whether there is eye recording
plt = 0;
X = [ch.zle ch.zre ch.yle ch.yre];  
Xstd = std(X(1:round(size(X,1)/3),:),1);
if any(Xstd < prs.rec_thresh) || NOREC == 1
    NOREC = 1;
    ch.zle(1:end) = 0;    ch.yle(1:end) = 0;
    ch.zre(1:end) = 0;    ch.yre(1:end) = 0;
    nanindx = zeros(length(ch.zle),1);
    blankindx = 1:length(ch.zle);
else
    NOREC = 0;
    %% Remove blinks and smooth
    % subtract offset to align blinks and detect shifts clearly, then replace
    Xoffset = [g{2}(2) g{4}(2) g{1}(2) g{3}(2)];
    Xscale = [g{2}(1) g{4}(1) g{1}(1) g{3}(1)];
    X = X - Xoffset;
    if plt
    figure;plot(X);grid on;ylim([-100 100])
    end
    [~,nanindx] = ReplaceWithNans(X, prs.blink_thresh_pos,prs.nanpadding,'normal');
%     [~,nanindx2] = ReplaceWithNans(X, prs.blink_thresh_vel,prs.nanpadding,'derivative');
%     nanindx = nanindx1 || nanindx2;
    % patch nans of 150 ms around blinks
    nanpatch = floor(.15/dt);
    indx_right = circshift(nanindx,nanpatch);
    indx_left = circshift(nanindx,-nanpatch);
    nanindx = nanindx|indx_right|indx_left;
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
    if isempty(G)
        X = X + Xoffset;
        newoffset = Xoffset;
    else
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
    if plt
    figure;plot(Xavg);grid on;hline(Xavgoffset),vline(2*filtwidth);
    end
    % find destroyed recording
    if sum(any(abs(Xavg) > prs.blink_thresh + abs(newoffset))) && any(any(Xmovstd < prs.rec_thresh,1))
       % 3 criteria (blink threshold, recording threshold, noise)
        blinkindx = sum(abs(Xavg) > prs.blink_thresh + abs(newoffset),2) > 0;
        recindx = sum(Xmovstd < prs.rec_thresh,2) > 0;
        noiseindx = sum(Xmovstd > 3*nanmean(Xmovstd,1),2) > 0;
  
        blankindx = blinkindx | recindx | noiseindx;
    else
        blankindx = [];
    end
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
    if plt
        figure;plot([ch.zle ch.yle ch.zre ch.yre]);grid on;ylim([-100 100]);
    end
    % remove part where recording is fucked
    if ~isempty(blankindx)
        ch.zle(blankindx,:) = 0;   ch.yle(blankindx,:) = 0;
        ch.zre(blankindx,:) = 0;   ch.yre(blankindx,:) = 0;
    end
end
blinks = nanindx;
%% Detect saccades
%% detect saccade times - based on Larsson et al., 2013
t.saccade = [];    t.sac_mag = [];    t.sac_vel = [];
if ~NOREC
    plt2 = 0;
   for i = 1:length(t.beg)
       if i > 1; Tstart = t.end(i-1); else; Tstart = t.beg(i)-1; end
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
% if ~NOREC
% fixateduration = prs.fixateduration; fixate_thresh = prctile(de_smooth,90); % set thresh to 90th prctile
% fixateduration_samples = round(fixateduration/dt);
% fixateindx = false(1,numel(ts) - round(2*fixateduration/dt));
% for i=1:(numel(ts) - round(2*fixateduration/dt))
%     if mean(de_smooth(i:i+fixateduration_samples)) < fixate_thresh && max(de_smooth(i:i+fixateduration_samples)) < 1.5*fixate_thresh
%         fixateindx(i) = true; 
%     end
% end
% fixation_switch = diff(fixateindx);
% t.fix = ts(fixation_switch>0);
% else
%     t.fix = [];
% end
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
% v_thresh = prs.v_thresh;
% v_time2thresh = prs.v_time2thresh;
% v = ch.v;
% for j=1:length(t.end)
%    % start-of-movement
%    if j==1, t.move(j) = t.beg(j); % first trial is special because there is no pre-trial period
%    else
%        indx = find(v(ts>t.end(j-1) & ts<t.end(j)) > v_thresh,1); % first upward threshold-crossing
%        if ~isempty(indx), t.move(j) = t.end(j-1) + indx*dt;
%        else, t.move(j) = t.beg(j); end % if monkey never moved, set movement onset = target onset
%    end
%    % end-of-movement
%    indx = find(v(ts>t.move(j) & ts<t.end(j)) < v_thresh,1); % first downward threshold-crossing
%    if ~isempty(indx), t.stop(j) = t.move(j) + indx*dt;
%    else, t.stop(j) = t.end(j); end % if monkey never stopped, set movement end = trial end
%    % if monkey stopped prematurely, set movement end = trial end
%    if (t.stop(j)<t.beg(j) || t.stop(j)-t.move(j)<0.5), t.stop(j) = t.end(j); end
% end
% 
%% extract trials (and downsample for storage)
% dt = dt*prs.factor_downsample;
for j=1:length(t.end)
    % define pretrial period
%     pretrial = prs.pretrial; % extract everything from "movement onset - pretrial" or "target onset - pretrial" - whichever is first
%     posttrial = prs.posttrial; % extract everything until "t_end + posttrial"
% fix start of trial
reset_detect = [0 ; diff(ch.ymp(ts>t.beg(j) & ts<t.end(j)))/dt];
reset_detect = find(reset_detect(1:833) >= 8000,1,'last');

if ~isempty(reset_detect)
    t.beg(j) = t.beg(j) + ts(reset_detect);
end

    for i=1:length(chnames)
        if ~any(strcmp(chnames{i},'mrk'))
            trl(j).continuous.(chnames{i}) = ch.(chnames{i})(ts>t.beg(j) & ts<t.end(j));
%             trl(j).continuous.(chnames{i}) = downsample(trl(j).continuous.(chnames{i}),prs.factor_downsample);
        end
    end
    trl(j).continuous.ts = (dt:dt:length(trl(j).continuous.(chnames{2}))*dt)';
    trl(j).continuous.firefly = trl(j).continuous.ts>=0 & trl(j).continuous.ts<(0+prs.fly_ONduration);
    trl(j).events.t_beg = t.beg(j);
    indx = find(t.shutter >= t.beg(j) & t.shutter <= t.end(j));
    trl(j).events.t_shutter = t.shutter(indx);
    trl(j).events.t_end = t.end(j);
%     trl(j).events.t_move = t.move(j);
%     trl(j).events.t_stop = t.stop(j);
    % saccade time
    sacindx = (t.saccade>t.beg(j)) & (t.saccade<t.end(j));
    trl(j).events.t_sac = t.saccade(sacindx);
    trl(j).events.sac_mag = t.sac_mag(sacindx);
    trl(j).events.sac_vel = t.sac_vel(sacindx);
%     trl(j).events.fix = t.fix(t.fix>(t.beg(j)) & t.fix<t.end(j));

    trl(j).continuous.blink = t.eyeblinks(ts > (t.beg(j)) & ts < t.end(j));
    if NOREC
        trl(j).events.NOREC =  blankindx(ts > (t.beg(j)) & ts < t.end(j));
        if ~isempty(trl(j).events.NOREC) % trials that are fine before recording got screwed
            trl(j).events.NOREC = 1;
        else
            trl(j).events.NOREC = 0;
        end
    else
        trl(j).events.NOREC = 0;
    end
    % reward time
%     if any(t.reward>t.beg(j) & t.reward<t.end(j))
%         trl(j).logical.reward = true;
%         trl(j).events.t_rew = t.reward(t.reward>t.beg(j) & t.reward<t.end(j));
%     else
%         trl(j).logical.reward = false;
%         trl(j).events.t_rew = nan;
%     end
     % ptb time
%     if any(t.ptb>t.beg(j) & t.ptb<t.end(j))
%         trl(j).logical.ptb = true;
%         trl(j).events.t_ptb = t.ptb(t.ptb>t.beg(j) & t.ptb<t.end(j));
%     else
%         trl(j).logical.ptb = false;
%         trl(j).events.t_ptb = nan;
%     end
end
%% remove trials that include moog crash (detect as jump in position)
if exist('spike2terminated')
    moog_crash = [];
    for i = 1:length(trl) 
            Y = trl(i).continuous.ymp;
            X = trl(i).continuous.xmp;
            scany = diff(Y(1:15:end));
            scanx = diff(X(1:15:end));
            [maxy,indy] = max(scany);
            [maxx,indx] = max(scanx);
            if [maxx maxy] > 60
                moog_crash = [moog_crash i];
            end
    end
    if ~isempty(moog_crash)
%         trl(moog_crash) = [];
        disp (['Trial ' num2str(moog_crash) ' MUST be removed because of possible moog crash !'])
    else
    disp('Moog Crash wasn''t found, make sure about reason of termination of Spike2')    
    end
else
    moog_crash = [];
end
%% define gaussian filter (done)
% use Median instead of Gaussian filter? No
sig = prs.filtwidth; %filter width
sz = prs.filtsize; %filter size
t2 = linspace(-sz/2, sz/2, sz);
h = exp(-t2.^2/(2*sig^2));
h = h/sum(h); % normalise filter to ensure area under the graph of the data is not altered
%% define Butterworth filter
SR = 1/dt;
[b,a] = butter(2,30/(SR/2),'low'); % for the sensors
%% filter position, velocity and sensor channels (done)
screwed = []; % keep track of screwed trials
for i=1:length(trl)
    %     if i == 13 % use this to fix temporarily wrong trials, if any
    %         keyboard;
    %     end
    for j = 1:length(chnames)
        if ~any(strcmp(chnames{j},{'tsi','ts','mrk','firefly','yle','yre','zle','zre'}))
            if isempty(trl(i).continuous.(chnames{j}))
                screwed = [screwed i];
                continue;
            end
            if any(strcmp(chnames{j},{'xac','yac','vrol','vyaw','vpit'}))
                trl(i).continuous.(chnames{j}) = filtfilt(b,a,trl(i).continuous.(chnames{j})); % butterowrth filter for the sensors
            else
                %           trl(i).continuous.(chnames{j}) = medfilt1(trl(i).continuous.(chnames{j}),5); % median filter
                trl(i).continuous.(chnames{j}) = [trl(i).continuous.(chnames{j})(1)*ones(1,100)';...
                    trl(i).continuous.(chnames{j});...
                    trl(i).continuous.(chnames{j})(end)*ones(1,100)']; % adding first and last value in edges
                trl(i).continuous.(chnames{j}) = conv(trl(i).continuous.(chnames{j})(1:end),h,'same'); % gaussian filter for data
                %           ch.(chnames{i}) = ch.(chnames{i})(sz/2+1:end);
                trl(i).continuous.(chnames{j}) = trl(i).continuous.(chnames{j})(1+100:end-100);
            end
        end
    end
end
if screwed
    screwed = unique(screwed);
    disp(['There are screwed trials No. ' num2str(screwed) '!!!!! Check AddSMRData'])
    keyboard;
end
%% Detect saccades
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
%% timestamps referenced relative to exp_beg
% exp_beg = t.events(1);
% exp_end = t.events(end);
% 
% for i=1:length(trl)
%     trl(i).events.t_beg = trl(i).events.t_beg - exp_beg;
%     trl(i).events.t_end = trl(i).events.t_end - exp_beg - trl(i).events.t_beg;    
%     trl(i).events.t_sac = trl(i).events.t_sac - exp_beg - trl(i).events.t_beg;
%     trl(i).events.t_move = trl(i).events.t_move - exp_beg - trl(i).events.t_beg;
%     trl(i).events.t_stop = trl(i).events.t_stop - exp_beg - trl(i).events.t_beg;
%     trl(i).events.t_ptb = trl(i).events.t_ptb - exp_beg - trl(i).events.t_beg;
%     trl(i).events.t_targ = 0;
% end

%% downsample continuous data
% for i=1:length(chnames)
%     if ~any(strcmp(chnames{i},'mrk'))
%         ch.(chnames{i}) = ch.(chnames{i})(ts>exp_beg & ts<exp_end);
%         ch.(chnames{i}) = downsample(ch.(chnames{i}),prs.factor_downsample);
%     end
% end
% ts = ts(ts>exp_beg & ts<exp_end) - exp_beg;
% ch.ts = downsample(ts,prs.factor_downsample); ch.ts = ch.ts(:);
ch.ntrls = length(trl);
ch.dt = dt;

%% debugging plots

% figure;plot(ch.ts,ch.v)
% hold on;plot(ch.ts,ch.tsi*100)
% vline(t.beg)
% vline(t.end,'k')

