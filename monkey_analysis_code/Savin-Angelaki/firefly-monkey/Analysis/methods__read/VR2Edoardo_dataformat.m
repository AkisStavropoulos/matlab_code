function trials = VR2Edoardo_dataformat(cont_data,disc_data,prs)

disp('converting data to old format...')

Ntrials = disc_data.trial_num(end);
dt = prs.dt;
scale = 100;
presac_t = prs.presac_t;
postsac_t = prs.postsac_t;
ts = cont_data.trial_time - cont_data.trial_time(1);

%% Detect saccades and fixations

if isfield(prs,'humanVR')
   LXnorm  = (cont_data.Gx - cont_data.Gx0);
   LYnorm  = (cont_data.Gy - cont_data.Gy0);
   LZnorm  = (cont_data.Gz - cont_data.Gz0);
   vec_mag = sqrt(LXnorm.^2 + LYnorm.^2 + LZnorm.^2);
   
   LXnorm = LXnorm./vec_mag;    RXnorm = LXnorm;
   LYnorm = LYnorm./vec_mag;    RYnorm = LYnorm;
   LZnorm = LZnorm./vec_mag;    RZnorm = LZnorm;
else
    LXnorm = cont_data.LXnorm;  RXnorm = cont_data.RXnorm;
    LYnorm = cont_data.LYnorm;  RYnorm = cont_data.RYnorm;
    LZnorm = cont_data.LZnorm;  RZnorm = cont_data.RZnorm;
    
end

    % detect saccade times
    yle = atan2d(LXnorm,LZnorm);  % double check
    zle = atan2d(LYnorm,LZnorm); 
    yre = atan2d(RXnorm,RZnorm); 
    zre = atan2d(RYnorm,RZnorm); 

    % take derivative of eye position = eye velocity
    if all(prs.eyechannels ~= 0)
        dze = diff(0.5*(zle + zre));
        dye = diff(0.5*(yle + yre));
    elseif prs.eyechannels(1) ~= 0
        dze = diff(zle);
        dye = diff(yle);
    else
        dze = diff(zre);
        dye = diff(yre);
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
    
    % apply threshold on eye speed
    saccade_thresh = prs.saccade_thresh;
    indx_thresh = de_smooth>saccade_thresh;
    dindx_thresh = diff(indx_thresh);
    t_saccade = find(dindx_thresh>0)*dt;
    
    % remove duplicates by applying a saccade refractory period
    min_isi = prs.min_intersaccade;
    t_saccade(diff(t_saccade)<min_isi) = [];
    t.saccade = t_saccade;
    
    % detect fixations
    fixateduration = prs.fixateduration; fixate_thresh = prctile(de_smooth,90); % set thresh to 90th prctile
    fixateduration_samples = round(fixateduration/dt);
    fixateindx = false(1,numel(ts) - round(2*fixateduration/dt));
    for j=1:(numel(ts) - round(2*fixateduration/dt))
        if mean(de_smooth(j:j+fixateduration_samples)) < fixate_thresh && max(de_smooth(j:j+fixateduration_samples)) < 1.5*fixate_thresh, fixateindx(j) = true; end
    end
    fixation_switch = diff(fixateindx);
    t.fix = ts(fixation_switch>0);

%% Set up data
for i = 1:Ntrials

    % detect start-of-movement
    v_thresh = 0.05; % m/s
    w_thresh = 3; % deg/s
    v_time2thresh = 0.05; % (s) approx time to go from zero to threshold or vice-versa
    v = cont_data.linear_velocity;
    if i==1, t_move(i) = disc_data.beginTime(i); % first trial is special because there is no pre-trial period
    else
        indx = find(v(ts>disc_data.endTime(i-1) & ts<disc_data.endTime(i)) > v_thresh,1); % first upward threshold-crossing
        if ~isempty(indx), t_move(i) = disc_data.endTime(i-1) + indx*dt;
        else, t_move(i) = disc_data.beginTime(i); end % if monkey never moved, set movement onset = target onset
    end
    
    % define pre-/post-trial period
    pretrial = max(disc_data.beginTime(i) - t_move(i),0) + prs.pretrial; % extract everything from "movement onset - pretrial" or "target onset - pretrial" - whichever is first
    posttrial = prs.posttrial; % extract everything until "t_end + posttrial"

    % select samples
    timeindx = find(cont_data.trial_time > disc_data.beginTime(i)-pretrial & cont_data.trial_time < disc_data.endTime(i)+posttrial); 
    beg_indx = find(cont_data.trial_time > disc_data.beginTime(i),1); % timeindx = timeindx(1:end-2); % cont_data.start_trial(i):cont_data.stop_trial(i);
    
    %% continuous
    trials(i).continuous.yle = yle(timeindx);
    trials(i).continuous.zle = zle(timeindx); 
    trials(i).continuous.yre = yre(timeindx); 
    trials(i).continuous.zre = zre(timeindx); 
    trials(i).continuous.ybe = 0.5*(trials(i).continuous.yle + trials(i).continuous.yre);
    trials(i).continuous.zbe = 0.5*(trials(i).continuous.zle + trials(i).continuous.zre);
    trials(i).continuous.xfp = disc_data.FFx(i)*ones(numel(timeindx),1)*scale; % same X as Jing? Is it relative to subject or world?
    trials(i).continuous.yfp = disc_data.FFz(i)*ones(numel(timeindx),1)*scale; % Z is Jing's Y
    trials(i).continuous.xmp = cont_data.posX(timeindx)*scale; % in meters!!!
    trials(i).continuous.ymp = cont_data.posZ(timeindx)*scale;    
    trials(i).continuous.v = cont_data.linear_velocity(timeindx)*scale;
    trials(i).continuous.w = cont_data.angular_velocity(timeindx);
    if isfield(prs,'humanVR') || isfield(cont_data,'clean_lin_vel')
        trials(i).continuous.v_clean = cont_data.clean_lin_vel(timeindx)*scale;
        trials(i).continuous.w_clean = cont_data.clean_ang_vel(timeindx);
        trials(i).continuous.v_ksi = cont_data.ksi_lin_vel(timeindx)*scale;
        trials(i).continuous.w_ksi = cont_data.ksi_ang_vel(timeindx);
        trials(i).continuous.v_eta = cont_data.eta_lin_vel(timeindx)*scale;
        trials(i).continuous.w_eta = cont_data.eta_ang_vel(timeindx);
        trials(i).continuous.yjs = cont_data.raw_lin_js(timeindx);
        trials(i).continuous.xjs = cont_data.raw_ang_js(timeindx);
    end
    trials(i).continuous.FFdraw = 5*(cont_data.trial_time(timeindx)-cont_data.trial_time(beg_indx) < disc_data.ff_duration(i) ); 
    trials(i).continuous.h1 = nan(numel(timeindx),1); % hand position
    trials(i).continuous.h2 = nan(numel(timeindx),1);
    trials(i).continuous.ts = cont_data.trial_time(timeindx)-cont_data.trial_time(beg_indx);
    trials(i).continuous.firefly = trials(i).continuous.ts>=0 & trials(i).continuous.ts<(0+disc_data.ff_duration(i)); % cont_data.trial_time(timeindx)-cont_data.trial_time(timeindx(1)) < disc_data.ff_duration(i);

    % corrected monkey position
    trials(i).continuous.xmp(end-floor(posttrial/dt)-2:end) = trials(i).continuous.xmp(end-floor(posttrial/dt)-3);    
    trials(i).continuous.ymp(end-floor(posttrial/dt)-2:end) = trials(i).continuous.ymp(end-floor(posttrial/dt)-3);    
    
    % correct target position at end of trial
    trials(i).continuous.xfp(end-floor(posttrial/dt)-3:end) = trials(i).continuous.xfp(end-floor(posttrial/dt)-4);
    trials(i).continuous.yfp(end-floor(posttrial/dt)-3:end) = trials(i).continuous.yfp(end-floor(posttrial/dt)-4);

    
    % set position values prior to target onset to nan
    chnames = fieldnames(trials(i).continuous);
    for n = 1:length(chnames)
        if any(strcmp(chnames{n},{'xfp','xmp','yfp','ymp','xmp1','ymp1','phi'})) % set position values prior to target onset to nan
            trials(i).continuous.(chnames{n})(1:floor(pretrial/dt)) = nan;
        end
    end
    
    %% events
    exp_beg = 0; % first marker that aligns the neural data, find out which!
    trials(i).events.t_beg = disc_data.beginTime(i) - exp_beg;
    trials(i).events.t_end = disc_data.endTime(i) - exp_beg - trials(i).events.t_beg;
    trials(i).events.t_move = t_move(i) - exp_beg - trials(i).events.t_beg;
    trials(i).events.t_stop = disc_data.checkTime(i) - exp_beg - trials(i).events.t_beg;
    trials(i).events.t_sac = t.saccade(t.saccade>(disc_data.beginTime(i)-pretrial) & t.saccade<disc_data.endTime(i))- exp_beg - trials(i).events.t_beg;   trials(i).events.t_sac = trials(i).events.t_sac(:)';
    trials(i).events.t_fix = t.fix(t.fix>(disc_data.beginTime(i)-3*pretrial) & t.fix<(disc_data.endTime(i)+3*posttrial)) - exp_beg - trials(i).events.t_beg; % wider search-range because fixation might be outside trial 
    t.fix(t.fix>(disc_data.beginTime(i)-3*pretrial) & t.fix<(disc_data.endTime(i)+3*posttrial)) = []; % remove from list
    
    if isfield(disc_data,'rewardTime')
    if disc_data.rewardTime(i) > 0
        trials(i).events.t_rew = disc_data.rewardTime(i) - exp_beg - trials(i).events.t_beg; % need to add or find out time of reward
    else; trials(i).events.t_rew = 0;
    end
    end
    
    try 
        if isfield(disc_data,'ptbSigma') % initial format Mehran implemented
        trials(i).events.t_ptb = disc_data.ptbStart(i); 
        else % final format consistent with other timestamps
        trials(i).events.t_ptb = disc_data.ptbStart(i) - exp_beg - trials(i).events.t_beg;
        end
    catch; trials(i).events.t_ptb = nan; 
    end
    try trials(i).events.t_microstim = disc_data.microstimTime(i); catch; trials(i).events.t_microstim = nan; end
    trials(i).events.t_targ = 0;
    trials(i).events.t_beg_correction = 0;
    trials(i).events.t_flyON = 0; % corrected time of FF appearance
    trials(i).events.t_flyON_minus_teleport = nan;
    
    trials(i) = get_sac_mag_dir(trials(i),presac_t,postsac_t);
    
    %% logical
    trials(i).logical.landmark_distance = false;
    trials(i).logical.landmark_angle = false;
    trials(i).logical.firefly_fullON = disc_data.ff_duration(i) > 0.4; % add binary column at the end for this
    trials(i).logical.replay = false;
    trials(i).logical.landmark_fixedground = false;
    trials(i).logical.joystick_gain = 1;
    trials(i).logical.reward = logical(disc_data.rewarded(i));
    try trials(i).logical.ptb = logical(disc_data.ptbTrial(i)); catch; trials(i).logical.ptb = false; end
    try trials(i).logical.microstim = logical(disc.microstimTrial(i)); catch; trials(i).logical.microstim = false; end
    trials(i).logical.spurioustarg = false;
    
    
    %% prs
    trials(i).prs.floordensity = disc_data.floordensity(i)*0.0001; 
    try trials(i).prs.ptb_linear = disc_data.ptbLinVel(i)*scale; catch; trials(i).prs.ptb_linear = 0; end
    try trials(i).prs.ptb_angular = disc_data.ptbAngVel(i); catch; trials(i).prs.ptb_angular = 0; end
    trials(i).prs.ptb_delay = 0;
    trials(i).prs.intertrial_interval = disc_data.ITI(i);
    trials(i).prs.stop_duration = 0; % distance checking interval (waitTime) but it's not used
    trials(i).prs.vmax = disc_data.maxV(i)*scale;
    trials(i).prs.wmax = disc_data.maxW(i);
    if isfield(disc_data,'rewardDur')
    trials(i).prs.reward_duration = disc_data.rewardDur(i); % include this
    end
    trials(i).prs.xfp = disc_data.FFx(i)*scale;
    trials(i).prs.yfp = disc_data.FFz(i)*scale;
    trials(i).prs.fly_duration = disc_data.ff_duration(i);
    
    if isfield(prs,'humanVR') || isfield(disc_data,'tau')
        trials(i).prs.tau = disc_data.tau(i);
        trials(i).prs.type = disc_data.type(i); % 0) discrete, 1)continuous , 2) none
        if trials(i).prs.type
            trials(i).prs.ntaus = nan;
        else
            trials(i).prs.ntaus = disc_data.ntaus(i);
        end
        trials(i).prs.mintau = disc_data.mintau(i);
        trials(i).prs.maxtau = disc_data.maxtau(i);
        trials(i).prs.x = disc_data.x(i);
        trials(i).prs.T = disc_data.T(i);
        trials(i).prs.vthresh = disc_data.vthresh(i);
        trials(i).prs.wthresh = disc_data.wthresh(i);
    end
    
    if isfield(prs,'humanVR') || isfield(disc_data,'noisetau')        
        trials(i).prs.noisetau = disc_data.noisetau(i);
        trials(i).prs.noisetautau = disc_data.noisetautau(i);
        trials(i).prs.gainw = disc_data.gainw(i);
        trials(i).prs.gainv = disc_data.gainv(i);
    end
    
end

disp('format converted')


if 0
    
    
    
    
    
end
