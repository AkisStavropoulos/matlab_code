function [trials,totaltrls] = AddLOGData_MKNOMC_EDOARDO(file)

count = 0;
fid = fopen(file, 'r');
eof=0; newline = 'nothingnew'; count=0;
% newline has to be up to 10 chars, because of a short word in the file early on
%% check if this data was generated by replaying the stimulus movie
% if contains(file,'replay'), replay_movie = true;
% else, replay_movie = false; end
%% check for fixed ground landmark
% while ~contains(newline,'Enable Lifetime')
%     newline = fgetl(fid);
% end
% fixed_ground = logical(1 - str2double(newline(18))); % if limited lifetime is enabled (1), fixed_ground is 0
%% check if this data was generated by replaying the stimulus movie
if strcmp(file(1:6),'replay'), replay_movie = true;
else, replay_movie = false; end

%% Get Eye offset and scale info
while ~contains(newline,'Monkey Name')
    newline = fgetl(fid);
end
if contains(newline,'Monkey Name')
    subj_name = split(newline,': ');
    subj_name = subj_name{2};
    newline = fgetl(fid);
    
end
while ~contains(newline,'DataFileName')
    newline = fgetl(fid);
end
if contains(newline,'DataFileName')
    file_name = split(newline,': ');
    file_name = file_name{2};
    newline = fgetl(fid);
end
while ~contains(newline,'Interocular')
    newline = fgetl(fid);
end
if contains(newline,'Interocular')
    eye_dist = split(newline,': ');
    eye_dist = str2double(eye_dist{2});
    newline = fgetl(fid);
end
% Left Eye
if contains(newline,'LeftEye Offset X')
    yle_offset = split(newline,': ');
    yle_offset = str2double(yle_offset{2});
    newline = fgetl(fid);
end
if contains(newline,'LeftEye Scale X')
    yle_scale = split(newline,': ');
    yle_scale = str2double(yle_scale{2});
    newline = fgetl(fid);
end
if contains(newline,'LeftEye Offset Y')
    zle_offset = split(newline,': ');
    zle_offset = str2double(zle_offset{2});
    newline = fgetl(fid);
end
if contains(newline,'LeftEye Scale Y')
    zle_scale = split(newline,': ');
    zle_scale = str2double(zle_scale{2});
    newline = fgetl(fid);
end
% Right Eye
if contains(newline,'RightEye Offset X')
    yre_offset = split(newline,': ');
    yre_offset = str2double(yre_offset{2});
    newline = fgetl(fid);
end
if contains(newline,'RightEye Scale X')
    yre_scale = split(newline,': ');
    yre_scale = str2double(yre_scale{2});
    newline = fgetl(fid);
end
if contains(newline,'RightEye Offset Y')
    zre_offset = split(newline,': ');
    zre_offset = str2double(zre_offset{2});
    newline = fgetl(fid);
end
if contains(newline,'RightEye Scale Y')
    zre_scale = split(newline,': ');
    zre_scale = str2double(zre_scale{2});
    newline = fgetl(fid);
end
%% Get head, eye and screen distance info
if contains(newline,'Monkey to Screen')
    screen_dist = split(newline,': ');
    screen_dist = str2double(screen_dist{2});
    newline = fgetl(fid);
end
if contains(newline,'Head Center X')
    headX = split(newline,': ');
    headX = str2double(headX{2});
    newline = fgetl(fid);
end
if contains(newline,'Head Center Y')
    headY = split(newline,': ');
    headY = str2double(headY{2});
    newline = fgetl(fid);
end
if contains(newline,'Head Center Z')
    headZ = split(newline,': ');
    headZ = str2double(headZ{2});
    newline = fgetl(fid);
end
if contains(newline,'Eye Offset X')
    eyeX = split(newline,': ');
    eyeX = str2double(eyeX{2});
    newline = fgetl(fid);
end
if contains(newline,'Eye Offset Y')
    eyeY = split(newline,': ');
    eyeY = str2double(eyeY{2});
    newline = fgetl(fid);
end
if contains(newline,'Eye Offset Z')
    eyeZ = split(newline,': ');
    eyeZ = str2double(eyeZ{2});
    newline = fgetl(fid);
end
while ~contains(newline,'Floor Origin Z')
    newline = fgetl(fid);
end
if contains(newline,'Floor Origin Z')
    viewHeight = split(newline,': ');
    viewHeight = str2double(viewHeight{2});
    newline = fgetl(fid);
end
while ~contains(newline,'Floor Density manipulation') && ~contains(newline,'Reward Bou')
    newline = fgetl(fid);
end
if contains(newline,'Floor Density manipulation')
    density_manipulation = split(newline,': ');
    density_manipulation = str2double(density_manipulation{2});
    newline = fgetl(fid);
else
    density_manipulation = 0;
end


%% Get firefly range & reward boundary info
while ~contains(newline,'Reward Bou')
    newline = fgetl(fid);
end
if contains(newline,'Reward Bou')
    rew_boundary = split(newline,': ');
    rew_boundary = str2double(rew_boundary{2});
    newline = fgetl(fid);
end
if contains(newline,'Linear Rew')
    lin_reward = split(newline,': ');
    lin_reward = str2double(lin_reward{2});
    newline = fgetl(fid);
end
if contains(newline,'Firefly Rad')
    ff_radius = split(newline,': ');
    ff_radius = str2double(ff_radius{2});
    newline = fgetl(fid);
end
if contains(newline,'Firefly On Dur')
    ff_dur = split(newline,': ');
    ff_dur = str2double(ff_dur{2});
    newline = fgetl(fid);
end
if contains(newline,'Firefly Range Radius')
    ff_dist_range = split(newline,': ');
    ff_dist_range = str2num(ff_dist_range{2});
    newline = fgetl(fid);
end
if contains(newline,'Firefly Range Angle ')
    ff_ang_range = split(newline,': ');
    ff_ang_range = str2num(ff_ang_range{2});
    newline = fgetl(fid);
end
%% Get Tau dynamics info
while ~contains(newline,'Discrete')
    newline = fgetl(fid);
end
while ~contains(newline,'Discrete(0) or'); newline = fgetl(fid); end
if contains(newline,'Discrete(0) or')
    tau_dyn_temp = split(newline,': ');
    tau_dyn_temp = str2double(tau_dyn_temp{2});
    if tau_dyn_temp
        tau_dynamics = 'continuous';
    else
        tau_dynamics = 'discrete';
    end
    newline = fgetl(fid);
end
if contains(newline,'Trial(0) or')
    timestep_temp = split(newline,': ');
    timestep_temp = str2double(timestep_temp{2});
    if timestep_temp
        timestep = 'time';
    else
        timestep = 'trial';
    end
    newline = fgetl(fid);
end
if contains(newline,'Tau Tau')
    tau_tau = split(newline,': ');
    tau_tau = str2double(tau_tau{2});
    newline = fgetl(fid);
end
if contains(newline,'Mean Firefly Dist')
    x = split(newline,': ');
    x = str2double(x{2});
    newline = fgetl(fid);
end
if contains(newline,'Mean Trial Dur')
    T = split(newline,': ');
    T = str2double(T{2});
    newline = fgetl(fid);
end
if contains(newline,'Gain of Max Ang')
    gain_wmax = split(newline,': ');
    gain_wmax = str2double(gain_wmax{2});
    newline = fgetl(fid);
else
    gain_wmax = 1;
end
if contains(newline,'Number of Taus')
    if strcmp(tau_dynamics,'discrete')
        num_taus = split(newline,': ');
        num_taus = str2double(num_taus{2});
    else
        num_taus = [];
    end
    newline = fgetl(fid);
end
if contains(newline,'MinTau')
    min_tau = split(newline,': ');
    min_tau = str2double(min_tau{2});
    newline = fgetl(fid);
end
if contains(newline,'MaxTau')
    max_tau = split(newline,': ');
    max_tau = str2double(max_tau{2});
    newline = fgetl(fid);
end
if contains(newline,'Gamma')
    gamma = split(newline,': ');
    gamma = str2double(gamma{2}); % check what gamma actually is
    newline = fgetl(fid);
end
if strcmp(tau_dynamics,'continuous')
    if contains(newline,'MuPhi')
        mu_phi = split(newline,': ');
        mu_phi = str2double(mu_phi{2});
        newline = fgetl(fid);
    end
    if contains(newline,'SigmaPhi')
        sig_phi = split(newline,': ');
        sig_phi = str2double(sig_phi{2});
        newline = fgetl(fid);
    end
    if contains(newline,'MuEta')
        mu_eta = split(newline,': ');
        mu_eta = str2double(mu_eta{2});
        newline = fgetl(fid);
    end
    
    if contains(newline,'SigmaEta')
        sig_eta = split(newline,': ');
        sig_eta = str2double(sig_eta{2});
        newline = fgetl(fid);
    end
else
    mu_phi = [];
    sig_phi = [];
    mu_eta = [];
    sig_eta = [];
    while ~contains(newline,'Trial Num#'); newline = fgetl(fid); end
end
while newline ~= -1
    %% get ground plane density, firefly position, stimulus type, vmax, wmax, tau
    while ~contains(newline,{'Trial Num#','Floor Dens','Firefly Po','Stim Type:','Joy Stick '})
        newline = fgetl(fid);
        if newline == -1, break; end
    end
    if newline == -1, break; end
    count = count+1;
    if contains(newline,'Trial Num#')
        totaltrls = split(newline,'# ');
        totaltrls = str2double(totaltrls{2});
        newline = fgetl(fid);
    end
    if contains(newline,'Floor Dens')
        floordensity = split(newline,': ');
        floordensity = str2double(floordensity{2});
        newline = fgetl(fid);
    end
    if contains(newline,'Inter-tria')
        iti = split(newline,': ');
        iti = str2double(iti{2});
        newline = fgetl(fid);
    end
    if contains(newline,'Firefly Full')
        FF_fullON = split(newline,': ');
        FF_fullON = str2double(FF_fullON{2});
        newline = fgetl(fid);
    else
        FF_fullON = 0;
    end
    if contains(newline,'Enable Lif')
        lifetime = split(newline,': ');
        lifetime = str2double(lifetime{2});
        newline = fgetl(fid);
    end
    if contains(newline,'Caught Fir')
        newline = fgetl(fid);
    end
    if contains(newline,'Reward Dur')
        rew_dur = split(newline,': ');
        rew_dur = str2double(rew_dur{2});
        newline = fgetl(fid);
    end
    if contains(newline,'Firefly')
        keyword = split(newline,': ');
        keyword = keyword{2};
        keyword = split(keyword,' ');
        temp1 = keyword{1};
        temp2 = keyword{2};
        trials(count).prs.xfp = str2double(temp1);
        trials(count).prs.yfp = str2double(temp2);
        newline = fgetl(fid);
    end
    if contains(newline,'Stim Type:')
        %         trials(count).prs.stimtype = str2double(newline(12));
        trials(count).prs.stimtype = nan;
        newline = fgetl(fid);
    end
    if contains(newline,'Joy Stick ')
        if contains(newline,'Joy Stick Tau')
            tau = split(newline,': ');
            tau = str2double(tau{2});
            newline = fgetl(fid);
        end
        if contains(newline,'Joy Stick Max Vel')
            vmax = split(newline,': ');
            vmax = str2double(vmax{2});
            newline = fgetl(fid);
        end
        if contains(newline,'Joy Stick Max Ang')
            wmax = split(newline,': ');
            wmax = str2double(wmax{2});
            newline = fgetl(fid);
        end
    end
    
    trials(count).prs.filename = file_name;
    trials(count).prs.eye_dist = eye_dist;
    trials(count).prs.yle_offset = yle_offset;  trials(count).prs.yle_scale = yle_scale;
    trials(count).prs.yre_offset = yre_offset;  trials(count).prs.yre_scale = yre_scale;
    trials(count).prs.zle_offset = zle_offset;  trials(count).prs.zle_scale = zle_scale;
    trials(count).prs.zre_offset = zre_offset;  trials(count).prs.zre_scale = zre_scale;
    
    trials(count).prs.screen_dist = screen_dist;
    trials(count).prs.headX = headX;
    trials(count).prs.headY = headY;
    trials(count).prs.headZ = headZ;
    trials(count).prs.eyeX = eyeX;
    trials(count).prs.eyeY = eyeY;
    trials(count).prs.eyeZ = eyeZ;
    trials(count).prs.viewHeight = viewHeight;
    
    trials(count).prs.floorlifetime = lifetime;
    trials(count).prs.floordensity = floordensity;
    trials(count).prs.density_manipulation = density_manipulation;
    trials(count).prs.intertrial_interval = iti;
    trials(count).prs.reward_boundary = rew_boundary;
    trials(count).prs.reward_duration = rew_dur;
    trials(count).prs.linear_reward = lin_reward;
    trials(count).prs.fly_radius = ff_radius;
    trials(count).prs.fly_duration = ff_dur;
    trials(count).prs.fly_dist_range = ff_dist_range;
    trials(count).prs.fly_ang_range = ff_ang_range;
    
    trials(count).logical.firefly_fullON = FF_fullON;
    trials(count).logical.replay = replay_movie;
    trials(count).logical.landmark_distance = false;
    trials(count).logical.landmark_angle = false;
    trials(count).logical.landmark_fixedground = false;
    trials(count).logical.joystick_gain = 1;
    trials(count).logical.ptb = 0;
    trials(count).logical.microstim = 0;
    trials(count).prs.ptb_linear = 0;
    trials(count).prs.ptb_angular = 0;
    trials(count).prs.ptb_delay = 0;
    
    
    trials(count).prs.tau = tau;
    trials(count).prs.v_max = vmax;
    trials(count).prs.w_max = wmax;
    trials(count).prs.x = x;
    trials(count).prs.T = T;
    trials(count).prs.gain_wmax = gain_wmax;
    
    trials(count).prs.tau_dynamics = tau_dynamics;
    trials(count).prs.num_taus = num_taus;
    trials(count).prs.timestep = timestep;
    trials(count).prs.tau_tau = tau_tau;
    trials(count).prs.gamma = gamma;
    trials(count).prs.min_tau = min_tau;
    trials(count).prs.max_tau = max_tau;
    trials(count).prs.mu_phi = mu_phi;
    trials(count).prs.sig_phi = sig_phi;
    trials(count).prs.mu_eta = mu_eta;
    trials(count).prs.sig_eta = sig_eta;
    
    % give subject name
    if ~isempty(subj_name)
        trials(count).prs.subject = [subj_name];
    else
        disp(['Subject name not obtained, check current folder name and AddTrials2Behaviour.m']);
    end
    
    if newline == -1, break; end
    % initialise
    %     trials(count).logical.landmark_distance = false;
    %     trials(count).logical.landmark_angle = false; % #$%^&&^&*^danger - change false to nan immediately (what if field missing from log file??)
    %     trials(count).logical.landmark_fixedground = fixed_ground;
    %     trials(count).prs.ptb_linear = 0;
    %     trials(count).prs.ptb_angular = 0;
    %     trials(count).prs.ptb_delay = 0;
    %     trials(count).prs.intertrial_interval = nan;
    %     trials(count).logical.firefly_fullON = nan;
    %     trials(count).prs.stop_duration = nan;
    %     trials(count).logical.replay = replay_movie;
    %% get landmark status, ptb velocities and ptb delay
    %     newline = fgetl(fid);
    %     if newline == -1, break; end
    %     if contains(newline(1:9),'Enable Di')
    %         trials(count).logical.landmark_distance = str2double(newline(26)); % 1=distance landmark was ON
    %         newline = fgetl(fid);
    %         trials(count).logical.landmark_angle = str2double(newline(25)); % 1=angular landmark was ON
    %         newline = fgetl(fid);
    %         trials(count).prs.ptb_linear = str2double(newline(35:end)); % amplitude of linear velocity ptb (cm/s)
    %         newline = fgetl(fid);
    %         trials(count).prs.ptb_angular = str2double(newline(37:end)); % amplitude of angular velocity ptb (deg/s)
    %         newline = fgetl(fid);
    %         trials(count).prs.ptb_delay = str2double(newline(31:end)); % time after trial onset at which to begin ptb
    %         newline = fgetl(fid);
    %     end
    %% get inter-trial interval and firefly status
    %     if newline == -1, break; end
    %     if contains(newline(1:10),'Firefly ON')
    %         trials(count).prs.intertrial_interval = str2double(newline(27:end)); % time between end of this trial and beg of next trial (s)
    %         newline = fgetl(fid);
    %         trials(count).logical.firefly_fullON = str2double(newline(18)); % 1=firefly was ON throughout the trial
    %         newline = fgetl(fid);
    %     end
    %% get stopping duration for reward
    %     if newline == -1, break; end
    %     if contains(newline(1:8),'Distance')
    %         trials(count).prs.stop_duration = str2double(newline(34:end))/1000; % wait duration after stopping before monkey is given feedback (s)
    %     end
end

fclose(fid);