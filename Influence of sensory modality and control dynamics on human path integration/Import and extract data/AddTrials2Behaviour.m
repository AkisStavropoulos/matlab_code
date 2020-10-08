function trials = AddTrials2Behaviour(prs)
% default_prs; RUN IF NEEDED
trials = []; % initialise
% DATA MUST BE IN THE CURRENT FOLDER !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%% list all files to read
flist_log=dir('*.log'); 
for j=1:length(flist_log)
    fnum_log(j) = str2num(flist_log(j).name(end-6:end-4)); 
end
flist_smr=dir('*.smr');
for j=1:length(flist_smr)
    fnum_smr(j) = str2num(flist_smr(j).name(end-6:end-4)); 
end
flist_cal = dir('eye calibration'); flist_cal(1:2) = [];    count = 1;  indx = [];
for j = 1:length(flist_cal)
   if strcmp(flist_cal(j).name(end-3:end),'.SMR') 
       indx = [indx j];
       fnum_cal(count) = str2num(flist_cal(j).name(end-6:end-4));
       count = count + 1;
   end
end
flist_cal = flist_cal(indx);
flist_mc=dir('*MC_Variables*');
fnum_mc = [];
for j=1:length(flist_mc)
    fnum_mc(j) = str2num(flist_mc(j).name(end-2:end)); 
end

% flist_mat=dir('*.mat');
% for i=1:length(flist_mat)
%     fnum_mat(i) = str2num(flist_mat(i).name(end-6:end-4)); 
% end
nfiles = length(flist_log);
%% read files
for j=1:nfiles
    fprintf(['... reading ' flist_log(j).name '\n']);
    % read .log file
    trials_log = AddLOGData(flist_log(j).name);
    % read all .smr files associated with this log file
    if j<nfiles 
        indx_smr = find(fnum_smr >= fnum_log(j) & fnum_smr < fnum_log(j+1));
    else
        indx_smr = find(fnum_smr >= fnum_log(j));
    end
    trials_smr = [];
    moog_crash = [];
    for j = indx_smr
        data_smr = ImportSMR(flist_smr(j).name);
        % find corresponding eye calibration file
        if j > 1
            ind = find((fnum_cal < fnum_smr(j)) .* (fnum_cal > fnum_smr(j-1)),1,'last');
        else
            ind = find(fnum_cal < fnum_smr(j),1,'last');
        end
        EyeCalFix;
        if ~isempty(ind)
            cd([pwd '\eye calibration']) % extract from eye calibration folder
            data_cal = ImportSMR(flist_cal(ind).name);
            disp(['Block: ' flist_smr(j).name ', EyeCal: ' flist_cal(ind).name]);
            curfol = pwd;
            cd(curfol(1:end-15)) % go back to subject folder
        else
            data_cal = [];
            disp(['Block: ' flist_smr(j).name ', EyeCal: -']);
        end
        % exporteye parameters
        gains = [];
        gains.eye_dist = trials_log(1).prs.eye_dist;    gains.screen_dist = trials_log(1).prs.screen_dist;
        gains.headY = trials_log(1).prs.headY;          gains.eyeY = trials_log(1).prs.eyeY;
        gains.yle_offset = trials_log(1).prs.yle_offset;    gains.yle_scale = trials_log(1).prs.yle_scale;
        gains.yre_offset = trials_log(1).prs.yre_offset;    gains.yre_scale = trials_log(1).prs.yre_scale;
        gains.zle_offset = trials_log(1).prs.zle_offset;    gains.zle_scale = trials_log(1).prs.zle_scale;
        gains.zre_offset = trials_log(1).prs.zre_offset;    gains.zre_scale = trials_log(1).prs.zre_scale;
        
        [trials_smr_temp,moog_crash_temp,~,t] = AddSMRData(data_smr,data_cal,gains,prs);
        trials_smr = [trials_smr trials_smr_temp];
        moog_crash = [moog_crash moog_crash_temp];
    end
    % read all MC_Variables files associated with this log file
    if j<nfiles
        indx_mc = find(fnum_mc >= fnum_log(j) & fnum_mc < fnum_log(j+1));
    else
        indx_mc = find(fnum_mc >= fnum_log(j));
    end
    trials_mc = [];
    if ~isempty(indx_mc)
        for k = indx_mc
            data_mc = dlmread(flist_mc(j).name);
            trials_mc = [trials_mc AddMCData(data_mc,t)];
        end
    else % missing MC file
        disp(['MC file is MISSING from block ' flist_smr(j).name(end-6:end-4) ' !!!'])
        fnum_mc = [fnum_mc(1:j-1) NaN fnum_mc(j:end)]; % fix upcoming data names first
        for n = length(flist_mc):-1:j
            flist_mc(n+1) = flist_mc(n);
        end
        for j = indx_smr
            data_mc = [];
            trials_mc = [trials_mc AddMCData(data_mc,t)];
        end
    end

    % merge contents of .log, .smr and mc files
    ntrls_log = length(trials_log); ntrls_smr = length(trials_smr); ntrls_mc = length(trials_mc);
    if ntrls_smr ~= ntrls_mc
        disp('Number of trials in SMR and MC files are not the same!!!!')
    end
    if ntrls_smr <= ntrls_log
        for n=1:ntrls_smr
            trials_temp(n) = catstruct(trials_smr(n),trials_log(n),trials_mc(n)) ;
        end
    end
%         else  % apply a very dirty fix if spike2 was not "stopped" on time (can happen when replaying stimulus movie)
%             for j=1:ntrls_log
%                 trials_temp(j) = catstruct(trials_smr(j),trials_log(j),trials_mc(j)) ;
%             end
%             dummy_trials_log = trials_log(1:ntrls_smr-ntrls_log);
%             for j=1:(ntrls_smr-ntrls_log)
%                 trials_temp(ntrls_log+j) = catstruct(trials_smr(ntrls_log+j),dummy_trials_log(j));
%             end
%         end
        

    % add contents of .mat file
    %     trials_temp = AddMATData(flist_mat(i).name,trials_temp);
% remove crashed moog trials
if ~isempty(moog_crash)
    trials_temp(moog_crash) = [];
    disp(['Trials ' num2str(moog_crash) ' from file ' flist_smr(j).name ' were removed because of possible moog crash !'])
end
% move NOREC info from events to prs field
for n = 1:length(trials_temp)
    trials_temp(n).prs.NOREC = trials_temp(n).events.NOREC;
    trials_temp(n).events = rmfield( trials_temp(n).events,'NOREC');    
end
% concatenate trials
trials = [trials trials_temp];
clear trials_temp;
fprintf(['... total trials = ' num2str(length(trials)) '\n']);
end
%% create a channel for the firefly position in monkey coordinates
% for i = 1:length(trials)
%     theta = trials(i).continuous.phi;
%     trials(i).continuous.xfp = trials(i).continuous.xmp;
%     trials(i).continuous.xfp(~isnan(trials(i).continuous.xfp)) = 1;
%     trials(i).continuous.xfp = trials(i).continuous.xfp*trials(i).prs.fireflyposx;
%     trials(i).continuous.xfp = trials(i).continuous.xfp - trials(i).continuous.xmp; % translate x
% 
%     trials(i).continuous.yfp = trials(i).continuous.ymp;
%     trials(i).continuous.yfp(~isnan(trials(i).continuous.yfp)) = 1;
%     trials(i).continuous.yfp = trials(i).continuous.yfp*trials(i).prs.fireflyposy;
%     trials(i).continuous.yfp = trials(i).continuous.yfp - trials(i).continuous.ymp; % translate y
%     
%     fposxy = [trials(i).continuous.xfp trials(i).continuous.yfp];
% 
%     for j = 1:length(theta)
%             R = [cosd(-theta(j)) -sind(-theta(j)); sind(-theta(j)) cosd(-theta(j))]; % rotation matrix, bringing the targets in the fov, hence -theta
%         fpos = R*fposxy(j,:)';
%         trials(i).continuous.xfp = fpos(1);
%         trials(i).continuous.yfp = fpos(2);
%     end
% end
%% get events from MC file
for j = 1:length(trials)
    trials(j).events.target_off = trials(j).mc.target_off(trials(j).mc.target_off>0);
    trials(j).events.button = trials(j).mc.flag(trials(j).mc.flag>0);
    
    % put button recordings in events field
    trials(j).events.push_on = trials(j).mc.push_on;
    trials(j).events.push_off = trials(j).mc.push_off;

    %     trials(i).mc = rmfield(trials(i).mc,{'push_on','push_off'});
mc = trials(j).mc;
mc = rmfield(mc,{'target_off','flag','push_on','push_off'});
trials(j).mc = mc;

end
%% initialize position and velocity
% remove position jumps at trial onset
% put velocity to 0 before the subject starts pushing the joystick
% you can also use the v= at +  bu equation to fix the onset of velocity based on the JS input

% initialize position
plt = 0;
if plt;figure;hold on; end
for j = 1:length(trials)
    indx = find(trials(j).continuous.firefly,1,'last');
    startx = trials(j).continuous.xmp(indx-1);
    starty = trials(j).continuous.ymp(indx-1);
    
    trials(j).continuous.xmp(1:indx-1) = startx;
    trials(j).continuous.ymp(1:indx-1) = starty;
    
    if plt; plot(startx,starty,'.b'); end
end
if plt; axis equal; title(trials(1).prs.subject); vline(0); hline(0); end

% initialize velocity
immobile = [];
for j = 1:length(trials)
    
    d_tr = sqrt((trials(j).continuous.xmp(end) - trials(j).continuous.xmp(1)).^2 + (trials(j).continuous.ymp(end) - trials(j).continuous.ymp(1)).^2);
    move_on_ind = 1+ find(abs(diff(trials(j).mc.JS_X_Raw)) > 0,1);
    move_onset = trials(j).mc.timestamp(move_on_ind);
    
  
    if  (d_tr <= 30) % remove trials were subject didn't move / make sure there is recorded JS input?
        immobile = [immobile j];
        % check:
        % plot(trials(i).continuous.xmp,trials(i).continuous.ymp);ylim([-50 50]);axis equal;vline(0);hline(0);title([num2str(trials(i).prs.subject) ' - trial ' num2str(i)]);
        
    elseif move_onset
        indx = find(trials(j).continuous.ts >= move_onset,1);
        trials(j).continuous.v(1:indx) = 0;
        count = indx;
        reset = 1;
        while reset
            diff0 = abs(trials(j).continuous.v(count) - trials(j).continuous.v(count + 1));
            if diff0 > 2 % keep adding zeros until the hump is gone
                trials(j).continuous.v(count + 1) = 0;
                count = count + 1;
            else reset = 0;
            end
        end
    end


end
if ~isempty(immobile)
    disp(['Trials ' num2str(immobile) ' must be thrown out because the subject didn''t move.'])
    % remove
    trials(immobile) = [];
    disp(['Trials ' num2str(immobile) ' were thrown out because the subject didn''t move.'])
end

%% Downsample SMR file
datanames = fieldnames(trials(1).continuous);
S2dt = trials(1).continuous.ts(3) - trials(1).continuous.ts(2);
MCdt = trials(1).mc.timestamp(2) - trials(1).mc.timestamp(1);

for j  = 1:length(trials)
    trials(j).prs.dt = MCdt;
    % check length of smr files first
    for n = 1:length(datanames)
        SMRlength(n) = length(trials(j).continuous.(datanames{n}));
    end
    if length(unique(SMRlength)) > 1
        minSMRlength = min(SMRlength);
        for n = 1:length(datanames)
            trials(j).continuous.(datanames{n}) = trials(j).continuous.(datanames{n})(1:minSMRlength);
        end
    end
    % interpolate
    datalength = [];
    MCts = [];
    MCts = 0:MCdt:trials(j).continuous.ts(end);
    S2ts = [];
    S2ts = S2dt:S2dt:trials(j).continuous.ts(end);
    
    for n = 1:length(datanames)
        if ~any(strcmp(datanames{n},{'firefly','blink'}))
            trials(j).continuous.(datanames{n}) = (interp1(S2ts,trials(j).continuous.(datanames{n}),MCts,'pchip'))';
            
        elseif any(strcmp(datanames{n},{'firefly','blink'}))
            trials(j).continuous.(datanames{n}) = (interp1(S2ts,double(trials(j).continuous.(datanames{n})),MCts,'nearest'))';
            nanind = find(isnan(trials(j).continuous.(datanames{n})));
            if nanind
                trials(j).continuous.(datanames{n})(nanind) = 0;
            end
            trials(j).continuous.(datanames{n}) = logical(trials(j).continuous.(datanames{n}));
        end
        datalength(n) = length(trials(j).continuous.(datanames{n}));
        
    end
    % check SMR and MC vectors have same length
    if ~isnan(trials(j).mc.timestamp)
        MClength = length(trials(j).mc.timestamp);
        S2length = length(trials(j).continuous.ts);
        
        if S2length > MClength
            d = S2length - MClength;
            for n = 1:length(datanames)
                trials(j).continuous.(datanames{n}) = trials(j).continuous.(datanames{n})(1:end-d);
            end
            
       end
    end
end
