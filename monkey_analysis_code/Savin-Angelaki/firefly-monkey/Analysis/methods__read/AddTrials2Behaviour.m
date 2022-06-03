function trials = AddTrials2Behaviour(prs)

trials = []; % initialise
cd(prs.filepath_behv)
%% list all files to read
flist_log=dir('*.log'); 

for i=1:length(flist_log)
    % get idx letter
    if contains(flist_log(i).name,'replay')
        splt = split(flist_log(i).name,'replay_');
        filename = splt{2};
    else
        filename = flist_log(i).name;
    end
    idx_letters = regexp(filename,'[a-zA-Z]');
    % first should be m and the last 3 should be log
    idx_letters = idx_letters(2:end-3);
    delim = filename(idx_letters);
    splt = split(filename,delim);
    splt = splt{2};
    splt = split(splt,'.');
    splt = splt{1};
    fnum_log(i) = str2num(splt); 
end
[fnum_log, indx] = sort(fnum_log); 

flist_smr=dir('*.smr');
for i=1:length(flist_smr)
    % get idx letter
    idx_letters = regexp(flist_smr(i).name,'[a-zA-Z]');
    % first should be m and the last 3 should be smr
    idx_letters = idx_letters(2:end-3);
    delim = flist_smr(i).name(idx_letters);
    splt = split(flist_smr(i).name,delim);
    splt = splt{2};
    splt = split(splt,'.');
    splt = splt{1};
    fnum_smr(i) = str2num(splt);
end

if prs.visual_acc_control || prs.multisensory_acc_control
flist_mc=dir('*.txt');
for i=1:length(flist_mc)
    % get idx letter
    idx_letters = regexp(flist_mc(i).name,'[a-zA-Z]');
    % first should be m and the last 3 should be smr
    idx_letters = idx_letters(2:end-3);
    delim = flist_mc(i).name(idx_letters);
    splt = split(flist_mc(i).name,delim);
    splt = splt{2};
    splt = split(splt,'.');
    splt = splt{1};
    fnum_mc(i) = str2num(splt);
end
   
end

% flist_mat=dir('*.mat');
% for i=1:length(flist_mat), fnum_mat(i) = str2num(flist_mat(i).name(end-6:end-4)); end
nfiles = length(flist_log);

%% read files
exp = strsplit(prs.comments{:},','); exp = exp{1};

cnt_tr = 0;
for jj=1:nfiles
    i = indx(jj); 
    fprintf(['... reading ' flist_log(i).name '\n']);
    % read .log file
    if prs.visual_acc_control
        trials_log = AddLOGData_MKNOMC_EDOARDO(flist_log(i).name);
    elseif prs.multisensory_acc_control
        trials_log = AddLOGData_MKMC_EDOARDO(flist_log(i).name);
    else
        trials_log = AddLOGData(flist_log(i).name,prs);
    end
    % read all .smr files associated with this log file
    [start_ID,end_ID] = regexp(flist_log(i).name,'m\d+s\d+.log$');
    file_ID = flist_log(i).name(start_ID:end_ID-4);
    
    if i<nfiles, indx_smr = find(fnum_smr >= fnum_log(i) & fnum_smr < fnum_log(i+1));
    else indx_smr = find(fnum_smr >= fnum_log(i)); end
    trials_smr = []; t = [];
    for j = indx_smr
        smr_name = flist_smr(j).name;
        fprintf(['... reading ' smr_name  '\n']);
        data_smr = ImportSMR(smr_name);
        if prs.visual_acc_control
            [trials_smr_temp,~,t_temp] = AddSMRData_MKNOMC_EDOARDO(data_smr,[],[],prs);
            trials_smr = [trials_smr trials_smr_temp]; 
            
            if j == indx_smr(1); t = [t t_temp]; tfields = fieldnames(t);
            else; t_temp = structfun(@(x) x + t.fstop + t_temp.fstart, t_temp,'un',0); 
                for n = 1:numel(tfields)
                    if ~any(strcmp(tfields{n},{'fstart','fstop','ntrls'})); t.(tfields{n}) = [t.(tfields{n}) ; t_temp.(tfields{n})]; end
                end
            end

        elseif prs.multisensory_acc_control
            [trials_smr_temp,~,~,t_temp] = AddSMRData_MKMC_EDOARDO(data_smr,[],[],prs);
            trials_smr = [trials_smr trials_smr_temp]; 
            t = [t t_temp];
        else
            trials_smr = [trials_smr AddSMRData(data_smr,prs)];
        end
    end
    
%     found = false;
%     for fh = 1:length(flist_smr)
%         if contains(flist_smr(fh).name, file_ID)
%             smr_name = flist_smr(fh).name;
%             found = true;
%             break
%         end
%     end
%     assert(found)
%     fprintf(['... reading ' smr_name  '\n']);
%     data_smr = ImportSMR(smr_name);
%     if prs.visual_acc_control
%         [trials_smr_temp,~,t] = AddSMRData_MKNOMC_EDOARDO(data_smr,[],[],prs);
%         trials_smr = [trials_smr trials_smr_temp];        
%     elseif prs.multisensory_acc_control
%         [trials_smr_temp,~,~,t] = AddSMRData_MKMC_EDOARDO(data_smr,[],[],prs);
%         trials_smr = [trials_smr trials_smr_temp];        
%     else
%         trials_smr = [trials_smr AddSMRData(data_smr,prs)];
%     end
    
    % read all .txt files associated with this log file (MC variables)
    trials_mc = [];
    if prs.visual_acc_control || prs.multisensory_acc_control
        found = false;
        for fh = 1:length(flist_mc)
            if contains(flist_mc(fh).name, file_ID)
                mc_name = flist_mc(fh).name;
                found = true;
                break
            end
        end
        assert(found)
        fprintf(['... reading ' mc_name  '\n']);
        
        if prs.visual_acc_control
            data_mc = dlmread(mc_name);
            trials_mc = [trials_mc AddJSData_EDOARDO(data_mc,t,prs)];
        elseif prs.multisensory_acc_control
            data_mc = dlmread(mc_name);
            trials_mc = [trials_mc AddJSData_EDOARDO(data_mc,t,prs)];            
        end
    end
    
    % merge contents of .log and .smr files (and .txt MC files for acceleration control)
    ntrls_log = length(trials_log); ntrls_smr = length(trials_smr); ntrls_mc = length(trials_mc);
    cnt_tr = cnt_tr + ntrls_log;
    
    if prs.visual_acc_control || prs.multisensory_acc_control
        
        if ntrls_smr <= ntrls_log
            for j=1:length(trials_smr), trials_temp(j) = catstruct(trials_smr(j),trials_log(j),trials_mc(j)) ; end
        else  % apply a very dirty fix if spike2 was not "stopped" on time (can happen when replaying stimulus movie)
            for j=1:ntrls_log, trials_temp(j) = catstruct(trials_smr(j),trials_log(j),trials_mc(j)) ; end
            dummy_trials_log = trials_log(1:ntrls_smr-ntrls_log); dummy_trials_mc = trials_mc(1:ntrls_smr-ntrls_mc);
            for j=1:(ntrls_smr-ntrls_log); trials_temp(ntrls_log+j) = catstruct(trials_smr(ntrls_log+j),dummy_trials_log(j),dummy_trials_mc(j)); end
        end
        
    else
        
        if ntrls_smr <= ntrls_log
            %         for j=1:length(trials_smr), trials_temp(j) = catstruct(trials_smr(j),trials_log(j)) ; end
            for j=1:length(trials_smr), trials_temp(j) = concat_smr_log(trials_smr(j),trials_log(j)) ; end
        else  % apply a very dirty fix if spike2 was not "stopped" on time (can happen when replaying stimulus movie)
            %         for j=1:ntrls_log, trials_temp(j) = catstruct(trials_smr(j),trials_log(j)) ; end
            for j=1:ntrls_log, trials_temp(j) = concat_smr_log(trials_smr(j),trials_log(j)) ; end
            dummy_trials_log = trials_log(1:ntrls_smr-ntrls_log);
            %         for j=1:(ntrls_smr-ntrls_log); trials_temp(ntrls_log+j) = catstruct(trials_smr(ntrls_log+j),dummy_trials_log(j)); end
            for j=1:(ntrls_smr-ntrls_log); trials_temp(ntrls_log+j) = concat_smr_log(trials_smr(ntrls_log+j),dummy_trials_log(j)); end
        end
        
    end
    % add contents of .mat file
%     trials_temp = AddMATData(flist_mat(i).name,trials_temp);
    trials = [trials trials_temp];
    clear trials_temp;
    fprintf(['... total trials = ' num2str(length(trials)) '\n']);
end

if 0
    
    dt = prs.dt;
    maxlength = max(arrayfun(@(x) numel(x.yle),continuous));
    for k = 1:length(trials)
        
        trials(k).continuous.phi = cumsum(trials(k).continuous.w,'omitnan')*dt;

        trials(k).continuous.xfp_rel = trials(k).continuous.xfp - trials(k).continuous.xmp;
        trials(k).continuous.yfp_rel = trials(k).continuous.yfp - trials(k).continuous.ymp;
        
        R = @(phi) [cosd(phi) -sind(phi); sind(phi) cosd(phi)];
        XY = cell2mat(arrayfun(@(phi,x,y) R(phi)*[x ; y], trials(k).continuous.phi, trials(k).continuous.xfp, trials(k).continuous.yfp, 'UniformOutput', false)');
        
        trials(k).continuous.xfp_rel = XY(1,:)'; trials(k).continuous.yfp_rel = XY(2,:)';
        % r_fly_rel = sqrt(ch.xfp.^2 + ch.yfp.^2); r_fly_rel(indx_stop+1:end) = nan;
        % theta_fly_rel = atan2d(ch.xfp_rel,ch.yfp_rel); theta_fly_rel(indx_stop+1:end) = nan;
        
        x = trials(k).continuous.xfp_rel;
        y = trials(k).continuous.yfp_rel; y(y < 0) = nan;
        delta = prs.interoculardist/2;
        z = -prs.height;
        
        x = [x ; nan(maxlength-numel(x),1)];
        y = [y ; nan(maxlength-numel(y),1)];
        ye = [0.5*(trials(k).continuous.yle+trials(k).continuous.yre) ; nan(maxlength-numel(trials(k).continuous.yre),1)];
        ze = [0.5*(trials(k).continuous.zle+trials(k).continuous.zre) ; nan(maxlength-numel(trials(k).continuous.zre),1)];

        [yle_targ,zle_targ,yre_targ,zre_targ] = world2eye(x,y,z,delta);
        trials(k).continuous.ver_mean_targ = nanmean([zle_targ , zre_targ],2);
        trials(k).continuous.hor_mean_targ = nanmean([yle_targ , yre_targ],2);
        
        trials(k).continuous.tte_hor = trials(k).continuous.hor_mean_targ - ye;
        trials(k).continuous.tte_ver = trials(k).continuous.ver_mean_targ - ze;
        trials(k).continuous.tte_mag = sqrt(trials(k).continuous.tte_hor.^2 + trials(k).continuous.tte_ver.^2);
        
    end
 
continuous = [trials.continuous];
% plot individual trials
figure;
subplot(1,3,1); hold on;
arrayfun(@(x) plot(x.tte_hor),continuous); title('horizontal'); xlabel('target-tracking error [deg]');
subplot(1,3,2); hold on;
arrayfun(@(x) plot(x.tte_ver),continuous); title('vertical'); xlabel('target-tracking error [deg]');
subplot(1,3,3); hold on;
arrayfun(@(x) plot(x.tte_mag),continuous); title('magnitude'); xlabel('target-tracking error [deg]');
suptitle('Spike2')
  
% plot mean of trials
ts = (1:maxlength)*dt;
figure;
subplot(2,3,1); hold on;
shadedErrorBar(ts,nanmean(abs([continuous.hor_mean_targ]),2),nanstd(abs([continuous.hor_mean_targ]),[],2)./sqrt(numel(trials))); title('horizontal'); ylabel('target position [deg]'); xlabel('time [s]'); xlim([0 3]);
subplot(2,3,2); hold on;
shadedErrorBar(ts,nanmean([continuous.tte_ver],2),nanstd([continuous.tte_ver],[],2)./sqrt(numel(trials))); title('vertical'); ylabel('target position [deg]'); xlabel('time [s]'); xlim([0 3]);

subplot(2,3,4); hold on;
shadedErrorBar(ts,nanmean(abs([continuous.tte_hor]),2),nanstd(abs([continuous.tte_hor]),[],2)./sqrt(numel(trials))); title('horizontal'); ylabel('target-tracking error [deg]'); xlabel('time [s]'); xlim([0 3]);
subplot(2,3,5); hold on;
shadedErrorBar(ts,nanmean(abs([continuous.tte_ver]),2),nanstd(abs([continuous.tte_ver]),[],2)./sqrt(numel(trials))); title('vertical'); ylabel('target-tracking error [deg]'); xlabel('time [s]'); xlim([0 3]);
subplot(2,3,6); hold on;
shadedErrorBar(ts,nanmean(abs([continuous.tte_mag]),2),nanstd(abs([continuous.tte_mag]),[],2)./sqrt(numel(trials))); title('magnitude'); ylabel('target-tracking error [deg]'); xlabel('time [s]'); xlim([0 3]);

suptitle('Spike2')
  
    
    
end