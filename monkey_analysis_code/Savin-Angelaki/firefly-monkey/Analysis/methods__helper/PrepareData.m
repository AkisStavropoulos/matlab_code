function PrepareData(this)
%% How-to?
% in default_prs, only enable fitGAM_coupled
% pause execution of experiments.AddSessions(...,{'behv','lfps','units','pop'}) at line 120 of
% AnalysePopulation.m and then run this script until line 19 (spike raster)
% then set unitcount appropriately: best case scenario is unitcount = nunits
% then continue execution to save dataset

%% visualise spike raster
% figure; imagesc(data_concat.Yt',[0 1]);

%% remove rows that have zero spikes -- bad bad spike sorting!!!
% unitcount = [1:40 50:56];
% data_concat.Yt = Yt(:,unitcount);
% data_concat.lfp_phase = lfp_phase(:,unitcount);
% units = units(unitcount);

%% convert to struct

for k=1:numel(this.sessions(end).units), units(k) = struct(this.sessions(end).units(k)); end
if ~all(size(this.sessions(end).lfps) == [0,0])
    for k=1:numel(this.sessions(end).lfps), lfps(k) = struct(this.sessions(end).lfps(k)); end
end

trials_behv= this.sessions(end).behaviours.trials;
behv_stats = this.sessions(end).behaviours.stats;

monk_id = this.sessions(end).monk_id;
sess_id = this.sessions(end).sess_id;
exportname = ['m',num2str(monk_id),'s',num2str(sess_id),'.mat'];
prs = default_prs(monk_id,sess_id);
file_path = fullfile(prs.filepath_neur,'Pre-processing X E');
fld_exist = exist(file_path,'dir');
if fld_exist == 0
    [status, msg, msgID] = mkdir(file_path);
    if status == 0  % Never create a file/folder without catching errors
        error(msgID, msg);
    end
end

cd(file_path)
%% Export after analysis

% for k=1:numel(experiments.sessions.units), units(k) = struct(experiments.sessions.units(k)); end
% for k=1:numel(experiments.sessions.lfps), lfps(k) = struct(experiments.sessions.lfps(k)); end
% trials_behv= experiments.sessions.behaviours.trials;
% behv_stats = experiments.sessions.behaviours.stats;
% prs= default_prs(experiments.sessions.monk_id,experiments.sessions.sess_id);
% 
% exportname = ['m',num2str(experiments.sessions.monk_id),'s',num2str(experiments.sessions.sess_id),'.mat'];

%% put concatenated data in a struct
if prs.addconcat
    data_concat.Yt = Yt;
    data_concat.Xt = xt;
    data_concat.Xt(:,7) = NaN; data_concat.Xt(:,12) = NaN;
    if ~all(size(this.sessions(end).lfps) == [0,0])
        lfp_phase = [];
        for k=1:nunits
            lfp_phase(:,k) = ConcatenateTrials(var_phase{k},[],{trials_spks_temp.tspk},{continuous_temp.ts},timewindow_full);
        end
        data_concat.lfp_phase = lfp_phase;
    end
    
    cd(prs.filepath_neur);
    if ~all(size(this.sessions(end).lfps) == [0,0])
        disp(['Saving exported data: ', exportname]);
        save(exportname,'behv_stats','data_concat','lfps','prs','trials_behv','units');
    
    else
        disp(['Saving exported data: ', exportname]);
        save(exportname,'behv_stats','data_concat','prs','trials_behv','units');
    end
%     save([prs.sess_date,'.mat'],'behv_stats','data_concat','lfps','prs','trials_behv','units');
else
    %cd(prs.filepath_neur);
    if prs.compute_spectrum && ~all(size(this.sessions(end).lfps) == [0,0])

        % flag that determines if it the phase or the analytic signal is extracted 
        is_phase = false;
        % extraxt analytic

        lfp_beta = struct();
        lfp_theta = struct();
        lfp_alpha = struct();
        lfps_new = struct();
        for ch = 1: length(lfps)
            lfp_beta(ch).trials = struct();
            lfp_theta(ch).trials = struct();
            lfp_alpha(ch).trials = struct();
            lfps_new(ch).trials = struct('lfp',[]);
            lfps_new(ch).channel_id = lfps(ch).channel_id;
            lfps_new(ch).electrode_id = lfps(ch).electrode_id;
            lfps_new(ch).brain_area = lfps(ch).brain_area;

            %lfps_new(ch).stationary = lfps(ch).stationary;
            %lfps_new(ch).mobile = lfps(ch).mobile;
            %lfps_new(ch).eyesfixed = lfps.eyesfixed;
            %lfps_new(ch).eyesfree = lfps.eyesfree;
            lfps_new(ch).stats = lfps.stats;
            for tr = 1:length(lfps(ch).trials)
                lfp_beta(ch).trials(tr).lfp_beta = single(lfps(ch).trials(tr).lfp_beta);
                lfp_theta(ch).trials(tr).lfp_theta = single(lfps(ch).trials(tr).lfp_theta);
                lfp_alpha(ch).trials(tr).lfp_alpha = single(lfps(ch).trials(tr).lfp_alpha);
                lfps_new(ch).trials(tr) = rmfield(lfps(ch).trials(tr),{'lfp_beta','lfp_alpha','lfp_theta'});

            end
        end

        struct_info = whos('lfp_beta');
        if struct_info.bytes > 2^31% 2*10^9
            is_phase = true;
            for ch = 1: length(lfps)
                for tr = 1:length(lfps(ch).trials)
                    lfp_beta(ch).trials(tr).lfp_beta = angle(lfp_beta(ch).trials(tr).lfp_beta);
                    lfp_theta(ch).trials(tr).lfp_theta = angle(lfp_theta(ch).trials(tr).lfp_theta);
                    lfp_alpha(ch).trials(tr).lfp_alpha = angle(lfp_alpha(ch).trials(tr).lfp_alpha);
                end
            end
        end
        struct_info = whos('lfp_beta');
        
        if struct_info.bytes > 2^31% 2*10^9

            disp('File saved with -v7.3');
            
            lfps = lfps_new;
            disp(['Saving exported data: ', exportname]);
            save(exportname,'behv_stats','lfps','prs','trials_behv','units','-v7.3');
            exportname2 = ['lfp_beta_','m',num2str(monk_id),'s',num2str(sess_id),'.mat'];
            save(exportname2,'lfp_beta','is_phase','-v7.3');
            exportname2 = ['lfp_alpha_','m',num2str(monk_id),'s',num2str(sess_id),'.mat'];
            save(exportname2,'lfp_alpha','is_phase','-v7.3');
            exportname2 = ['lfp_theta_','m',num2str(monk_id),'s',num2str(sess_id),'.mat'];
            save(exportname2,'lfp_theta','is_phase','-v7.3');

        else
            
            lfps = lfps_new;
            disp(['Saving exported data: ', exportname]);
            save(exportname,'behv_stats','lfps','prs','trials_behv','units');
            exportname2 = ['lfp_beta_','m',num2str(monk_id),'s',num2str(sess_id),'.mat'];
            save(exportname2,'lfp_beta','is_phase');
            exportname2 = ['lfp_alpha_','m',num2str(monk_id),'s',num2str(sess_id),'.mat'];
            save(exportname2,'lfp_alpha','is_phase');
            exportname2 = ['lfp_theta_','m',num2str(monk_id),'s',num2str(sess_id),'.mat'];
            save(exportname2,'lfp_theta','is_phase');
        end
        
   elseif all(size(this.sessions(end).lfps) == [0,0])
        disp(['Saving exported data: ', exportname]);
        try
        save(exportname,'behv_stats','prs','trials_behv','units');
        catch
        save(exportname,'behv_stats','prs','trials_behv'); % in case units were not extracted 
        end
   else
        disp(['Saving exported data: ', exportname]);
        
        save(exportname,'behv_stats','lfps','prs','trials_behv','units');
    end
%     disp(['Saving exported data: ', exportname]);
%     save(exportname,'behv_stats','lfps','prs','trials_behv','units','-v7.3');
end
end
