%% 8ar analysis script

%% Extract data
clr = input('Clear everything?? ');
if clr; clear; end

choose_monkey = input('Jimmy (1) or Viktor (2)? ');
if choose_monkey==1; monkeyname='Jimmy'; 
elseif choose_monkey==2; monkeyname='Viktor'; 
end

main_dir = 'C:\Data\firefly_analysis\8ar';
saved_data_path = fullfile(main_dir,monkeyname);
cd(saved_data_path);
file_exist = dir(saved_data_path);
file_exist = file_exist(~ismember({file_exist.name},{'.','..', '.DS_Store'}));
file_exist = file_exist(~isfolder({file_exist.name}));

% choose files to load
files2load = arrayfun(@(x) ismember(x.name,{'*single_units.mat','*single_units_analysis_trialtypes.mat'}),file_exist);

% choose action
extract_data = input('Extract data? 1 for yes, 0 for no: ');
if ~extract_data && ~isempty(file_exist)
    disp('Loading Data ......... ')
    ind = find(files2load);
    for i = 1:numel(ind); load(file_exist(ind(i)).name);    disp(['......... Loaded ' file_exist(ind(i)).name]); end
else
    disp('...........Extracting data')
    
    save_file = 1;
    old_sess = 0;
    unit_filter = 'singleunit';

    su = [];
    % load saved data
    if ~isempty(file_exist) && any(files2load)
        cd(saved_data_path);
        ind = find(files2load);
        for i = 1:numel(ind); load(file_exist(ind(i)).name); disp(['......... Loaded ' file_exist(ind(i)).name]); end
        cnt = numel(su);
        old_sess = arrayfun(@(x) split(x(end).session,{'s','.'}), su,'un',0);
        old_sess = cellfun(@(x) str2double(x{2}), old_sess);
    else
        disp('No save data file found.')
        disp('...........Extracting data')
        cnt = 0;
    end
    
    % extract and add new data to struct
    if strcmp(monkeyname,'Jimmy')
    mk_dir = 'Z:\Data\Monkey2_newzdrive\Jimmy\U-probe'; cd(mk_dir);
    elseif strcmp(monkeyname,'Viktor')
    mk_dir = 'Z:\Data\Monkey2_newzdrive\Viktor\Utah Array'; cd(mk_dir);
    end
    sess_dir = dir(mk_dir);
    sess_dir = sess_dir(~contains({sess_dir.name},{'.','Sorted', '_', 'template'}));
    if strcmp(monkeyname,'Viktor')
        filter_date = datetime('Nov 30 2021','inputformat','MM dd yyyy');
        sess_date = datetime({sess_dir.name},'inputformat','MM dd yyyy');
        sessindx = sess_date > filter_date;
        sess_dir = sess_dir(sessindx);
    end
    neur_fld = 'neural data';
    data_fld = 'Pre-processing X E';
    

    for k = 1:numel(sess_dir)
        f_path = fullfile(mk_dir,sess_dir(k).name,neur_fld);
        fld = dir(f_path);
        
        % check whether session has been processed
        if any(ismember({fld.name},{data_fld}))
            file_path = fullfile(mk_dir,sess_dir(k).name,neur_fld,data_fld);
            
            temp = dir(file_path);
            if strcmp(monkeyname,'Jimmy'); mkspec = 'm73s'; elseif strcmp(monkeyname,'Viktor'); mkspec = 'm71s'; end
            indx = arrayfun(@(x) contains(x.name,mkspec) & ~contains(x.name,'lfp'),temp);
            sess_name = temp(indx).name;
            
            sessnum = strsplit(sess_name,{'s','.'}); sessnum = str2double(sessnum{2});
            
            % choose sessions
            if any(sessnum ~= old_sess)
                % load data
                disp(['............loading ' sess_name])
                A = load(fullfile(file_path,sess_name));
                cnt = cnt+1;
                
                % add unityVR indicator if missing
                if sessnum <12 && strcmp(monkeyname,'Jimmy'); A.prs.unityVR = 0;
                elseif sessnum <1000 && strcmp(monkeyname,'Viktor'); A.prs.unityVR = 0; end
                
                % add polar position from origin
                [A.behv_stats.pos_abs.r_monk, A.behv_stats.pos_abs.theta_monk] = ...
                    cellfun(@(x,y) cart2polarY(x,y), A.behv_stats.pos_abs.x_monk, A.behv_stats.pos_abs.y_monk,'un',0);
               
                % add peri-saccade window 
                fldnames = fieldnames(A.prs.ts);
                if ~any(cellfun(@(x) strcmp(x,'saccade'),fldnames))
                    A.prs.ts.saccade = [-0.5:0.02:0.5];
                end
                
                % add binrange for saccade-related tuning functions
                fldnames = fieldnames(A.prs.binrange);
                if ~any(cellfun(@(x) strcmp(x,'sac_dir'),fldnames))
                    A.prs.binrange.sac_dir = [-180 180]; % deg
                end
                if ~any(cellfun(@(x) strcmp(x,'sac_mag'),fldnames))
                    A.prs.binrange.sac_mag = [-1 35]; % deg
                end

                % Add saccade direction and magnitude
                fldnames = fieldnames(A.trials_behv(1).events);
                if ~any(cellfun(@(x) strcmp(x,'sac_mag'),fldnames))
                    presac_t = 0.1;  postsac_t = 0.15; % s
                    [A.trials_behv] = get_sac_mag_dir(A.trials_behv,presac_t,postsac_t);
                end
                
                % add target position on screen
                delta = A.prs.interoculardist/2;
                z = -A.prs.height;
                for i = 1:numel(A.trials_behv)
                    x = A.behv_stats.pos_rel.x_targ{i};
                    y = A.behv_stats.pos_rel.y_targ{i}; y(y < 0) = nan;
                    [yle_targ,zle_targ,yre_targ,zre_targ] = world2eye(x,y,z,delta);
                    ver_mean_targ{i} = nanmean([zle_targ , zre_targ],2);
                    hor_mean_targ{i} = nanmean([yle_targ , yre_targ],2);
                    
                    A.trials_behv(i).continuous.z_targ = ver_mean_targ{i};
                    A.trials_behv(i).continuous.y_targ = hor_mean_targ{i};
                end
                
                % select units
                unit_indx = arrayfun(@(x) strcmp(x.type,unit_filter),A.units);
                su(cnt).units = A.units(unit_indx);
                su(cnt).behv_stats = A.behv_stats;
                su(cnt).trials_behv = A.trials_behv;
                su(cnt).prs = A.prs;
                su(cnt).session = sess_name;
                clear A
            end
        end
    end
    
    % sort by session first before saving
    sessnum = arrayfun(@(x) split(x(end).session,{'s','.'}), su,'un',0);
    sessnum = cellfun(@(x) str2double(x{2}), sessnum);
    [~,I] = sort(sessnum);
    su = su(I);
    
    % save extracted data
    if save_file
        disp(['........saving file ' monkeyname '_single_units.mat']);
        cd(saved_data_path);
        exportname = [monkeyname,'_single_units.mat'];

        save(exportname,'su','-v7.3');
        disp('File saved!');
    end
end
%% Raster over time
LineFormat.LineWidth = 2;
figure;
for s = 1:numel(su)
    for i = 1:numel(su(s).units)
        spks = {su(s).units(i).trials.tspk}';
        spks = cellfun(@(x) x',spks,'un',0);
        % nanify empty trials
        nanindx = find(cellfun(@(x) sum(numel(x))==0, spks));
        if ~isempty(nanindx)
            for n = 1:numel(nanindx); spks{nanindx(n)} = nan; end
        end
        % t_end aligned trials
        Ntrials = numel(spks);
        t_end = arrayfun(@(x) x.events.t_end, su(s).trials_behv);
        spks_end = arrayfun(@(n) spks{n} - t_end(n), 1:Ntrials,'un',0);        
        % sort based on trial length
        trl_length = cellfun(@(x) numel(x),su(s).behv_stats.time);
        [short2long1,I] = sort(trl_length);
        short2long1 = spks(I);
        plt_var1 = short2long1;
        plt_var2 = spks_end(I);
        
        clf;
        subplot(1,2,1);
        plotSpikeRaster(plt_var1,'LineFormat',LineFormat); xlim([-0.5 3.5]); title([su(s).session ', unit ' num2str(i)])
        vline([0 0.3],{'r','k'},{'t_{on}','t_{off}'});
        xlabel('time [s]'); ylabel('trials');
        
        subplot(1,2,2); hold on;
        plotSpikeRaster(plt_var2,'LineFormat',LineFormat);xlim([-3.5 0.5]);
        vline([0 -0.3],{'r','k'},{'t_{end}','t_{rew}'});
        xlabel('time from trial end [s]'); ylabel('trials');
        pause(2);
    end
end


%% Raster over distance
LineFormat.LineWidth = 0.1;
figure;
for s = 1:numel(su)
    for i = 1:numel(su(s).units)
        spks = {su(s).units(i).trials.tspk}';
        spks = cellfun(@(x) x',spks,'un',0);
        ts = su(s).behv_stats.time';
        Ntrials = numel(ts);
        % nanify empty trials
        nanindx = find(cellfun(@(x) sum(numel(x))==0, spks));
        if ~isempty(nanindx)
            for n = 1:numel(nanindx); spks{nanindx(n)} = nan; end
        end

        % find distances of spike times
        spks2startpos = []; spks2endpos = []; timeindx = [];
        for n = 1:Ntrials
            delta_t = arrayfun(@(x) x - ts{n},spks{n},'un',0);
            [~,timeindx{n}] = cellfun(@(x) min(abs(x)), delta_t); % find time indices that correspond to spike times
            
            spks2startpos.x_monk{n,1} = su(s).behv_stats.pos_abs.x_monk{n}(timeindx{n})'; % get the position at those time indices
            spks2startpos.y_monk{n,1} = su(s).behv_stats.pos_abs.y_monk{n}(timeindx{n})';
            spks2startpos.r_monk{n,1} = su(s).behv_stats.pos_abs.r_monk{n}(timeindx{n})';
            spks2startpos.theta_monk{n,1} = su(s).behv_stats.pos_abs.theta_monk{n}(timeindx{n})';
            
            spks2endpos.x_monk{n,1} = su(s).behv_stats.pos_abs.x_monk{n}(timeindx{n})' - su(s).behv_stats.pos_abs.x_monk{n}(end);
            spks2endpos.y_monk{n,1} = su(s).behv_stats.pos_abs.y_monk{n}(timeindx{n})' - su(s).behv_stats.pos_abs.y_monk{n}(end);
            spks2endpos.r_monk{n,1} = su(s).behv_stats.pos_abs.r_monk{n}(timeindx{n})' - su(s).behv_stats.pos_abs.r_monk{n}(end);
            spks2endpos.theta_monk{n,1} = su(s).behv_stats.pos_abs.theta_monk{n}(timeindx{n})' - su(s).behv_stats.pos_abs.theta_monk{n}(end);
        end
        
        % sort based on distance
        final_y = cellfun(@(x) x(end),su(s).behv_stats.pos_abs.y_monk);
        final_x = cellfun(@(x) x(end),su(s).behv_stats.pos_abs.x_monk);
        final_r = cellfun(@(x) x(end),su(s).behv_stats.pos_abs.r_monk);
        final_th = cellfun(@(x) x(end),su(s).behv_stats.pos_abs.theta_monk);
        
        [short2long_y,I] = sort(final_y);        y_start = spks2startpos.y_monk(I);     y_end = spks2endpos.y_monk(I);
        [short2long_x,I] = sort(final_x);        x_start = spks2startpos.x_monk(I);     x_end = spks2endpos.x_monk(I);
        [short2long_r,I] = sort(final_r);        r_start = spks2startpos.r_monk(I);     r_end = spks2endpos.r_monk(I);
        [short2long_th,I] = sort(final_th);        th_start = spks2startpos.theta_monk(I);     th_end = spks2endpos.theta_monk(I);
        hist_y_start = [y_start{:}];
        hist_x_start = [x_start{:}]; 
        hist_r_start = [r_start{:}];
        hist_th_start = [th_start{:}]; 
                
        clf;
        subplot(2,4,1);
        plotSpikeRaster(y_start,'SpikeDuration',1,'LineFormat',LineFormat);
        xlim([-50 500]); title(['session ' num2str(s) ', unit ' num2str(i)])
        vline(0,'r','start pos'); 
        xlabel('y distance from origin [cm]'); ylabel('trials');

        subplot(2,4,2);
        plotSpikeRaster(x_start,'SpikeDuration',1,'LineFormat',LineFormat);
        xlim([-200 200]); title(['session ' num2str(s) ', unit ' num2str(i)])
        vline(0,'r','start pos');
        xlabel('x distance from origin [cm]'); ylabel('trials');
        
        subplot(2,4,3);
        plotSpikeRaster(y_end,'SpikeDuration',1,'LineFormat',LineFormat);
        xlim([-500 50]); title(['session ' num2str(s) ', unit ' num2str(i)])
        vline(0,'r','end pos');
        xlabel('y distance from origin [cm]'); ylabel('trials');
        
        subplot(2,4,4);
        plotSpikeRaster(x_end,'SpikeDuration',1,'LineFormat',LineFormat);
        xlim([-200 200]); title(['session ' num2str(s) ', unit ' num2str(i)])
        vline(0,'r','end pos'); 
        xlabel('x distance from origin [cm]'); ylabel('trials');


        subplot(2,4,5);
        plotSpikeRaster(r_start,'SpikeDuration',1,'LineFormat',LineFormat);
        xlim([-50 500]); title(['session ' num2str(s) ', unit ' num2str(i)])
        vline(0,'r','start pos'); 
        xlabel('r distance from origin [cm]'); ylabel('trials');

        subplot(2,4,6);
        plotSpikeRaster(th_start,'SpikeDuration',1,'LineFormat',LineFormat);
        xlim([-40 40]); title(['session ' num2str(s) ', unit ' num2str(i)])
        vline(0,'r','start pos'); 
        xlabel('\theta distance from origin [cm]'); ylabel('trials');

        subplot(2,4,7);
        plotSpikeRaster(r_end,'SpikeDuration',1,'LineFormat',LineFormat);
        xlim([-500 50]); title(['session ' num2str(s) ', unit ' num2str(i)])
        vline(0,'r','end pos'); 
        xlabel('r distance from origin [cm]'); ylabel('trials');

        subplot(2,4,8);
        plotSpikeRaster(th_end,'SpikeDuration',1,'LineFormat',LineFormat);
        xlim([-40 40]); title(['session ' num2str(s) ', unit ' num2str(i)])
        vline(0,'r','end pos'); 
        xlabel('\theta distance from origin [cm]'); ylabel('trials');
        
        pause(2)
    end
end

%% Raster aligned to saccade onset
LineFormat.LineWidth = 2;
figure;
for s = 1:numel(su)
    prs = su(s).prs;
    for i = 1:numel(su(s).units)
        trials_spks = su(s).units(i).trials';
        % saccade-aligned trials
        Ntrials = numel(trials_spks);
        % select
        t_sac = arrayfun(@(x) x.events.t_sac, su(s).trials_behv,'un',0);
        spks_saccade = ShiftSpikes4saccades(trials_spks,t_sac);
        cnt = 0;
        tspk = [];
        for k=1:Ntrials
            tspk_temp = spks_saccade(k).tspk;
            for n = 1:numel(t_sac{k})
                if ~isempty(tspk_temp)
                    cnt = cnt+1;
                    tspk{cnt}(1,:) = tspk_temp(:,n);
                end
            end
        end
        plt_var1 = tspk;
%         % sort based on trial length
%         trl_length = cellfun(@(x) numel(x),su(s).behv_stats.time);
%         [short2long1,I] = sort(trl_length);
%         short2long1 = spks(I);
%         plt_var1 = short2long1;
%         plt_var2 = spks_end(I);
        
        
        clf;
        plotSpikeRaster(plt_var1,'LineFormat',LineFormat); xlim([-0.5 3.5]); 
        xlim([min(prs.ts.saccade)*2 max(prs.ts.saccade)*2]); vline(0,'k'); xlabel('Time from saccade (s)'); ylabel('Saccade number');
        title([su(s).session ', area: ' su(s).prs.area{:} ', unit ' num2str(i)]);
        
        pause(2);
    end
end

%% Analyze units using Kaushik's function
for s = 1:numel(su)
    disp(' '); disp(['.....Analyzing session ' su(s).session '(' num2str(s) '/' num2str(numel(su)) ')']);
    for i = 1:numel(su(s).units)
        disp(['Unit ' num2str(i)]); 
        
        trials_spks = su(s).units(i).trials;
        trials_behv = su(s).trials_behv;
        behv_stats = su(s).behv_stats;
        lfps = [];
        prs = su(s).prs;
        prs.ts.saccade = -0.5:0.02:0.5;
        
        stats = AnalyseUnit(trials_spks,trials_behv,behv_stats,lfps,prs);
        
        su_analysis(s).units(i).stats = stats;
        su_analysis(s).session = su(s).session;
        
    end
end

% save
disp('Saving analyzed units.......');
cd('C:\Data\firefly_analysis\8ar');
save('single_units_analysis_trialtypes.mat','su_analysis','-v7.3');
disp('........Saved!')


%% Plot tuning curves
plot_type = {'rate_move','rate_targ','rate_stop','rate_rew','rate_sac'};
Nrow = 2;
Ncol = ceil(numel(plot_type)/Nrow);

for s = 1:numel(su)
   for i = 1:numel(su(s).units)
       figure;
       
       prs = su(s).prs;
       behv.stats = su(s).behv_stats;
       behv.trials = su(s).trials_behv;
       unit.trials = su(s).units.trials;
       unit.stats = su_analysis(s).units(i).stats;
       trial_type = 'all';
       
       for n = 1:numel(plot_type)
           subplot(Nrow,Ncol,n); hold on;
           PlotUnit_Akis(prs,behv,unit,plot_type{n},trial_type);
       end
       sgtitle([su(s).session ', \bfarea:\rm ' su(s).prs.area{:} ', \bfunit:\rm ' num2str(i)]);
       
       figure;
       PlotUnit_Akis(prs,behv,unit,'tuning_continuous',trial_type);
       sgtitle([su(s).session ', \bfarea:\rm ' su(s).prs.area{:} ', \bfunit:\rm ' num2str(i)]);

   end
end

%% Event-aligned firing rates, individual units ********************************
varnames = fieldnames(su_analysis(1).units(1).stats.trialtype.all.events);
ncol = ceil(numel(varnames)/2);
for s = 1:numel(su)
    figure;

    subplot(2,ncol,1); hold on;
    rate_move = cell2mat(arrayfun(@(x) x.stats.trialtype.all.events.move.rate(:), su_analysis(s).units,'un',0));
    ts = su_analysis(s).units(1).stats.trialtype.all.events.move.time;
    plot(ts,rate_move); vline(0,'k'); xlabel('time from movement onset [s]'); ylabel('firing rate [Hz]'); title('movement onset');

    subplot(2,ncol,2); hold on;
    rate_targ = cell2mat(arrayfun(@(x) x.stats.trialtype.all.events.target.rate(:), su_analysis(s).units,'un',0));
    ts = su_analysis(s).units(1).stats.trialtype.all.events.target.time;
    plot(ts,rate_targ); vline(0,'k'); xlabel('time from target onset [s]'); ylabel('firing rate [Hz]'); title('target onset');

    subplot(2,ncol,3); hold on;
    rate_stop = cell2mat(arrayfun(@(x) x.stats.trialtype.all.events.stop.rate(:), su_analysis(s).units,'un',0));
    ts = su_analysis(s).units(1).stats.trialtype.all.events.stop.time;
    plot(ts,rate_stop); vline(0,'k'); xlabel('time from stop [s]'); ylabel('firing rate [Hz]'); title('stop');

    subplot(2,ncol,4); hold on;
    rate_rew = cell2mat(arrayfun(@(x) x.stats.trialtype.all.events.reward.rate(:), su_analysis(s).units,'un',0));
    ts = su_analysis(s).units(1).stats.trialtype.all.events.reward.time;
    plot(ts,rate_rew); vline(0,'k'); xlabel('time from reward [s]'); ylabel('firing rate [Hz]'); title('reward');

    subplot(2,ncol,5); hold on;
    rate_sac = cell2mat(arrayfun(@(x) x.stats.trialtype.all.events.saccade.rate(:), su_analysis(s).units,'un',0));
    ts = su_analysis(s).units(1).stats.trialtype.all.events.saccade.time;
    plot(ts,rate_sac); vline(0,'k'); xlabel('time from saccade [s]'); ylabel('firing rate [Hz]'); title('saccade');

    sgtitle([su(s).session ', \bfarea:\rm ' su(s).prs.area{:} ', \bfall units\rm (' num2str(numel(su(s).units)) ')']);

end

%% Event-aling firing rates, group units based on setup (8ar or 7a) *******************************
unityVR = arrayfun(@(x) x.prs.unityVR, su);
recording_area = '8ar';
brain_area = arrayfun(@(x) strcmp(recording_area,x.prs.area), su);

selectindx = [{unityVR & brain_area} {~unityVR & brain_area}];
set_up = {'unity','spike2'};
colr = brewermap(3,'Set1');
cond = 'all';
indivunits = 0;

figure;

rate_move = []; rate_targ = []; rate_stop = []; rate_rew = []; rate_sac = [];
for k = 1:numel(selectindx)
    
    nconds = numel(su_analysis(1).units(1).stats.trialtype.(cond));
    condcheck = [];
    for s = 1:numel(su); condcheck(s) = numel(su_analysis(s).units(1).stats.trialtype.(cond)); end
    if ~strcmp('all',cond); condcheck = condcheck==2; else; condcheck = condcheck==1; end

    analysis_temp = su_analysis(selectindx{k} & condcheck);
    units = [analysis_temp.units];

    for m = 1:nconds
    
    rate_move.(set_up{k}).rate = cell2mat(arrayfun(@(x) x.stats.trialtype.(cond)(m).events.move.rate(:), units,'un',0));
    rate_move.(set_up{k}).ts = cell2mat(arrayfun(@(x) x.stats.trialtype.(cond)(m).events.move.time(:), units,'un',0));
    
    rate_targ.(set_up{k}).rate = cell2mat(arrayfun(@(x) x.stats.trialtype.(cond)(m).events.target.rate(:), units,'un',0));
    rate_targ.(set_up{k}).ts = cell2mat(arrayfun(@(x) x.stats.trialtype.(cond)(m).events.target.time(:), units,'un',0));
    
    rate_stop.(set_up{k}).rate = cell2mat(arrayfun(@(x) x.stats.trialtype.(cond)(m).events.stop.rate(:), units,'un',0));
    rate_stop.(set_up{k}).ts = cell2mat(arrayfun(@(x) x.stats.trialtype.(cond)(m).events.stop.time(:), units,'un',0));
    
    rate_rew.(set_up{k}).rate = cell2mat(arrayfun(@(x) x.stats.trialtype.(cond)(m).events.reward.rate(:), units,'un',0));
    rate_rew.(set_up{k}).ts = cell2mat(arrayfun(@(x) x.stats.trialtype.(cond)(m).events.reward.time(:), units,'un',0));
    
    rate_sac.(set_up{k}).rate = cell2mat(arrayfun(@(x) x.stats.trialtype.(cond)(m).events.saccade.rate(:), units,'un',0));
    rate_sac.(set_up{k}).ts = cell2mat(arrayfun(@(x) x.stats.trialtype.(cond)(m).events.saccade.time(:), units,'un',0));
    
    
    % plot
    if indivunits % individual units
        subplot(2,ncol,1); hold on;
        plot(rate_move.(set_up{k}).ts,rate_move.(set_up{k}).rate,'color',colr(k,:).*(1/m));
        vline(0,'k'); xlabel('time from movement onset [s]'); ylabel('firing rate [Hz]'); title('movement onset');
        
        subplot(2,ncol,2); hold on;
        plot(rate_targ.(set_up{k}).ts,rate_targ.(set_up{k}).rate,'color',colr(k,:).*(1/m));
        vline(0,'k'); xlabel('time from target onset [s]'); ylabel('firing rate [Hz]'); title('target onset');
        
        subplot(2,ncol,3); hold on;
        plot(rate_stop.(set_up{k}).ts,rate_stop.(set_up{k}).rate,'color',colr(k,:).*(1/m));
        vline(0,'k'); xlabel('time from stop [s]'); ylabel('firing rate [Hz]'); title('stop');
        
        subplot(2,ncol,4); hold on;
        plot(rate_rew.(set_up{k}).ts,rate_rew.(set_up{k}).rate,'color',colr(k,:).*(1/m));
        vline(0,'k'); xlabel('time from reward [s]'); ylabel('firing rate [Hz]'); title('reward');
        
        subplot(2,ncol,5); hold on;
        plot(rate_sac.(set_up{k}).ts,rate_sac.(set_up{k}).rate,'color',colr(k,:).*(1/m));
        vline(0,'k'); xlabel('time from saccade [s]'); ylabel('firing rate [Hz]'); title('saccade');
        
    else % average
        
        subplot(2,ncol,1); hold on;
        shadedErrorBar(nanmean(rate_move.(set_up{k}).ts,2),nanmean(rate_move.(set_up{k}).rate,2),nanstd(rate_move.(set_up{k}).rate,[],2)./sqrt(numel(units)),'lineprops',{'color',colr(k,:).*(1/m)});
        vline(0,'k'); xlabel('time from movement onset [s]'); ylabel('firing rate [Hz]'); title('movement onset');
        
        subplot(2,ncol,2); hold on;
        shadedErrorBar(nanmean(rate_targ.(set_up{k}).ts,2),nanmean(rate_targ.(set_up{k}).rate,2),nanstd(rate_targ.(set_up{k}).rate,[],2)./sqrt(numel(units)),'lineprops',{'color',colr(k,:).*(1/m)});
        vline(0,'k'); xlabel('time from target onset [s]'); ylabel('firing rate [Hz]'); title('target onset');
        
        subplot(2,ncol,3); hold on;
        shadedErrorBar(nanmean(rate_stop.(set_up{k}).ts,2),nanmean(rate_stop.(set_up{k}).rate,2),nanstd(rate_stop.(set_up{k}).rate,[],2)./sqrt(numel(units)),'lineprops',{'color',colr(k,:).*(1/m)});
        vline(0,'k'); xlabel('time from stop [s]'); ylabel('firing rate [Hz]'); title('stop');
        
        subplot(2,ncol,4); hold on;
        shadedErrorBar(nanmean(rate_rew.(set_up{k}).ts,2),nanmean(rate_rew.(set_up{k}).rate,2),nanstd(rate_rew.(set_up{k}).rate,[],2)./sqrt(numel(units)),'lineprops',{'color',colr(k,:).*(1/m)});
        vline(0,'k'); xlabel('time from reward [s]'); ylabel('firing rate [Hz]'); title('reward');
        
        subplot(2,ncol,5); hold on;
        shadedErrorBar(nanmean(rate_sac.(set_up{k}).ts,2),nanmean(rate_sac.(set_up{k}).rate,2),nanstd(rate_sac.(set_up{k}).rate,[],2)./sqrt(numel(units)),'lineprops',{'color',colr(k,:).*(1/m)});
        vline(0,'k'); xlabel('time from saccade [s]'); ylabel('firing rate [Hz]'); title('saccade');
    end 
    end
end
sgtitle(['\bfarea:\rm ' recording_area ' -- \bf\color{red}Unity \color{black}vs \color{blue}Spike2\rm']);


%% Tuning functions, individual units
nrows = 3;
varnames = fieldnames(su_analysis(1).units(1).stats.trialtype.all.continuous);
varnames = varnames(~ismember(varnames,{'a','alpha','vheye','vhtarg','vhtte','rtheta','sac_magdir'}));
ncol = ceil(numel(varnames)/nrows);

figure;
for s = 1:numel(su)
    clf;
    subplot(nrows,ncol,1); hold on;
    rate_heye.mu = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.heye.tuning.rate.mu(:), su_analysis(s).units,'un',0));
    rate_heye.se = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.heye.tuning.rate.sem(:), su_analysis(s).units,'un',0));
    stim_heye = su_analysis(s).units(1).stats.trialtype.all.continuous.heye.tuning.stim.mu;
    plot(stim_heye,rate_heye.mu); vline(0,'k'); xlabel('hor. eye position [deg]'); ylabel('firing rate [Hz]'); title('horizontal eye position');

    subplot(nrows,ncol,2); hold on;
    rate_veye.mu = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.veye.tuning.rate.mu(:), su_analysis(s).units,'un',0));
    rate_veye.se = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.veye.tuning.rate.sem(:), su_analysis(s).units,'un',0));
    stim_veye = su_analysis(s).units(1).stats.trialtype.all.continuous.veye.tuning.stim.mu;
    plot(stim_veye,rate_veye.mu); vline(0,'k'); xlabel('ver. eye position [deg]'); ylabel('firing rate [Hz]'); title('vertical eye position');
    
    subplot(nrows,ncol,3); hold on;
    rate_htarg.mu = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.htarg.tuning.rate.mu(:), su_analysis(s).units,'un',0));
    rate_htarg.se = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.htarg.tuning.rate.sem(:), su_analysis(s).units,'un',0));
    stim_htarg = su_analysis(s).units(1).stats.trialtype.all.continuous.htarg.tuning.stim.mu;
    plot(stim_htarg,rate_htarg.mu); vline(0,'k'); xlabel('hor. target position [deg]'); ylabel('firing rate [Hz]'); title('hor. target position');

    subplot(nrows,ncol,4); hold on;
    rate_vtarg.mu = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.vtarg.tuning.rate.mu(:), su_analysis(s).units,'un',0));
    rate_vtarg.se = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.vtarg.tuning.rate.sem(:), su_analysis(s).units,'un',0));
    stim_vtarg = su_analysis(s).units(1).stats.trialtype.all.continuous.vtarg.tuning.stim.mu;
    plot(stim_vtarg,rate_vtarg.mu); vline(0,'k'); xlabel('ver. target position [deg]'); ylabel('firing rate [Hz]'); title('ver. target position');

    subplot(nrows,ncol,5); hold on;
    rate_htte.mu = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.htte.tuning.rate.mu(:), su_analysis(s).units,'un',0));
    rate_htte.se = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.htte.tuning.rate.sem(:), su_analysis(s).units,'un',0));
    stim_htte = su_analysis(s).units(1).stats.trialtype.all.continuous.htte.tuning.stim.mu;
    plot(stim_htte,rate_htte.mu); vline(0,'k'); xlabel('hor. tracking error [deg]'); ylabel('firing rate [Hz]'); title('hor. tracking error');

    subplot(nrows,ncol,6); hold on;
    rate_vtte.mu = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.vtte.tuning.rate.mu(:), su_analysis(s).units,'un',0));
    rate_vtte.se = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.vtte.tuning.rate.sem(:), su_analysis(s).units,'un',0));
    stim_vtte = su_analysis(s).units(1).stats.trialtype.all.continuous.vtte.tuning.stim.mu;
    plot(stim_vtte,rate_vtte.mu); vline(0,'k'); xlabel('ver. tracking error [deg]'); ylabel('firing rate [Hz]'); title('ver. tracking error');

    subplot(nrows,ncol,7); hold on;
    rate_ttemag.mu = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.ttemag.tuning.rate.mu(:), su_analysis(s).units,'un',0));
    rate_ttemag.se = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.ttemag.tuning.rate.sem(:), su_analysis(s).units,'un',0));
    stim_ttemag = su_analysis(s).units(1).stats.trialtype.all.continuous.ttemag.tuning.stim.mu;
    plot(stim_ttemag,rate_ttemag.mu); vline(0,'k'); xlabel('tracking error mag. [deg]'); ylabel('firing rate [Hz]'); title('tracking error magnitude');

    subplot(nrows,ncol,8); hold on;
    rate_sac_dir.mu = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.sac_dir.tuning.rate.mu(:), su_analysis(s).units,'un',0));
    rate_sac_dir.se = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.sac_dir.tuning.rate.sem(:), su_analysis(s).units,'un',0));
    stim_sac_dir = su_analysis(s).units(1).stats.trialtype.all.continuous.sac_dir.tuning.stim.mu;
    plot(stim_sac_dir,rate_sac_dir.mu); vline(0,'k'); xlabel('saccade direction [deg]'); ylabel('firing rate [Hz]'); title('saccade direction');

    subplot(nrows,ncol,9); hold on;
    rate_sac_mag.mu = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.sac_mag.tuning.rate.mu(:), su_analysis(s).units,'un',0));
    rate_sac_mag.se = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.sac_mag.tuning.rate.sem(:), su_analysis(s).units,'un',0));
    stim_sac_mag = su_analysis(s).units(1).stats.trialtype.all.continuous.sac_mag.tuning.stim.mu;
    plot(stim_sac_mag,rate_sac_mag.mu); vline(0,'k'); xlabel('saccade magnitude [deg]'); ylabel('firing rate [Hz]'); title('saccade magnitude');

    subplot(nrows,ncol,10); hold on;
    rate_r.mu = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.r.tuning.rate.mu(:), su_analysis(s).units,'un',0));
    rate_r.se = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.r.tuning.rate.sem(:), su_analysis(s).units,'un',0));
    stim_r = su_analysis(s).units(1).stats.trialtype.all.continuous.r.tuning.stim.mu;
    plot(stim_r,rate_r.mu); vline(0,'k'); xlabel('radial distance [cm]'); ylabel('firing rate [Hz]'); title('radial distance');

    subplot(nrows,ncol,11); hold on;
    rate_theta.mu = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.theta.tuning.rate.mu(:), su_analysis(s).units,'un',0));
    rate_theta.se = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.theta.tuning.rate.sem(:), su_analysis(s).units,'un',0));
    stim_theta = su_analysis(s).units(1).stats.trialtype.all.continuous.theta.tuning.stim.mu;
    plot(stim_theta,rate_theta.mu); vline(0,'k'); xlabel('angular distance [deg]'); ylabel('firing rate [Hz]'); title('angular distance');

    subplot(nrows,ncol,12); hold on;
    rate_d.mu = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.d.tuning.rate.mu(:), su_analysis(s).units,'un',0));
    rate_d.se = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.d.tuning.rate.sem(:), su_analysis(s).units,'un',0));
    stim_d = su_analysis(s).units(1).stats.trialtype.all.continuous.d.tuning.stim.mu;
    plot(stim_d,rate_d.mu); vline(0,'k'); xlabel('distance traveled [cm]'); ylabel('firing rate [Hz]'); title('distance traveled');
    
    subplot(nrows,ncol,13); hold on;
    rate_phi.mu = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.phi.tuning.rate.mu(:), su_analysis(s).units,'un',0));
    rate_phi.se = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.phi.tuning.rate.sem(:), su_analysis(s).units,'un',0));
    stim_phi = su_analysis(s).units(1).stats.trialtype.all.continuous.phi.tuning.stim.mu;
    plot(stim_phi,rate_phi.mu); vline(0,'k'); xlabel('heading [deg]'); ylabel('firing rate [Hz]'); title('heading');

    subplot(nrows,ncol,14); hold on;
    rate_v.mu = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.v.tuning.rate.mu(:), su_analysis(s).units,'un',0));
    rate_v.se = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.v.tuning.rate.sem(:), su_analysis(s).units,'un',0));
    stim_v = su_analysis(s).units(1).stats.trialtype.all.continuous.v.tuning.stim.mu;
    plot(stim_v,rate_v.mu); vline(0,'k'); xlabel('lin. velocity [cm/s]'); ylabel('firing rate [Hz]'); title('lin. velocity');
    
    subplot(nrows,ncol,15); hold on;
    rate_w.mu = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.w.tuning.rate.mu(:), su_analysis(s).units,'un',0));
    rate_w.se = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.w.tuning.rate.sem(:), su_analysis(s).units,'un',0));
    stim_w = su_analysis(s).units(1).stats.trialtype.all.continuous.w.tuning.stim.mu;
    plot(stim_w,rate_w.mu); vline(0,'k'); xlabel('ang. velocity [deg/s]'); ylabel('firing rate [Hz]'); title('ang. velocity');

    sgtitle([su(s).session ', \bfarea:\rm ' su(s).prs.area{:} ', \bfall units\rm (' num2str(numel(su(s).units)) ')']);

    pause(4);
end


%% Tuning functions, group units based on setup (8ar or 7a)
unityVR = arrayfun(@(x) x.prs.unityVR, su);
recording_area = '8ar';
brain_area = arrayfun(@(x) strcmp(recording_area,x.prs.area), su);

selectindx = [{unityVR & brain_area} {~unityVR & brain_area}];
set_up = {'unity','spike2'};
colr = brewermap(3,'Set1');
cond = 'all';
indivunits = 0;

figure;
varnames = fieldnames(su_analysis(1).units(1).stats.trialtype.all.continuous);
varnames = varnames(~ismember(varnames,{'a','alpha','vheye','vhtarg','vhtte','rtheta','sac_magdir'}));
nrows = 3;
ncol = ceil(numel(varnames)/nrows);

vars = [];
for k = 1:numel(selectindx)
    
    nconds = numel(su_analysis(1).units(1).stats.trialtype.(cond));
    condcheck = [];
    for s = 1:numel(su); condcheck(s) = numel(su_analysis(s).units(1).stats.trialtype.(cond)); end
    if ~strcmp('all',cond); condcheck = condcheck==2; else; condcheck = condcheck==1; end

    analysis_temp = su_analysis(selectindx{k} & condcheck);
    units = [analysis_temp.units];

    for m = 1:nconds
    
        for n = 1:length(varnames)
            vars.(varnames{n}).(set_up{k}).rate = cell2mat(arrayfun(@(x) x.stats.trialtype.(cond)(m).continuous.(varnames{n}).tuning.rate.mu(:), units,'un',0));
            vars.(varnames{n}).(set_up{k}).stim = cell2mat(arrayfun(@(x) x.stats.trialtype.(cond)(m).continuous.(varnames{n}).tuning.stim.mu(:),units,'un',0));
        end
        
    % plot
    if indivunits % individual units
        
        for n = 1:length(varnames)
            subplot(nrows,ncol,n); hold on;
            plot(vars.(varnames{n}).(set_up{k}).stim,vars.(varnames{n}).(set_up{k}).rate,'color',colr(k,:).*(1/m));
            vline(0,'k'); xlabel(varnames{n}); ylabel('firing rate [Hz]'); title(varnames{n});
        end
        
    else % average
        
        for n = 1:length(varnames)
        subplot(nrows,ncol,n); hold on;
        shadedErrorBar(nanmedian(vars.(varnames{n}).(set_up{k}).stim,2),nanmean(vars.(varnames{n}).(set_up{k}).rate,2),nanstd(vars.(varnames{n}).(set_up{k}).rate,[],2)./sqrt(numel(units)),'lineprops',{'color',colr(k,:).*(1/m)});
        vline(0,'k'); xlabel(varnames{n}); ylabel('firing rate [Hz]'); title(varnames{n});
        end
    end 
    end
end
sgtitle(['\bfarea:\rm ' recording_area ' -- \bf\color{red}Unity \color{black}vs \color{blue}Spike2\rm']);

%% Compare fraction tuned across setups (Unity vs Spike2)
unityVR = arrayfun(@(x) x.prs.unityVR, su)';
brain_areas = unique(arrayfun(@(x) x.prs.area, su));
setups = {'unity','spike2'};

rec_area = cell2mat(arrayfun(@(x) strcmp(brain_areas,x.prs.area)', su,'un',0))';

selectindx = [{unityVR & rec_area(:,1)} {~unityVR & rec_area(:,1)} {unityVR & rec_area(:,2)} {~unityVR & rec_area(:,2)}];
vals = {'unity+7a','spike2+7a','unity+8ar','spike2+8ar'};

fr_tuned = [];
for k = 1:numel(selectindx)  

    analysis_temp = su_analysis(selectindx{k});
    units = [analysis_temp.units];

    for n = 1:length(varnames)
        fr_tuned(k).val = vals{k};
        fr_tuned(k).vars.(varnames{n}).pval = cell2mat(arrayfun(@(x) x.stats.trialtype.all.continuous.(varnames{n}).tuning.pval, units,'un',0));
        fr_tuned(k).vars.(varnames{n}).tuned =  fr_tuned(k).vars.(varnames{n}).pval < 0.05;
    end

end

% plot all
nrows = 3;
ncol = ceil(numel(varnames)/nrows);
colr = brewermap(2,'Set1'); colr = repmat(colr,2,1); colr = diag([0.5 0.5 1 1])*colr;

figure;
fr = [];
for n = 1:length(varnames)
    subplot(nrows,ncol,n); hold on;
    for k = 1:numel(fr_tuned)
        fr.(varnames{n})(k) = sum(fr_tuned(k).vars.(varnames{n}).tuned)/numel(fr_tuned(k).vars.(varnames{n}).tuned);
        
        bar(k,fr.(varnames{n})(k),'facecolor',colr(k,:));
    end
    xticks(1:numel(fr_tuned)); xticklabels({fr_tuned.val}); xtickangle(90);
    ylabel('% tuned'); title(varnames{n});
end
sgtitle('Fraction tuned comparison across setups');

% compare fractions between setups (X^2 test)
chi2stat = []; p = [];
for n = 1:numel(varnames)
    cnt=1;
    for k = 1:2:numel(fr_tuned)
        % Observed data
        n1 = sum(fr_tuned(k).vars.(varnames{n}).tuned); N1 = numel(fr_tuned(k).vars.(varnames{n}).tuned);
        n2 = sum(fr_tuned(k+1).vars.(varnames{n}).tuned); N2 = numel(fr_tuned(k+1).vars.(varnames{n}).tuned);
        % Pooled estimate of proportion
        p0 = (n1+n2) / (N1+N2);
        % Expected counts under H0 (null hypothesis)
        n10 = N1 * p0;
        n20 = N2 * p0;
        % Chi-square test, by hand
        observed = [n1 N1-n1 n2 N2-n2];
        expected = [n10 N1-n10 n20 N2-n20];
        chi2stat.(varnames{n})(cnt) = sum((observed-expected).^2 ./ expected);
        p.(varnames{n})(cnt) = 1 - chi2cdf(chi2stat.(varnames{n})(cnt),1);
        
        cnt=cnt+1;
    end
end
chi2significant = structfun(@(x) x < 0.05, p,'un',0)

%% Separate units whose firing increases/decreases with TTE magnitude
brain_area = '8ar';
sessindx = arrayfun(@(x) strcmp(x.prs.area,brain_area),su); sessindx(end) = []; % temporary fix for s27 without eye tracking

str_temp = su_analysis(sessindx);

units_area = [str_temp.units];

ttemag.stim = arrayfun(@(x) x.stats.trialtype.all.continuous.ttemag.tuning.stim.mu,units_area,'un',0);
ttemag.stim = nanmedian(cell2mat(ttemag.stim'));
ttemag.rate = arrayfun(@(x) x.stats.trialtype.all.continuous.ttemag.tuning.rate.mu,units_area,'un',0);
ttemag.rate_filt = arrayfun(@(x) movmean(x.stats.trialtype.all.continuous.ttemag.tuning.rate.mu,5,'omitnan'),units_area,'un',0);

incrindx = cellfun(@(x) x(end)-x(1) > 0,ttemag.rate_filt);
decrindx = cellfun(@(x) x(end)-x(1) < 0,ttemag.rate_filt);

figure; subplot(1,3,1); hold on;
plot(ttemag.stim,cell2mat(ttemag.rate(incrindx)'),'b'); 
xlabel('tracking-error magnitude [deg]'); ylabel('firing rate [Hz]'); title('increasing activity');
subplot(1,3,2); hold on;
plot(ttemag.stim,cell2mat(ttemag.rate(decrindx)'),'r');
xlabel('tracking-error magnitude [deg]'); ylabel('firing rate [Hz]'); title('decreasing activity');
subplot(1,3,3); hold on;
shadedErrorBar(ttemag.stim,nanmean(cell2mat(ttemag.rate(incrindx)')),nanstd(cell2mat(ttemag.rate(incrindx)'))./sqrt(numel(units_area)),'lineprops',{'color','b'}); 
shadedErrorBar(ttemag.stim,nanmean(cell2mat(ttemag.rate(decrindx)')),nanstd(cell2mat(ttemag.rate(decrindx)'))./sqrt(numel(units_area)),'lineprops',{'color','r'});
xlabel('tracking-error magnitude [deg]'); ylabel('firing rate [Hz]'); title('average');

sgtitle(brain_area)
%% Compare peri-saccadic activity between saccades that decrease/increase target-tracking error
disp('Comparison of peri-saccadic activity between decrease/increase of target-tracking error');

presac = floor(0.3/prs.dt);
postsac = floor(0.3/prs.dt);

for s = 1:numel(su)
    disp(' '); disp(['.....Session ' su(s).session '(' num2str(s) '/' num2str(numel(su)) ')']);
%     figure;
    for i = 1:numel(su(s).units)
        disp(['Unit ' num2str(i)]); 
        
        trials_spks = su(s).units(i).trials;
        trials_behv = su(s).trials_behv;
        behv_stats = su(s).behv_stats;
        lfps = [];
        prs = su(s).prs;
        prs.ts.saccade = -0.5:0.02:0.5;
        
        % compute tracking error
        ts = arrayfun(@(x) x.continuous.ts, trials_behv,'un',0);
        
        veye = arrayfun(@(x) nanmean([x.continuous.zle(:)  x.continuous.zre(:)],2), trials_behv,'un',0);
        heye = arrayfun(@(x) nanmean([x.continuous.yle(:)  x.continuous.yre(:)],2), trials_behv,'un',0);
        vtarg = arrayfun(@(x) x.continuous.z_targ, trials_behv,'un',0);
        htarg = arrayfun(@(x) x.continuous.y_targ, trials_behv,'un',0);
        
        vtte = cellfun(@(x,y) x - y, vtarg,veye,'un',0);
        htte = cellfun(@(x,y) x - y, htarg,heye,'un',0);
        ttemag = cellfun(@(x1,x2,y1,y2) sqrt((x1-x2).^2 + (y1-y2).^2), vtarg,veye,htarg,heye,'un',0);
        
        % compare TTE before and after saccades
        t_targOFF = arrayfun(@(x) x.events.t_targ+0.3, trials_behv,'un',0);
        t_stop = arrayfun(@(x) x.events.t_stop, trials_behv,'un',0);
        t_sac = arrayfun(@(x) x.events.t_sac, trials_behv,'un',0);
        sac_mag = arrayfun(@(x) x.events.sac_mag, trials_behv,'un',0);
        sac_dir = arrayfun(@(x) x.events.sac_dir, trials_behv,'un',0);
        keepindx = cellfun(@(x,tmin,tmax) (x>=tmin & x<=tmax), t_sac,t_targOFF,t_stop,'un',0); 
        t_sac = cellfun(@(x,i) x(i), t_sac,keepindx,'un',0);  
        
        sac_mag = cellfun(@(x,i) x(i), sac_mag,keepindx,'un',0);  
        sac_dir = cellfun(@(x,i) x(i), sac_dir,keepindx,'un',0);  
        
        % first and last saccades always ruin things, remove them
        if~isempty(t_sac{1});   t_sac{1}(1) = [];  sac_mag{1}(1) = []; sac_dir{1}(1) = []; end
        if~isempty(t_sac{end}); t_sac{end}(end) = []; sac_mag{end}(end) = []; sac_dir{end}(end) = []; end 
        
        % find saccade time points
        sacindx = cellfun(@(x,t) find(hist(x,t)),t_sac,ts,'UniformOutput',false);
        
        pre_tte = []; post_tte = [];
        for n = 1:numel(trials_behv)
        pre_tte{n} = arrayfun(@(t) nanmean(ttemag{n}(t-presac:t-2)),sacindx{n});  
        post_tte{n} = arrayfun(@(t) nanmean(ttemag{n}(t+2:t+presac)),sacindx{n}); 
%         nanindx = isnan(pre_tte{n}) | isnan(post_tte{n});
%         if nanindx; post_tte{n}(nanindx) = []; pre_tte{n}(nanindx) = []; end
        end
        
        decrindx = cellfun(@(x,y) x > y, pre_tte,post_tte,'un',0); 
        incrindx = cellfun(@(x,y) x < y, pre_tte,post_tte,'un',0); 
        
        temp_trials_de = trials_behv;
        temp_trials_in = trials_behv;
        
        for n = 1:length(trials_behv)
            if ~isempty(t_sac{n})
            % keep saccades that decreased error
            temp_trials_de(n).events.t_sac = t_sac{n}(decrindx{n});
            temp_trials_de(n).events.sac_mag = sac_mag{n}(decrindx{n});
            temp_trials_de(n).events.sac_dir = sac_dir{n}(decrindx{n});
            % keep saccades that increased error
            temp_trials_in(n).events.t_sac = t_sac{n}(incrindx{n});
            temp_trials_in(n).events.sac_mag = sac_mag{n}(incrindx{n});
            temp_trials_in(n).events.sac_dir = sac_dir{n}(incrindx{n});
            else
            temp_trials_de(n).events.t_sac = [];
            temp_trials_de(n).events.sac_mag = [];
            temp_trials_de(n).events.sac_dir = [];
            
            temp_trials_in(n).events.t_sac = [];
            temp_trials_in(n).events.sac_mag = [];
            temp_trials_in(n).events.sac_dir = [];
            end
        end
        
        stats_de = AnalyseUnit(trials_spks,temp_trials_de,behv_stats,lfps,prs);
        stats_in = AnalyseUnit(trials_spks,temp_trials_in,behv_stats,lfps,prs);
        
        su_sac_tte(s).units(i).stats.decrease = stats_de;
        su_sac_tte(s).units(i).stats.increase = stats_in;
        
        if 0
        % plot saccade rate
        behv.stats = su(s).behv_stats;
        behv.trials = temp_trials_de;
        unit.trials = su(s).units.trials;
        unit.stats = stats_de;
        trial_type = 'all';        
        PlotUnit_Akis(prs,behv,unit,'rate_sac',trial_type); hold on;% decreased error
        
        behv.trials = temp_trials_in;
        unit.stats = stats_in;
        PlotUnit_Akis(prs,behv,unit,'rate_sac',trial_type); % increased error
        
        legend({'decreased TTE','increased TTE'});
        sgtitle([su(s).session ', \bfarea:\rm ' su(s).prs.area{:} ', \bfunit:\rm ' num2str(i)]);

        % plot tuning functions
        figure;
        behv.trials = temp_trials_de;
        unit.stats = stats_de;
        PlotUnit_Akis(prs,behv,unit,'tuning_continuous',trial_type); % decreased error
        
        behv.trials = temp_trials_in;
        unit.stats = stats_in;
        PlotUnit_Akis(prs,behv,unit,'tuning_continuous',trial_type); % increased error
        
        legend({'decreased TTE','increased TTE'});
        sgtitle([su(s).session ', \bfarea:\rm ' su(s).prs.area{:} ', \bfunit:\rm ' num2str(i)]);
        end      
    end
end

disp('........saving file su_saccade_effect_tte.mat');
cd('C:\Data\firefly_analysis\8ar');
save('su_saccade_effect_tte.mat','su_sac_tte','-v7.3');
disp('File saved!');

%% Fit GAM model (PYTHON)


%% Get tuning functions from Python fits
name_pattern = 'm71s*_fit_info_results.mat';
fits_path = 'G:\My Drive\MATLAB\Code\Savin-Angelaki\PGAM\fitting';
pgam = GetGAMfits(fits_path,name_pattern);

nrows = 5;

% Plot units
Nsessions = numel(pgam);
for n = 1:Nsessions
    
    if 1
        Nunits = numel(pgam(n).units);
        for k = 1:Nunits
            
            if strcmp(pgam(n).units(k).info(1).unit_type,'singleunit') && strcmp(pgam(n).units(k).info(1).brain_area,'8ar') 
                neuron = pgam(n).units(k).info(1).neuron;
                
                figure('Name',['session: ' pgam(n).units(k).info(1).session ', unit: ' num2str(neuron)],'NumberTitle','off');
                set(gcf, 'Position', [674 511 560 420]) % [1322 511 560 420]
                
                Nvar = numel(pgam(n).units(k).info);
                ncol = ceil(Nvar/nrows);
                for v = 1:Nvar-1 % skip spike_hist
                    if pgam(n).units(k).info(v).pval < 0.001; colr = 'r'; else; colr = 'k'; end
                    
                    subplot(nrows,ncol,v); hold on;
                    %             plot(pgam(n).units(k).info(v).x_tuning,pgam(n).units(k).info(v).tuning_function,colr)

                    varname = pgam(n).units(k).info(v).variable;
                    tuning_function = pgam(n).units(k).info(v).tuning_function;
                    fplus = pgam(n).units(k).info(v).fplus;
                    fminus = pgam(n).units(k).info(v).fminus;
                    len = numel(tuning_function);
                    
                    % temp fix
                    if contains(varname,'t_')
%                         offset0 = -2:.001:2;
%                         fval = arrayfun(@(x) fix_pgam_CI(x,tuning_function,fplus,fminus),offset0);
%                         [~,ind] = min(fval);
%                         offset = offset0(ind);
%                         tuning_function = tuning_function + offset;
                        
                        mvpos = 0.5*(fplus(round(len/2)) + fminus(round(len/2)));
                        offset = mvpos - tuning_function(round(len/2));
                        tuning_function = tuning_function + offset;
                    end
                    fplus = fplus - tuning_function;
                    fminus = fminus - tuning_function;
                    
                shadedErrorBar(pgam(n).units(k).info(v).x_tuning,tuning_function,...
                    [fplus ; -fminus],'lineprops',{'color',colr})
                title(pgam(n).units(k).info(v).variable);
                end
                suptitle(['session: ' pgam(n).units(k).info(v).session ', unit: ' num2str(neuron)]); 
                pause(1.5);
            end
        end
    end
end

%% Compare fraction tuned across setups **************************
unityVR = arrayfun(@(x) x.prs.unityVR, su)';
brain_areas = unique(arrayfun(@(x) x.prs.area, su));
setups = {'unity','spike2'};

rec_area = cell2mat(arrayfun(@(x) strcmp(brain_areas,x.prs.area)', su,'un',0))';

selectindx = [{unityVR & rec_area(:,1)} {~unityVR & rec_area(:,1)} {unityVR & rec_area(:,2)} {~unityVR & rec_area(:,2)}];
vals = {'unity+7a','spike2+7a','unity+8ar','spike2+8ar'};

varnames = {pgam(1).units(1).info.variable}; varnames(end) = []; % remove spike_hist

name_pattern = 'm*s*_fit_info_results_tte_mag.mat';
fits_path = 'G:\My Drive\MATLAB\Code\Savin-Angelaki\PGAM\fitting';
pgam = GetGAMfits(fits_path,name_pattern);

fr_tuned = [];
for k = 1:numel(selectindx)  

    pgam_temp = pgam(selectindx{k});
    units = [pgam_temp.units];
        
    for n = 1:length(varnames)
        fr_tuned(k).val = vals{k};
        
        nvar = arrayfun(@(x) numel(x.info),units);
        
        keepindx = nvar == max(nvar);
        fr_tuned(k).vars.(varnames{n}).pval = cell2mat(arrayfun(@(x) x.info(n).pval, units(keepindx),'un',0));
        fr_tuned(k).vars.(varnames{n}).tuned =  fr_tuned(k).vars.(varnames{n}).pval < 0.001;
    end

end

% plot all
nrows = 3;
ncol = ceil(numel(varnames)/nrows);
colr = brewermap(2,'Set1'); colr = repmat(colr,2,1); colr = diag([0.5 0.5 1 1])*colr;

figure;
fr = [];
for n = 1:length(varnames)
    subplot(nrows,ncol,n); hold on;
    for k = 1:numel(fr_tuned)
        fr.(varnames{n})(k) = sum(fr_tuned(k).vars.(varnames{n}).tuned)/numel(fr_tuned(k).vars.(varnames{n}).tuned);
        
        bar(k,fr.(varnames{n})(k),'facecolor',colr(k,:));
    end
    xticks(1:numel(fr_tuned)); xticklabels({fr_tuned.val}); xtickangle(90);
    ylabel('% tuned'); title(varnames{n});
end
sgtitle('Fraction tuned comparison across setups');

% compare fractions between setups (X^2 test)
chi2stat = []; p = [];
for n = 1:numel(varnames)
    cnt=1;
    for k = 1:2:numel(fr_tuned)
        % Observed data
        n1 = sum(fr_tuned(k).vars.(varnames{n}).tuned); N1 = numel(fr_tuned(k).vars.(varnames{n}).tuned);
        n2 = sum(fr_tuned(k+1).vars.(varnames{n}).tuned); N2 = numel(fr_tuned(k+1).vars.(varnames{n}).tuned);
        % Pooled estimate of proportion
        p0 = (n1+n2) / (N1+N2);
        % Expected counts under H0 (null hypothesis)
        n10 = N1 * p0;
        n20 = N2 * p0;
        % Chi-square test, by hand
        observed = [n1 N1-n1 n2 N2-n2];
        expected = [n10 N1-n10 n20 N2-n20];
        chi2stat.(varnames{n})(cnt) = sum((observed-expected).^2 ./ expected);
        p.(varnames{n})(cnt) = 1 - chi2cdf(chi2stat.(varnames{n})(cnt),1);
        
        cnt=cnt+1;
    end
end
chi2significant = structfun(@(x) x < 0.05, p,'un',0) % 1st column: 7a, 2nd column: 8ar

%% Tuning analysis ************************************
name_pattern = 'm71s*_fit_info_results_tte_mag.mat';
fits_path = 'G:\My Drive\MATLAB\Code\Savin-Angelaki\PGAM\fitting';
pgam = GetGAMfits(fits_path,name_pattern);
Nsessions = numel(pgam);

varnames = {pgam(1).units(1).info.variable}; % rmindx = cellfun(@(x) strcmp(x,'spike_hist'),varnames); varnames(rmindx) = [];
Nvars = numel(varnames);

var_label.sensorimotor = {'rad_vel','ang_vel','rad_acc','ang_acc','t_move','t_stop','t_flyOFF'};
var_label.latent = {'rad_target','ang_target','rad_path','ang_path','tte_mag'};
var_label.lfp = {'lfp_theta','lfp_alpha','lfp_beta'};
var_label.other = {'t_reward','t_sac','eye_vert','eye_hori'};
    
units = [pgam.units];
Nunits = numel(units);

% correct PPC label to 7a
brain_area = arrayfun(@(x) x.info(1).brain_area,units,'un',0);
for k = 1:numel(brain_area)
    if strcmp(brain_area{k},'PPC')
        for i = 1:numel(units(k).info); units(k).info(i).brain_area = '7a';  end
    elseif strcmp(brain_area{k},'PFC')
        for i = 1:numel(units(k).info); units(k).info(i).brain_area = 'dlPFC';  end
    end
end

areas = unique(arrayfun(@(x) x.brain_area, tuning_analysis.var_tuned.unit,'un',0));

a = 0.001; % statistical significance level

tuning_analysis = [];
% number of variables a unit is tuned to
for n = 1:Nunits
    for v = 1:Nvars
        
        varind = arrayfun(@(x) strcmp(x.variable,varnames{v}), units(n).info);
        
        if any(varind)
            tuning_analysis.var_tuned.unit(n).val(v) = double(units(n).info(varind).pval < a);
            tuning_analysis.var_tuned.unit(n).var(v) = {units(n).info(varind).variable};
        else
            tuning_analysis.var_tuned.unit(n).val(v) = nan;
            tuning_analysis.var_tuned.unit(n).var(v) = varnames(v);
        end
        
    end
    tuning_analysis.var_tuned.unit(n).brain_area = units(n).info(1).brain_area;
    tuning_analysis.var_tuned.unit(n).unit_type = units(n).info(1).unit_type;    
    tuning_analysis.var_tuned.unit(n).session = units(n).info(1).session;

end


% fraction of units tuned to a variable

for n = 1:Nunits
    for v = 1:Nvars
        
        unitindx = strcmp({units(n).info.variable},varnames{v});
        
        unitvars = {units(n).info.variable};
        varindx = cellfun(@(x) strcmp(x,varnames{v}),unitvars);
        
        if sum(varindx)==1
            tuning_analysis.fr_tuned.(varnames{v}).val(n) = double(units(n).info(varindx).pval < a);
            tuning_analysis.fr_tuned.(varnames{v}).brain_area{n} = units(n).info(varindx).brain_area;
            tuning_analysis.fr_tuned.(varnames{v}).unit_type{n} = units(n).info(varindx).unit_type;
            tuning_analysis.fr_tuned.(varnames{v}).session{n} = units(n).info(varindx).session;
        else
            tuning_analysis.fr_tuned.(varnames{v}).val(n) = nan;
            tuning_analysis.fr_tuned.(varnames{v}).brain_area{n} = nan;
            tuning_analysis.fr_tuned.(varnames{v}).unit_type{n} = nan;
            tuning_analysis.fr_tuned.(varnames{v}).session{n} = nan;
        end
        
    end
end

% plots
figure;
nx = 0:18;

varlabels = fieldnames(var_label);
for k = 1:numel(varlabels)
    subplot(2,2,k); hold on;
    
    varindx = ismember(varnames,var_label.(varlabels{k}));
        
    for n = 1:numel(areas)
        indx_area = arrayfun(@(x) strcmp(x.brain_area,areas{n}), tuning_analysis.var_tuned.unit);
        
        temp_area = arrayfun(@(x) sum(x.val(varindx)), tuning_analysis.var_tuned.unit(indx_area));
        [hist_area{n},nx] = hist(temp_area,nx); hist_area{n} = hist_area{n}/sum(hist_area{n});
    end
    hist_area = cellfun(@(x) x',hist_area,'un',0);
    bar(nx,[hist_area{:}]); xlabel('# tuned variables'); ylabel('probability'); xlim([-1 8]); title(varlabels{k});

end
legend(areas); suptitle('Number of variables a unit is tuned to - Viktor (Unity)');%suptitle('Number of variables a unit is tuned to');


figure;
colr = brewermap(numel(areas),'*Set1');

for k = 1:numel(varlabels)
    subplot(2,2,k); hold on;
    
    varindx = find(ismember(varnames,var_label.(varlabels{k})));
    
    groupvars = var_label.(varlabels{k});
    rmindx = [];
    for j = 1:numel(groupvars)
        jj = cellfun(@(x) strcmp(groupvars(j),x), varnames);
        
        if any(jj)
            
            for n = 1:numel(areas)
                
                indx_area = arrayfun(@(x) strcmp(x,areas{n}), tuning_analysis.fr_tuned.(varnames{jj}).brain_area);
                
                temp_area = tuning_analysis.fr_tuned.(varnames{jj}).val(indx_area);
                                
                plot(j+0.1*(n-.5*numel(areas)),sum(temp_area)/numel(temp_area),'wo','markerfacecolor',colr(n,:));
            end
        end
        ylabel('fraction tuned'); title(varlabels{k});
    end
    xticks(1:numel(groupvars)); xticklabels(var_label.(varlabels{k})); xtickangle(90); axis([0 numel(groupvars)+1 0 1]);
end
legend(areas); suptitle('Fraction of neurons tuned to each variable - Viktor (Unity)');%suptitle('Fraction of neurons tuned to each variable');

%% Find sessions with > 15 sua+mua for GPFA analysis

fits_path = 'G:\My Drive\MATLAB\Code\Savin-Angelaki\PGAM\fitting';
cd(fits_path)
fitfile = dir(fullfile(fits_path,'m71s*_fit_info_results*.mat'));
sessnum = arrayfun(@(x) strsplit(x.name,{'s','_'}),fitfile,'un',0); sessnum = cellfun(@(x) str2double(x{2}),sessnum);
[~,I] = sort(sessnum);
fitfile = fitfile(I);
Nsessions = numel(fitfile);

cnt=1;
qualifying_session = [];
for k = 1:Nsessions
    temp = importdata(fitfile(k).name);
    units = unique([temp.neuron]);
    Nunits(k) = numel(units);
    Nvar = numel(temp)/Nunits(k);
    unit_type = {temp.unit_type};
    brain_area = {temp.brain_area};
    areas = unique(brain_area);
    
    for n = 1:numel(areas)
        
        sessname = strsplit(fitfile(k).name,{'_'}); sessname = sessname{1};
        singleunits = sum(cellfun(@(x,y) strcmp(x,'singleunit') & strcmp(y,areas{n}), unit_type, brain_area))/Nvar;
        multiunits = sum(cellfun(@(x,y) strcmp(x,'multiunit') & strcmp(y,areas{n}), unit_type, brain_area))/Nvar;
        
        totalunits = sum(singleunits + multiunits);
        if totalunits >= 15
            qualifying_session(cnt).session = sessname;
            qualifying_session(cnt).totalunits = totalunits;
            qualifying_session(cnt).area = areas{n};
            
            cnt=cnt+1;
            
        end
        
    end
end

for n = 1:numel(areas)
indx_area = arrayfun(@(x) strcmp(x.area,areas{n}),qualifying_session);
sessions_area = cellfun(@(x,y,z) sprintf('%s, session: %s, Nunits: %d ',x,y,z), {qualifying_session(indx_area).area}, {qualifying_session(indx_area).session}, {qualifying_session(indx_area).totalunits},'un',0)'
end

%% Compare pseudo-R^2 between fits with and without t_sac
fits_path = 'G:\My Drive\MATLAB\Code\Savin-Angelaki\PGAM\fitting';
cd(fits_path)
% saccade fit files
fitfile_sac = dir(fullfile(fits_path,'m*s*_fit_info_results_t_sac.mat'));
sessnum_sac = arrayfun(@(x) strsplit(x.name,{'s','_'}),fitfile_sac,'un',0); sessnum_sac = cellfun(@(x) str2double(x{2}),sessnum_sac);
[~,I] = sort(sessnum_sac);
fitfile_sac = fitfile_sac(I);
Nsessions = numel(fitfile_sac);

% no-saccade fit files
fitfile = dir(fullfile(fits_path,'m*s*_fit_info_results.mat'));
sessnum = arrayfun(@(x) strsplit(x.name,{'s','_'}),fitfile,'un',0); sessnum = cellfun(@(x) str2double(x{2}),sessnum);
sessindx = ismember(sessnum,sessnum_sac);
sessnum = sessnum(sessindx);
fitfile = fitfile(sessindx);
[~,I] = sort(sessnum);
fitfile = fitfile(I);

figure;
pseudo_r2 = [];
for n = 1:Nsessions
    % saccade fits
    temp_sac = importdata(fitfile_sac(n).name);
    sessname = strsplit(fitfile_sac(n).name,{'_'}); sessname = sessname{1};
    units_sac = unique([temp_sac.neuron]);
    Nunits_sac = numel(units_sac);
    Nvar_sac = numel(temp_sac)/Nunits_sac;
    unit_type_sac = {temp_sac.unit_type};
    brain_area = unique([temp_sac.brain_area]);
    pseudo_r2_temp = [temp_sac.pseudo_r2]; pseudo_r2_temp = pseudo_r2_temp(1:Nvar_sac:end);
    
    pseudo_r2.t_sac.val = pseudo_r2_temp;
    pseudo_r2.t_sac.brain_area = brain_area;
    pseudo_r2.t_sac.session = sessname;
    pseudo_r2.t_sac.sua = cellfun(@(x) strcmp(x,'singleunit'), unit_type_sac(1:Nvar_sac:end));
    pseudo_r2.t_sac.mua = cellfun(@(x) strcmp(x,'multiunit'), unit_type_sac(1:Nvar_sac:end));
    
    % no-saccade fits
    temp = importdata(fitfile(n).name);
    units = unique([temp.neuron]);
    Nunits = numel(units);
    Nvar = numel(temp)/Nunits;
    unit_type = {temp.unit_type};
    brain_area = unique([temp.brain_area]);
    pseudo_r2_temp = [temp.pseudo_r2]; pseudo_r2_temp = pseudo_r2_temp(1:Nvar:end);

    pseudo_r2.nosac.val = pseudo_r2_temp;
    pseudo_r2.nosac.brain_area = brain_area;
    pseudo_r2.nosac.session = sessname;
    pseudo_r2.nosac.sua = cellfun(@(x) strcmp(x,'singleunit'), unit_type(1:Nvar:end));
    pseudo_r2.nosac.mua = cellfun(@(x) strcmp(x,'multiunit'), unit_type(1:Nvar:end));

    % match units before plotting
    unitindx_sac = ismember(units_sac,units);
    unitindx = ismember(units,units_sac);
    
    % Plot
    subplot(1,2,1); hold on;
    plot(pseudo_r2.nosac.val(pseudo_r2.nosac.sua & unitindx),pseudo_r2.t_sac.val(pseudo_r2.t_sac.sua & unitindx_sac),'o'); plot(-1:1,-1:1,'k--');
    xlabel('pseudo-R^2 (no t_{sac})'); ylabel('pseudo-R^2 (with t_{sac})'); title('single units'); axis([-0.05 0.6 -0.05 0.6]);
    
    subplot(1,2,2); hold on;
    plot(pseudo_r2.nosac.val(pseudo_r2.nosac.mua & unitindx),pseudo_r2.t_sac.val(pseudo_r2.t_sac.mua & unitindx_sac),'o'); plot(-1:1,-1:1,'k--');
    xlabel('pseudo-R^2 (no t_{sac})'); ylabel('pseudo-R^2 (with t_{sac})'); title('multi units'); axis([-0.05 0.6 -0.05 0.6]);
end


%%
for s = 1:numel(su)
    
    su(s).prs.tuning_continuous = {'v','w','r_targ','theta_targ','d','phi','eye_ver','eye_hor','eye_verhor','r','theta','rtheta','phase','h1','h2','r_accel','theta_accel','sacmag','sacdir','sacmagdir','targ_ver','targ_hor','targ_verhor','tte_ver','tte_hor','tte_verhor'};
end



