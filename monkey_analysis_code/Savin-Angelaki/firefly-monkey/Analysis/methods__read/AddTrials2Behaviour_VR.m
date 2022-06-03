function trials = AddTrials2Behaviour_VR(prs)

%% Importing and nonsense

% Path stuff
addpath('Z:\Data\Monkey2_newzdrive\Jimmy\U-probe'); 
addpath(genpath('G:\My Drive\MATLAB\Code\VRcode\monkeyVR\basicFF'));
addpath('G:\My Drive\MATLAB\Code\VRcode\monkeyVR\');

% Directory stuff
monkey_directory = prs.filepath_behv;
cd(monkey_directory);

monkey_dir = dir(monkey_directory); % monkey_dir = dir(fullfile(monkey_directory, monkey_name, date_name));
monkey_dir = monkey_dir(contains({monkey_dir.name},{'continuous','discontinuous'})); % getting rid of nonesense

flist_cont = dir('continuous*');
flist_disc = dir('discontinuous*');
if numel(flist_cont) ~= numel(flist_disc); error('Number of continuous and discontinuous files is not the same.'); end

trials = []; 
for i = 1:size(flist_cont, 1)
    data = []; dis_data = []; % clear d data dis dis_data
    
    %%
    if ~strcmp(flist_cont(i).name(end-15:end),flist_disc(i).name(end-15:end))
        error('Continuous and Discontinuous data files don''t match!')
    else
        disp(['.............reading ' flist_cont(i).name(end-15:end)]);
    end
    
    if prs.visual_acc_control
    data = import_continous_monkeyVR_Tau(fullfile(monkey_directory,flist_cont(i).name));
    data = organize_continous_monkeyVR_Tau(data);
    
    dis_data = import_discontinous_monkeyVR_Tau(fullfile(monkey_directory,flist_disc(i).name));
    dis_data = organize_discontinous_monkeyVR_Tau(dis_data);        
    elseif prs.basicFF
    data = import_continous_monkey_VR(fullfile(monkey_directory,flist_cont(i).name));
    data = organize_continous_monkey_VR(data);
    
    dis_data = import_discontinous_monkey_VR(fullfile(monkey_directory,flist_disc(i).name));
    dis_data = organize_discontinous_monkey_VR(dis_data);
    end
    
    % Get events
    data.rewarded = dis_data.rewarded;
    data = addEvents(data);
    
    % Resample data to a fixed sampling rate
    dt_original = mean(diff(data.trial_time)); % 0.005;
    dt_new = prs.dt; % matching Kaushik's code 0.006s, real dt is 1/200;
    t_beg1 = dis_data.beginTime(1);
    if numel(unique(round(diff(data.trial_time),10)))~=1
        data = ReSample2FixedDt(data,dt_new,t_beg1);
    end

    
    % check sizes
    [data, dis_data] = check_sizes_monkey_VR(data, dis_data);
    data.rewarded = dis_data.rewarded;
           
    % Transform data to match Edoardo's code format

    trials_temp = VR2Edoardo_dataformat(data,dis_data,prs);
    trials = [trials trials_temp];
    
    if 0
    %% Sanity check
    
    % trajectories
    figure;
    for j = 1:size(dis_data.beginTime, 1)
        
        if 0
        timeindx = find(data.trial_time > dis_data.beginTime(j) & data.trial_time < dis_data.endTime(j));  timeindx = timeindx(1:end-1);
        
        plot(data.posX(timeindx), data.posZ(timeindx), '--k'); hold on;
        plot(data.posX(timeindx(1:floor((dis_data.ff_duration(j))/dt_original))), data.posZ(timeindx(1:floor((dis_data.ff_duration(j))/dt_original))), 'k', 'LineWidth',4); hold on;
        plot(data.FFx(timeindx), data.FFz(timeindx), '--r'); hold on;
        plot(data.FFx(timeindx(1:floor((dis_data.ff_duration(j))/dt_original))), data.FFz(timeindx(1:floor((dis_data.ff_duration(j))/dt_original))), 'r', 'LineWidth',4); hold on;
        scatter(data.FFx(timeindx(1)), data.FFz(timeindx(1)), 100, 'r', 'filled'); hold on;
        xlim([-2 2]);
        ylim([0 6]);
        hold off;
        pause(0.25);
        
        else
        
        timeindx = trials(j).continuous.ts >=0 & trials(j).continuous.ts < trials(j).events.t_end;
        
        plot(trials(j).continuous.xmp(timeindx), trials(j).continuous.ymp(timeindx), '--k'); hold on;
        plot(trials(j).continuous.xfp(timeindx), trials(j).continuous.yfp(timeindx), '--r'); hold on;
        scatter(trials(j).continuous.xfp(find(timeindx,1)), trials(j).continuous.yfp(find(timeindx,1)), 100, 'r', 'filled'); hold on;
        xlim([-200 200]);
        ylim([0 600]);
        hold off;
        pause(0.25);
            
            
        end
    end
    
    
    % scatterplot
    
    
    % eye position
    N = 1300;
    zle = arrayfun(@(x) x.continuous.zle, trials,'un',0);
    zre = arrayfun(@(x) x.continuous.zre, trials,'un',0);
    yle = arrayfun(@(x) x.continuous.yle, trials,'un',0);
    yre = arrayfun(@(x) x.continuous.yre, trials,'un',0);
    zle = cellfun(@(x) [x ; nan(N-numel(x),1)], zle,'un',0);
    zre = cellfun(@(x) [x ; nan(N-numel(x),1)], zre,'un',0);
    yle = cellfun(@(x) abs([x ; nan(N-numel(x),1)]), yle,'un',0);
    yre = cellfun(@(x) abs([x ; nan(N-numel(x),1)]), yre,'un',0);
    zle = [zle{:}];
    zre = [zre{:}];
    yle = [yle{:}];
    yre = [yre{:}];
    figure;plot(nanmedian(zle,2)); hold on;
    plot(nanmedian(zre,2));
    plot(nanmedian(yle,2));
    plot(nanmedian(yre,2));
    
    end
    
end

end

