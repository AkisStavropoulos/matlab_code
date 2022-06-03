%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Jean-Paul Noel - March 2021
%
% First try at the new Monkey VR with Unity (Josh, Nastaran, Panos)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Importing and nonsense

% housekeeping
close all; clear all; 

% Path stuff
addpath('Z:\Users\Akis\Process noise tests\Monkey setup'); % addpath('Y:\Projects\Monkey_Unity_VR\data');
addpath('G:\My Drive\MATLAB\Code\VRcode\monkeyVR\acceleration control and process noise');

% extraction parameters
default_prs1;

% color stuff
reds = cbrewer('seq', 'Reds', 10);

% Directory stuff
monkey_directory = 'Z:\Users\Akis\Process noise tests\Monkey setup\behavioural data'; % 'Y:\Projects\Monkey_Unity_VR\data'; % 'Z:\Data\Monkey2_newzdrive\Jimmy\U-probe\'; 

%% THINGS TO CHANGE

cd(monkey_directory);

monkey_dir = dir(monkey_directory); % monkey_dir = dir(fullfile(monkey_directory, monkey_name, date_name));
monkey_dir = monkey_dir(contains({monkey_dir.name},{'continuous','discontinuous'})); % getting rid of nonesense

flist_cont = dir('continuous*');
flist_disc = dir('discontinuous*');
if numel(flist_cont) ~= numel(flist_disc); error('Number of continuous and discontinuous files is not the same.'); end

% monkey_name = 'testData';
% date_name = '20210730\behavioural data'; % '20210713\behavioural data';
% monkey_dir = dir(fullfile(monkey_directory, monkey_name, date_name)); % dir(fullfile(monkey_directory, date_name));
% monkey_dir = monkey_dir(~ismember({monkey_dir.name},{'.','..', '.DS_Store'})); % getting rid of nonesense

% Booleans
boolean_plot =0;
boolean_save_movie =0;
boolean_save = 0; 

%%
trials = []; 
for i = 1:size(monkey_dir, 1)
    data = []; dis_data = []; % clear d data dis dis_data
    
    %%
    
%     d = dir(fullfile(monkey_dir(i).folder, monkey_dir(i).name, 'continuous*'));
    data = import_continous_monkeyVR_Tau(fullfile(monkey_directory,flist_cont(i).name));
    data = organize_continous_monkeyVR_Tau(data);
    
    
%     dis = dir(fullfile(monkey_dir(i).folder, monkey_dir(i).name, 'discontinuous*'));
    dis_data = import_discontinous_monkeyVR_Tau(fullfile(monkey_directory,flist_disc(i).name));
    dis_data = organize_discontinous_monkeyVR_Tau(dis_data);
    
    % Resample data to a fixed sampling rate
    dt = 0.006; % matching Kaushik's code, real dt is 1/200;
    t_beg1 = dis_data.beginTime(1);
    if numel(unique(round(diff(data.trial_time),10)))~=1
        data = ReSample2FixedDt(data,dt,t_beg1);
    end
    % Get events
    data.rewarded = dis_data.rewarded; 
    data = addEvents(data);
    
    
    % check sizes
    [data, dis_data] = check_sizes_monkeyVR_Tau(data, dis_data);
    data.rewarded = dis_data.rewarded;
    
    
    %% For plotting in case we want to see what is happening...
    if boolean_plot
        figure;
        for j = 1:size(dis_data.beginTime, 1)
            
            timeindx = find(data.trial_time > dis_data.beginTime(j) & data.trial_time < dis_data.endTime(j));  timeindx = timeindx(1:end-1);
            
            plot(data.posX(timeindx), data.posZ(timeindx), '--k'); hold on;
            plot(data.posX(timeindx(1:floor(dis_data.ff_duration(j)/dt))), data.posZ(timeindx(1:floor(dis_data.ff_duration(j)/dt))), 'k', 'LineWidth',4); hold on;
            plot(data.FFx(timeindx), data.FFz(timeindx), '--r'); hold on;
            plot(data.FFx(timeindx(1:floor(dis_data.ff_duration(j)/dt))), data.FFz(timeindx(1:floor(dis_data.ff_duration(j)/dt))), 'r', 'LineWidth',4); hold on;
            scatter(data.FFx(timeindx(1)), data.FFz(timeindx(1)), 100, 'r', 'filled'); hold on;
            xlim([-2 2]);
            ylim([0 6]);
            
            % Scatter plot for initial position of FF
            % scatter plot for end position of human (color coded by
            % response)
            
            hold off;
            pause(0.25);
            if boolean_save_movie
                F(j) = getframe(gcf) ;
            end
        end
    end
    
    if boolean_save_movie
        % create the video writer with 1 fps
        fname = [monkey_name, '_', date_name, '_', num2str(i)];
        writerObj = VideoWriter(fname, 'MPEG-4');
        writerObj.FrameRate = 10;
        open(writerObj);
        for j=1:length(F)
            % convert the image to a frame
            frame = F(j) ;
            writeVideo(writerObj, frame);
        end
        % close the writer object
        close(writerObj);
    end
           
    % Transform data to match Edoardo's code format

    trials_temp = VR2Edoardo_dataformat(data,dis_data,prs);
    trials = [trials trials_temp];
    
    if 1
        figure;
        for j = 1:numel(trials)
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
end

%% SAVING
if boolean_save
    save_name = fullfile('Y:\Projects\Monkey_Unity_VR\processed_data', [monkey_name, '_', date_name, '.mat']);
    
    try
        save(save_name, 'trials', '-v7.3');
    catch
        mkdir('Y:\Projects\Monkey_Unity_VR\processed_data');
        save(save_name, 'trials', '-v7.3');
    end
end




