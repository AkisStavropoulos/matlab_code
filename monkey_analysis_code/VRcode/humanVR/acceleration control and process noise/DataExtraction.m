%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Akis Stavropoulos - July 2021
%
% First try of the acceleration control and process noise paradigms
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Importing and nonsense

% housekeeping
close all; clear all; 

% Path stuff
% addpath('Y:\Projects\Monkey_Unity_VR\data');
addpath('Z:\Users\Akis\Process noise tests');
addpath(genpath('G:\My Drive\MATLAB\Code\monkeyVR\acceleration control and process noise'));

% color stuff
reds = cbrewer('seq', 'Reds', 10);

% Directory stuff
% human_directory = 'Y:\Projects\Monkey_Unity_VR\data';
human_directory = 'Z:\Users\Akis';

%% THINGS TO CHANGE

monkey_name = 'Process noise tests';
% date_name = '20210629';
date_name = '';

% Booleans
boolean_plot =0;
boolean_save_movie =0;
boolean_save = 0; 


%%
monkey_dir = dir(fullfile(human_directory, monkey_name, date_name));
monkey_dir = monkey_dir(~ismember({monkey_dir.name},{'.','..', '.DS_Store'})); % getting rid of nonesense
monkey_dir(~arrayfun(@(x) x.isdir,monkey_dir)) = [];

for i = 2%:size(monkey_dir, 1)
    clear d data dis dis_data
    
    % containers
    
    master_position = [];
    master_polar = [];
    master_cart = [];
    %%
    
    d = dir(fullfile(monkey_dir(i).folder, monkey_dir(i).name, 'continuous*'));
    data = import_continuous_human_VR(fullfile(d.folder, d.name));
    data = organize_continuous_human_VR(data);
    
    
    dis = dir(fullfile(monkey_dir(i).folder, monkey_dir(i).name, 'discontinuous*'));
    dis_data = import_discontinuous_human_VR(fullfile(dis.folder, dis.name));
    dis_data = organize_discontinuous_human_VR(dis_data);
    
    data.rewarded = dis_data.rewarded;
    data = addEvents(data);
    
    
    [data, dis_data] = check_sizes_monkey_VR(data, dis_data);
    data.rewarded = dis_data.rewarded;
    
    
    
    %% For plotting in case we want to see what is happening...
    if boolean_plot
        figure;
        for j = 1:size(data.start_trial, 1)
            plot(data.posX( data.start_trial(j):data.stop_trial(j)), data.posZ(data.start_trial(j):data.stop_trial(j)), '--k'); hold on;
            plot(data.posX( data.start_trial(j):data.start_trial(j)+floor(dis_data.ff_duration(j)*90)), data.posZ(data.start_trial(j):data.start_trial(j)+floor(dis_data.ff_duration(j)*90)), 'k', 'LineWidth',4); hold on;
            scatter(dis_data.FFx(j), dis_data.FFz(j),100, 'r','filled'); hold on;
            xlim([-2 2]);
            ylim([0 6]);
            
            % Scatter plot for initial position of FF
            % scatter plot for end position of human (color codded by
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
    
    %% Resample data to a fixed sampling rate
    dt = 1/90; % mean(diff(data.trial_time))
    t_beg1 = dis_data.beginTime(1);
    data = ReSample2FixedDt_human(data,dt,t_beg1);
    
    %% Transform data to match Edoardo's code format
    default_prs1;
    prs.dt = dt;
    prs.humanVR = 1;
    trials = VR2Edoardo_dataformat(data,dis_data,prs);
end

%% SAVING
if boolean_save
    save_name = fullfile('Y:\Projects\Monkey_Unity_VR\processed_data', [monkey_name, '_', date_name, '.mat']);
    
    try
        save(save_name, 's', '-v7.3');
    catch
        mkdir('Y:\Projects\Monkey_Unity_VR\processed_data');
        save(save_name, 's', '-v7.3');
    end
end




