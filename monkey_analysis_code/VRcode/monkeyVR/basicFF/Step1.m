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
addpath('Z:\Data\Monkey2_newzdrive\Jimmy\U-probe'); % addpath('Y:\Projects\Monkey_Unity_VR\data');
addpath('G:\My Drive\MATLAB\Code\VRcode\monkeyVR\basicFF');
addpath('G:\My Drive\MATLAB\Code\VRcode\monkeyVR\basicFF\cbrewer');

% color stuff
reds = cbrewer('seq', 'Reds', 10);

% Directory stuff
human_directory = 'Y:\Projects\Monkey_Unity_VR\data'; % human_directory = 'Y:\Projects\Monkey_Unity_VR\data';

%% THINGS TO CHANGE


monkey_name = 'testData';
date_name = '20210716\behavioural data';

% Booleans
boolean_plot =0;
boolean_save_movie =0;
boolean_save = 1; 



%%
monkey_dir = dir(fullfile(human_directory,monkey_name, date_name)); % monkey_dir = dir(fullfile(human_directory, monkey_name, date_name));
monkey_dir = monkey_dir(~ismember({monkey_dir.name},{'.','..', '.DS_Store'})); % getting rid of nonesense


for i = 1%:size(monkey_dir, 1)
    clear d data dis dis_data
    
    % containers
    
    master_position = [];
    master_polar = [];
    master_cart = [];
    %%
    
    d = dir(fullfile(monkey_dir(i).folder, monkey_dir(i).name, 'continuous*'));
    data = import_continous_monkey_VR(fullfile(d.folder, d.name));
    data = organize_continous_monkey_VR(data);
    
    
    dis = dir(fullfile(monkey_dir(i).folder, monkey_dir(i).name, 'discontinuous*'));
    dis_data = import_discontinous_monkey_VR(fullfile(dis.folder, dis.name));
    dis_data = organize_discontinous_monkey_VR(dis_data);
    
    data.rewarded = dis_data.rewarded;
    data = addEvents(data);
    
    
    [data, dis_data] = check_sizes_monkey_VR(data, dis_data);
    data.rewarded = dis_data.rewarded;
    
    
    
    %% For plotting in case we want to see what is happening...
    if boolean_plot
        figure;
        for j = 1:size(data.start_trial, 1)
            plot(data.posX( data.start_trial(j):data.stop_trial(j)), data.posZ(data.start_trial(j):data.stop_trial(j)), '--k'); hold on;
            plot(data.posX( data.start_trial(j):floor(dis_data.ff_duration(j)*90)-data.start_trial(j)), data.posZ(data.start_trial(j):floor(dis_data.ff_duration(j)*90)-data.start_trial(j)), 'k', 'LineWidth',4); hold on;
            plot(data.FFx( data.start_trial(j):data.stop_trial(j)), data.FFz(data.start_trial(j):data.stop_trial(j)), '--r'); hold on;
            plot(data.FFx( data.start_trial(j):floor(dis_data.ff_duration(j)*90)-data.start_trial(j)), data.FFz(data.start_trial(j):floor(dis_data.ff_duration(j)*90)-data.start_trial(j)), 'r', 'LineWidth',4); hold on;
            scatter(data.FFx( data.start_trial(j)), data.FFz(data.start_trial(j)), 100, 'r', 'filled'); hold on;
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
    
    %%

    % for TARGET and RESPONSE
    to_del = isnan(data.stop_trial) | isnan(data.start_trial);
    data.stop_trial(to_del) = []; 
    data.start_trial(to_del) = []; 
    data.end_trial(to_del) = []; 
    
    x =  data.FFx(data.stop_trial) -  data.FFx(data.start_trial);
    y =  data.posX(data.stop_trial) - data.FFx(data.start_trial);
    dur = dis_data.ff_duration;
    vel = dis_data.ffv/10; % this
    dur_trial = (data.end_trial - data.start_trial)/90; % and this is hardcoded again, should be changed...
    for k = 1:size(data.start_trial, 1)
        vel_observer(k) = nanmean(data.linear_velocity(data.start_trial(k):data.end_trial(k)));
    end
    
    
    master_position = vertcat(master_position, horzcat(x, y, dur, vel, dur_trial, vel_observer'));
    
    
    % theta and rho stuff...
    [theta_resp,rho_resp] = cart2pol(data.posX(data.stop_trial), data.posZ(data.stop_trial) - data.posZ(data.start_trial));
    [theta_targ,rho_targ] = cart2pol(data.FFx(data.stop_trial), data.FFz(data.stop_trial) - data.posZ(data.start_trial));
    [theta_targ_start,rho_targ_start] = cart2pol(data.FFx(data.start_trial), data.FFz(data.start_trial));
    
    master_polar = vertcat(master_polar, horzcat(theta_resp, rho_resp, theta_targ, rho_targ, theta_targ_start, rho_targ_start));
    
    
    % posX start, posZ start, posX end, posZ end, FFx start, FFz start, FFx end, FFz end
    master_cart = vertcat(master_cart, ...
        horzcat(data.posX(data.start_trial), data.posZ(data.start_trial),...
        data.posX(data.stop_trial), data.posZ(data.stop_trial),...
        data.FFx(data.start_trial), data.FFz(data.start_trial),...
        data.FFx(data.stop_trial), data.FFz(data.stop_trial)));
    
    
    
    % v_est = v_a (sin(theta) - cos(theta) z_l,ini / z_d)
    
    %     va = vel_observer';
    %     sin_theta = data.FFx(data.stop_trial);
    %     cos_theta = data.FFz(data.start_trial);
    %     zl = data.FFx(data.start_trial);
    %     zd = data.FFz(data.start_trial);
    %
    %     v_est = sin_theta - ((cos_theta .* zl)./zd);
    
    
    %% (Added By Nastaran) other info
    
    % need a vector for on or off trials
    OnOff= dis_data.ff_duration;
    
    %% Save to s
    
    s(i).master_position = master_position;
    s(i).master_polar = master_polar;
    s(i).master_cart = master_cart;
    s(i).continuous = data; % I typically don't save this, as it's big and I really just care about endpoints, but still...
    s(i).OnOff = OnOff;
    
    
    clear master_position master_polar master_cart vel_observer
    
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

%% Resample data to a fixed sampling rate
dt = 1/200;
t_beg1 = dis_data.beginTime(1);
data = ReSample2FixedDt(data,dt,t_beg1);


%% Transform data to match Edoardo's code format

trials = VR2Edoardo_dataformat(data,dis_data);



