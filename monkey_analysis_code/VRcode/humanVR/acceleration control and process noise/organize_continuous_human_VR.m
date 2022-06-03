function [outdata] = organize_continuous_human_VR(data)

scale = 10; % this is hard coded for now, but must be changed

% remove data before first trial starts
data(data(:, 1) == 0, :) = []; 

outdata.trial_num =         data(:, 1); 
outdata.trial_time =        data(:, 2) - data(1,2);
outdata.phase =             data(:, 3); 
outdata.on_off =            data(:, 4); 
outdata.posX =              data(:, 5)/scale; % in meters!!!
outdata.posY =              data(:, 6)/scale; 
outdata.posZ =              data(:, 7)/scale; 
outdata.rotX =              data(:, 8); 
outdata.rotY =              data(:, 9); 
outdata.rotZ =              data(:, 10); 
outdata.rotW =              data(:, 11); 
outdata.clean_lin_vel =     data(:, 12)/scale; 
outdata.clean_ang_vel =     data(:, 13); 
outdata.FFx =               data(:, 14)/scale; 
outdata.FFy =               data(:, 15)/scale; 
outdata.FFz =               data(:, 16)/scale;
outdata.FFvel =             data(:, 17)/scale; 

outdata.Gx =                data(:, 18); % left/right
outdata.Gy =                data(:, 19); % up/down
outdata.Gz =                data(:, 20); % back/forth
outdata.Gx0 =               data(:, 21); % left/right
outdata.Gy0 =               data(:, 22); % up/down
outdata.Gz0 =               data(:, 23); % back/forth
outdata.hitX =              data(:, 24)/scale;  % hit = gaze + distance
outdata.hitY =              data(:, 25); 
outdata.hitZ =              data(:, 26)/scale; 
outdata.conv_dist =         data(:, 27)/scale;
outdata.LPD =               data(:, 28); 
outdata.RPD =               data(:, 29); 
outdata.Lopen =             data(:, 30); 
outdata.Ropen =             data(:, 31); 

outdata.ksi_lin_vel =       data(:, 32)/scale; % first filtered process noise velocity
outdata.eta_lin_vel =       data(:, 33)/scale; % re-filtered process noise velocity
outdata.ksi_ang_vel =       data(:, 34);
outdata.eta_ang_vel =       data(:, 35);
outdata.linear_velocity =   data(:, 36)/scale; % clean+process noise velocity
outdata.angular_velocity =  data(:, 37)/scale;
outdata.raw_lin_js =        data(:, 38); % raw joystick input
outdata.raw_ang_js =        data(:, 39);

end

