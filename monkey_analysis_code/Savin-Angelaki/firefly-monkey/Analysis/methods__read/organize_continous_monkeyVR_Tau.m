function [outdata] = organize_continous_monkeyVR_Tau(data)

scale = 10; % this is hard coded for now, but must be changed

% data(data(:, 1) == 0, :) = []; 

outdata.trial_num =         data(:, 1); 
outdata.trial_time =        data(:, 2);
outdata.phase =             data(:, 3); 
outdata.on_off =            data(:, 4); 
outdata.posX =              data(:, 5)/scale; % in meters!!!
outdata.posY =              data(:, 6)/scale; 
outdata.posZ =              data(:, 7)/scale; 
outdata.rotX =              data(:, 8); 
outdata.rotY =              data(:, 9); 
outdata.rotZ =              data(:, 10); 
outdata.rotW =              data(:, 11); 
outdata.FFx =               data(:, 12)/scale; 
outdata.FFy =               data(:, 13)/scale; 
outdata.FFz =               data(:, 14)/scale;
outdata.FFvel =             data(:, 15)/scale; 

outdata.eyetype =           data(:, 16);
outdata.confidence =        data(:, 17);
outdata.Gx =                data(:, 18); % left/right
outdata.Gy =                data(:, 19); % up/down
outdata.Gz =                data(:, 20); % back/forth
outdata.conv_dist =         data(:, 21)/scale; % double check this

outdata.RXctr =             data(:, 22);
outdata.RYctr =             data(:, 23);
outdata.RZctr =             data(:, 24);
outdata.LXctr =             data(:, 25);
outdata.LYctr =             data(:, 26);
outdata.LZctr =             data(:, 27);

outdata.RXnorm =            data(:, 28);
outdata.RYnorm =            data(:, 29);
outdata.RZnorm =            data(:, 30);
outdata.LXnorm =            data(:, 31);
outdata.LYnorm =            data(:, 32);
outdata.LZnorm =            data(:, 33);

outdata.ksi_lin_vel =       data(:, 34)/scale; % first filtered process noise velocity
outdata.eta_lin_vel =       data(:, 35)/scale; % re-filtered process noise velocity
outdata.ksi_ang_vel =       data(:, 36);
outdata.eta_ang_vel =       data(:, 37);
outdata.linear_velocity =   data(:, 38)/scale; % clean+process noise velocity (i.e. final velocity)
outdata.angular_velocity =  data(:, 39);
outdata.clean_lin_vel =     data(:, 40)/scale; % clean+process noise velocity (i.e. final velocity)
outdata.clean_ang_vel =     data(:, 41);
outdata.raw_lin_js =        data(:, 42); % raw joystick input
outdata.raw_ang_js =        data(:, 43);


% outdata.hitX =              data(:, 24)/scale;  % hit = gaze + distance
% outdata.hitY =              data(:, 25); 
% outdata.hitZ =              data(:, 26)/scale; 
% outdata.conv_dist =         data(:, 27)/scale;
% outdata.LPD =               data(:, 28); 
% outdata.RPD =               data(:, 29); 
% outdata.Lopen =             data(:, 30); 
% outdata.Ropen =             data(:, 31); 

end

