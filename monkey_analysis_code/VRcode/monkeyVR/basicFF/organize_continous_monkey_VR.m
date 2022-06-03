function [outdata] = organize_continous_monkey_VR(data)

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
outdata.linear_velocity =   data(:, 12)/scale; 
outdata.angular_velocity =  data(:, 13); 
outdata.FFx =               data(:, 14)/scale; 
outdata.FFy =               data(:, 15)/scale; 
outdata.FFz =               data(:, 16)/scale;
outdata.FFvel =             data(:, 17)/scale; 

outdata.eyetype =           data(:, 18);
outdata.confidence =        data(:, 19);
outdata.Gx =                data(:, 20); % left/right
outdata.Gy =                data(:, 21); % up/down
outdata.Gz =                data(:, 22); % back/forth
outdata.conv_dist =         data(:, 23)/scale; % double check this

outdata.RXctr =             data(:, 24);
outdata.RYctr =             data(:, 25);
outdata.RZctr =             data(:, 26);
outdata.LXctr =             data(:, 27);
outdata.LYctr =             data(:, 28);
outdata.LZctr =             data(:, 29);

outdata.RXnorm =            data(:, 30);
outdata.RYnorm =            data(:, 31);
outdata.RZnorm =            data(:, 32);
outdata.LXnorm =            data(:, 33);
outdata.LYnorm =            data(:, 34);
outdata.LZnorm =            data(:, 35);


% outdata.hitX =              data(:, 24)/scale;  % hit = gaze + distance
% outdata.hitY =              data(:, 25); 
% outdata.hitZ =              data(:, 26)/scale; 
% outdata.conv_dist =         data(:, 27)/scale;
% outdata.LPD =               data(:, 28); 
% outdata.RPD =               data(:, 29); 
% outdata.Lopen =             data(:, 30); 
% outdata.Ropen =             data(:, 31); 

end

