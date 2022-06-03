function [outdata] = organize_discontinuous_human_VR(data)

scale = 10; 


outdata.trial_num =         data(:, 1); 

outdata.maxV =              data(:, 2)/scale; 
outdata.maxW =              data(:, 3); 

outdata.ffv =               data(:, 4)/scale; 
outdata.ff_duration =       data(:, 5); 
% outdata.aswer =             data(:, 6); 

outdata.PosXo =             data(:, 7)/scale; 
outdata.PosYo =             data(:, 8)/scale; 
outdata.PosZo =             data(:, 9)/scale; 

outdata.RotXo =             data(:, 10); 
outdata.RotYo =             data(:, 11); 
outdata.RotZo =             data(:, 12); 
outdata.RotWo =             data(:, 13); 

outdata.FFx =               data(:, 14)/scale; 
outdata.FFy =               data(:, 15)/scale; 
outdata.FFz =               data(:, 16)/scale; 

outdata.pcheckX =           data(:, 17); % position at stop
outdata.pcheckY =           data(:, 18); 
outdata.pcheckZ =           data(:, 19); 

outdata.rcheckX =           data(:, 20);
outdata.rcheckY =           data(:, 21); % heading at stop 
outdata.rcheckZ =           data(:, 22); 
outdata.rcheckW =           data(:, 23); 

outdata.distToFF =          data(:, 24)/scale; 
outdata.rewarded =          data(:, 25); 
outdata.timeout =           data(:, 26); 
outdata.beginTime =         data(:, 27); % trial start
outdata.checkTime =         data(:, 28); % monkey stop time
outdata.endTime =           data(:, 29); 
outdata.waitTime =          data(:, 30); % distance checking interval (not used)
outdata.ITI =               data(:, 31);

outdata.tau =               data(:, 32);
outdata.type =              data(:, 34); % Type: 0) discrete, 1) continuous, 2) none
outdata.noisetautau =       data(:, 35);
outdata.noisetau =          data(:, 36);
outdata.gainw =             data(:, 37);
outdata.gainv =             data(:, 38); % /scale?
outdata.ntaus =             data(:, 39);
outdata.mintau =            data(:, 40);
outdata.maxtau =            data(:, 41);
outdata.x =                 data(:, 42)/scale;
outdata.T =                 data(:, 43);
outdata.vthresh =           data(:, 44)/scale;
outdata.wthresh =           data(:, 45);

end

