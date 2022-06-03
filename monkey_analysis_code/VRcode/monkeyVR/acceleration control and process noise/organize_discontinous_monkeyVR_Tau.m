function [outdata] = organize_discontinous_monkeyVR_Tau(data)

scale = 10; 


outdata.trial_num =         data(:, 1); 

outdata.maxV =              data(:, 2)/scale; 
outdata.maxW =              data(:, 3); 

outdata.ffv =               data(:, 4)/scale; 
outdata.ff_duration =       data(:, 5); 
outdata.floordensity =      data(:, 6); 

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
outdata.rewardDur =         data(:, 27); 
outdata.beginTime =         data(:, 28); % trial start
outdata.checkTime =         data(:, 29); % monkey stop time
outdata.rewardTime =        data(:, 30); 
outdata.endTime =           data(:, 31); 
outdata.waitTime =          data(:, 32); % distance checking interval (not used)
outdata.ITI =               data(:, 33);


outdata.tau =               data(:, 34);
outdata.type =              data(:, 35); % Type: 0) discrete, 1) continuous, 2) none
outdata.TauTau =            data(:, 36);
outdata.NoiseTau =          data(:, 37);
outdata.gainw =             data(:, 38);
outdata.gainv =             data(:, 39); % /scale?
outdata.ntaus =             data(:, 40);
outdata.mintau =            data(:, 41);
outdata.maxtau =            data(:, 42);
outdata.x =                 data(:, 43)/scale;
outdata.T =                 data(:, 44);
outdata.vthresh =           data(:, 45)*scale; % check why this needs to be multipied and not divided by scale....
outdata.wthresh =           data(:, 46);


end

