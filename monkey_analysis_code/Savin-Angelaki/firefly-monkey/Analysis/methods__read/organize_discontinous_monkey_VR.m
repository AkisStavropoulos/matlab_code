function [outdata] = organize_discontinous_monkey_VR(data)

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


% perturbations
try
outdata.ptbStart =          data(:, 37); % start of perturbation
outdata.ptbMu =             data(:, 38); % perturbation temporal center
outdata.ptbSigma =          data(:, 39); % perturbation width (not duration)
outdata.ptbLinVel =         data(:, 40)/scale; % max ptb linear velocity
outdata.ptbAngVel =         data(:, 50); % max ptb angular velocity
outdata.ptb =               data(:, 54); % ??? or 55?? logical for perturbation or not in the trial
catch

end

