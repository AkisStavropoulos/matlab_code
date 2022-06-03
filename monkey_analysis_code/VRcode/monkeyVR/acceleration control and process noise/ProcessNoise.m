function [finalLinVel,finalAngVel,cleanLinVel,cleanAngVel,linEta,angEta,linKsi,angKsi,ProcessNoiseMag,ts] = ...
    ProcessNoise(JSx,JSy,x,theta,T,tau,noiseTau,linGain,angGain)
%% ideal parameters
% N = 500;
% JSx = [ones(1,N/2)/40 zeros(1,N/2)]; % zeros(1,N);
% JSy = [ones(1,N/2)/2 zeros(1,N/2)];
% x = 250;
% theta = 80;
% T = 1.25;
% tau = 3;
% noiseTau = 0.25;
% linGain = 1;
% angGain = 1;

%% Process Noise design

dt = 1/90;
N = numel(JSx);

vmax = findvmax(x,T,tau);
wmax = findvmax(theta,T,tau);
alpha = exp(-dt/tau);
beta = 1-alpha;

gamma = exp(-dt/noiseTau);
delta = 1-gamma;

linKsi = zeros(1,N); linEta = zeros(1,N); angKsi = zeros(1,N); angEta = zeros(1,N);
cleanLinVel = zeros(1,N); cleanAngVel = zeros(1,N);
for t = 2:N
    % linear component
    linKsi(t) = gamma*linKsi(t-1) + delta*linGain*randn;
    linEta(t) = gamma*linEta(t-1) + delta*linKsi(t);
    cleanLinVel(t) = alpha*cleanLinVel(t-1) + beta*vmax*JSy(t);
    
    % angular component
    angKsi(t) = gamma*angKsi(t-1) + delta*angGain*randn;
    angEta(t) = gamma*angEta(t-1) + delta*angKsi(t);
    cleanAngVel(t) = alpha*cleanAngVel(t-1) + beta*wmax*JSx(t);
    
    ProcessNoiseMag(t) = sqrt((cleanLinVel(t)/vmax)^2 + (cleanAngVel(t)/wmax)^2);
    
    if abs(linEta(t)*ProcessNoiseMag(t)) > vmax
        finalLinVel(t) = cleanLinVel(t) + sign(linEta(t))*cleanLinVel(t);
    else
        finalLinVel(t) = cleanLinVel(t) + linEta(t)*ProcessNoiseMag(t)*vmax;
    end
    
    if abs(linEta(t)*ProcessNoiseMag(t)) > wmax
        finalAngVel(t) = cleanAngVel(t) + sign(angEta(t))*cleanAngVel(t);
    else
        finalAngVel(t) = cleanAngVel(t) + angEta(t)*ProcessNoiseMag(t)*wmax;
    end
end

% calculate displacement
ts = (1:N)*dt;
[cleanPosX, cleanPosY, ~, ~, ~, cleanPosPhi] = gen_traj(cleanAngVel, cleanLinVel, ts);

[finalPosX, finalPosY, ~, ~, ~, finalPosPhi] = gen_traj(finalAngVel, finalLinVel, ts);

