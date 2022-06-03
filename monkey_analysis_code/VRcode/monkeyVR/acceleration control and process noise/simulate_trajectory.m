function [finalPosX, finalPosY, finalPosPhi, cleanPosX, cleanPosY, cleanPosPhi, finalLinVel, finalAngVel, cleanLinVel, cleanAngVel, linEta, angEta, ProcessNoiseMag, jsx, jsy, ts] = ...
    simulate_trajectory(x_tar,y_tar,x,theta,T,tau,noiseTau,linGain,angGain,plt)

%% simulate trajectories under Control Dynamics and Process Noise

dt = 1/90;

% plan trajectory for given parameters
[av, aw, v, w, jsy, jsx, t] = plan_traj_tau(x_tar, y_tar, x, theta, T, tau, dt); % [a, b, v, w, t] = plan_traj(x_tar, y_tar, vmax, dt, nr);

% apply process noise
[finalLinVel,finalAngVel,cleanLinVel,cleanAngVel,linEta,angEta,linKsi,angKsi,ProcessNoiseMag,ts] = ...
    ProcessNoise(jsx,jsy,x,theta,T,tau,noiseTau,linGain,angGain);

% calculate displacement
[cleanPosX, cleanPosY, ~, ~, ~, cleanPosPhi] = gen_traj(cleanAngVel, cleanLinVel, ts);

[finalPosX, finalPosY, ~, ~, ~, finalPosPhi] = gen_traj(finalAngVel, finalLinVel, ts);

vmax = findvmax(x,T,tau);
wmax = findvmax(theta,T,tau);


%% Plots
if plt
figure(999);clf;
subplot(2,4,1); hold on;
plot(cleanPosX,cleanPosY,'k'); plot(finalPosX,finalPosY,'r');
title('trajectory'); %axis([-1000 1000 -1000 1000]);

subplot(2,4,2); hold on;
plot(ts,cleanLinVel,'k'); plot(ts,finalLinVel,'r');
title('linear velocity'); xlabel('time [s]');

subplot(2,4,3); hold on;
plot(ts,cleanAngVel,'k'); plot(ts,finalAngVel,'r');
title('angular velocity'); xlabel('time [s]');

subplot(2,4,5); hold on;
plot(ts,linKsi,'b'); plot(ts,linEta,'m');
title('linear noise'); xlabel('time [s]'); ylim([-0.5 0.5]);

subplot(2,4,6); hold on;
plot(ts,angKsi,'b'); plot(ts,angEta,'m');
title('angular noise'); xlabel('time [s]'); ylim([-0.5 0.5]);

subplot(2,4,7); hold on;
plot(ts,(finalLinVel - cleanLinVel)./cleanLinVel,'m'); plot(ts,(finalAngVel - cleanAngVel)./cleanAngVel,'b');
title({'effective noise','(normalized)'}); xlabel('time [s]'); ylim([-0.5 0.5]); legend({'linear','angular'})
end
