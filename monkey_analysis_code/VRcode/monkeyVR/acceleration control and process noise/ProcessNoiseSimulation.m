%% ProcessNoiseSimulation
rng(0)
plt = 0;
%% generate target positions
N = 1000;

% target position
r_tar = 100 + 300*rand(1,N);
th_tar = 60*rand(1,N) - 30;

% Control Dynamics Parameters
x = 250;
theta = 60;
T = 2;
tau = 3;

% Process Noise parameters
noiseTau = 0.25;
linGain = 1;
angGain = 1;

for n = 1:N

[x_tar,y_tar] = polar2cartY(r_tar(n),th_tar(n));

[finalPosX, finalPosY, finalPosPhi, cleanPosX, cleanPosY, cleanPosPhi, finalLinVel, finalAngVel, cleanLinVel, cleanAngVel, linEta, angEta, ProcessNoiseMag, jsx, jsy, ts] = ...
    simulate_trajectory(x_tar,y_tar,x,theta,T,tau,noiseTau,linGain,angGain,plt);

[r_sub(n),th_sub(n)] = cart2polarY(finalPosX(end),finalPosY(end));

if plt
subplot(2,4,1);hold on; plot(x_tar,y_tar,'go','markersize',7,'markerfacecolor','g')
suptitle(['tau=' num2str(tau) ', noiseTau=' num2str(noiseTau) ', G_{lin}=' num2str(linGain) ', G_{ang}=' num2str(angGain)]);
end

if numel(tau) > 1
vmax = findvmax(x,T,tau);
noisevar(k) = std((finalLinVel - cleanLinVel)/vmax);
PNmag(k) = max(ProcessNoiseMag);

subplot(2,4,4); hold on;
plot(taus,noisevar);xlabel('tau [s]');ylabel('noise variability');title('noise variability dependence on tau');

subplot(2,4,8); hold on;
plot(taus,PNmag);xlabel('tau [s]');ylabel('noise magnitude');title('noise magnitude dependence on tau');
end

end


%% plot responses vs targets

figure; subplot(1,2,1); hold on;
plot(r_tar,r_sub,'r.'); plot(0:600,0:600,'k--'); axis equal; axis([0 600 0 600]);
xlabel('target distance [cm]'); ylabel('response distance [cm]'); title('radial');

subplot(1,2,2); hold on;
plot(th_tar,th_sub,'r.'); plot(-50:50,-50:50,'k--'); axis equal; axis([-50 50 -50 50]); vline(0,'k'); hline(0,'k');
xlabel('target angle [cm]'); ylabel('response angle [cm]'); title('angular');

suptitle(['G_{lin}='  num2str(linGain) ', G_{ang}=' num2str(angGain) ])

r_sd = std(r_sub - r_tar);
th_sd = std(th_sub - th_tar);
disp(['...... SD of uncompensated process noise (G_lin,G_ang) = (' num2str(linGain) ',' num2str(angGain) '):'])
disp(['radial = ' num2str(r_sd)])
disp(['angular = ' num2str(th_sd)])

