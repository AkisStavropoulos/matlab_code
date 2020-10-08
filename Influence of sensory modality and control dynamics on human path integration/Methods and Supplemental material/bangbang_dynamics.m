%% example bang-bang control dynamics

x = 400;
T = 8.5;
tau = 2;
s = switchtime(tau,T);
Tpulse = s; 
Tbrake = T-s; 
Tsim = 12;
[j_x, vel, ts, vmax, dist] = tau2vmax(tau, Tpulse, Tbrake, Tsim, x, T, 0);
[~, vel2] = tau2vmax(tau, Tpulse, Tsim-Tpulse, Tsim, x, T, 0);

figure;

subplot(3,1,1);plot(ts,[0 j_x(1:end-1)],'b','linewidth',1);
hline(0,'k');ylabel('control');xlabel('time');

subplot(3,1,2);plot(ts,vel,'b','linewidth',1);
hline([0 vmax -vmax],'k');vline([s T],'k');ylabel('velocity');xlabel('time');
hold on; plot(ts,vel2,'b:','linewidth',1);

subplot(3,1,3);plot(ts,cumsum(vel)*dt,'r','linewidth',1);
hline(dist,'k');vline([s T],'k');ylabel('position');xlabel('time');

suptitle('bang-bang control')

%% Vmax scaling
% vmax scaling as a function of tau/T
tau = 2;
x = 400;
T = linspace(1,500,1000); % T = 8.5;

v_mu = x./T;
s = [];
for i = 1:length(T)
    s(i) = switchtime(tau,T(i));
end
Tpulse = s; 
Tbrake = T-s; 
Tsim = 600;

vmax = [];
for i = 1:length(T)
    [~, ~, ~, vmax(i)] = tau2vmax(tau, Tpulse(i), Tbrake(i), Tsim, x, T(i), 0);
end

figure; hold on;
plot(tau./T,vmax./v_mu);
axis([0 1.5 0 6]); xlabel('\tau/T'); ylabel('v_{max}/v_{\mu}');

% vmax for tau much larger than T
plot(tau./T,4*tau./T,'g')

% Vmax for tau much smaller than T (velocity control)
tau = .1;
x = 400;
T = 8.5;

v_mu = x./T;
s = switchtime(tau,T);
Tpulse = s; 
Tbrake = T-s; 
Tsim = 600;

[~, ~, ~, vmax] = tau2vmax(tau, Tpulse, Tbrake, Tsim, x, T, 0);

hline(vmax./v_mu,'k--')
