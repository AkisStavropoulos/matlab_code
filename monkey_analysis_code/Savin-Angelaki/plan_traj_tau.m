function [av, aw, v, w, jsy, jsx, t] = plan_traj_tau(x_tar, y_tar, x, th, T, tau, dt)
%% Plan a circular trajectory towards a given target
% under Control Dynamics

R = (x_tar^2 + y_tar^2) / (2*abs(x_tar));
r = sqrt(x_tar^2 + y_tar^2);

theta = atan2d(abs(y_tar), abs(x_tar));
phi = atan2d(r*sind(theta), R - r*cosd(theta));
phi = phi*pi/180; % transform to radians
traj_length = phi*R;
dist = sqrt(x_tar^2 + y_tar^2);

vmax = findvmax(x,T,tau);
wmax = findvmax(th,T,tau);
alpha = exp(-dt/tau);
beta = 1-alpha;
traveltime = bangbangdur(traj_length,tau,vmax);
sw = switchtime(tau, traveltime);

t_d = (R*phi) / vmax;
t_d = sw;

St = round(t_d/dt);
Nt = round(traveltime/dt);
v = zeros(1,St);
% acceleration
for t = 2:St
    v(t) = alpha*v(t-1) + beta*vmax;
end
% braking
for t = St+1:Nt
    v(t) = alpha*v(t-1) - beta*vmax;
end
v(v <= 0) = 0;

w = sign(x_tar) * (v/R);
w = w*180/pi;               % transform to angles from radians
av = 10*diff(v)/dt;                        % m per tau^2 scaled by 10 so network output need not be too small
aw = 10*diff(w)/dt;
t = linspace(0, (dt*(length(v))), dt) + dt;
t = t(1:end-1);

% generate joystick input
jsy = zeros(1,Nt); jsx = zeros(1,Nt);
for t = 2:Nt
    jsx(t) = (w(t)-alpha*w(t-1)) / (beta*wmax);
    jsy(t) = (v(t)-alpha*v(t-1)) / (beta*vmax);
end


