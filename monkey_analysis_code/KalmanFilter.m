%% KalmanFilter
% t0: given state
% t1: predicted next state
% t1_f: state after filter
%% given state
X = -100:0.1:100;
m_p = 10;
sig_p = 1;
pos = normpdf(X,m_p,sig_p);
%
m_v = 2;
sig_v = 2;
vel = normpdf(X,m_v,sig_v);
%
Xt0 = [m_p ; m_v] % given state
Sig_t0 = [sig_p*sig_p  0; ...
          0  sig_v*sig_v] % sig_p*sig_v are 0 because they're independent
%% LOOP
sensor_p = [];
sensor_v = [];
curr_X = [];
pred_X = [];
final_X = [];
measurements = [];
Error = [];
for t = 1:10
%%  predicted state
dt = 1; % 1 sec
F = [1  dt;...
    0  1]; % Prediction matrix
Xt1 = F*Xt0;
Sig_t1 = F*Sig_t0*F';
% 
pos1 = normpdf(X,Xt1(1),sqrt(Sig_t1(1,1)));
vel1 = normpdf(X,Xt1(2),sqrt(Sig_t1(2,2)));
figure;plot(X,pos1);title('pred. state position before/after ext. input & pred. measurement');hold on; % blue
%% External influence
% acceleration 
a = [1 0;0 1]; %    1 m/s^2
B = a*[dt^2/2 ; dt];
Q = [(t-dt)/2 0 ; 0 1]; % noise of external influence, affects only velocity in first step
Xt1 = Xt1 + B;
Sig_t1 = Sig_t1 + Q;
pos1 = normpdf(X,Xt1(1),sqrt(Sig_t1(1,1)));
vel1 = normpdf(X,Xt1(2),sqrt(Sig_t1(2,2)));
plot(X,pos1); % orange
%% predicted measurement
H = [1 0;...
     0 1]; % sensor expectation matrix
Zexp = H*Xt1;
Sig_exp = H*Sig_t1*H';
% 
pos_exp = normpdf(X,Zexp(1),sqrt(Sig_exp(1,1)));
vel_exp = normpdf(X,Zexp(2),sqrt(Sig_exp(2,2)));
plot(X,pos_exp,'g'); % green
%% real measurement
meas_noise = 0.8 + rand*0.4; % random noise between 0.8 and 1.2
new_z_p = meas_noise*Xt1(1);
new_z_v = meas_noise*Xt1(2);
sensor_p = [sensor_p new_z_p];
sensor_v = [sensor_v new_z_v];
z_p = sensor_p(t); % means of p and v (actual measurements)
z_v = sensor_v(t);
z = [z_p; z_v];
Sig_z = [2*(t-dt) 0 ;0 2*(t-dt)]; % sensor's measurement noise matrix, [2*(t-dt)]
pos_meas = normpdf(X,z(1),sqrt(Sig_z(1,1)));
vel_meas = normpdf(X,z(2),sqrt(Sig_z(2,2)));
plot(X,pos_meas, 'k'); % black
%% The filter
K = (Sig_t1*H')*inv(Sig_exp + Sig_z)
E = z - Zexp; % ERROR
Xt1_f = Xt1 + K*E
Sig_t1_f = Sig_t1 - K*H*Sig_t1
pos_f = normpdf(X,Xt1_f(1),sqrt(Sig_t1_f(1,1)));
vel_f = normpdf(X,Xt1_f(2),sqrt(Sig_t1_f(2,2)));
plot(X,pos_f, 'r');hold off; % red
%% record keeping
curr_X = [curr_X Xt0];
pred_X = [pred_X Zexp];
measurements = [measurements z];
Error = [Error E];
%% update
Xt0 = Xt1_f;
Sig_t0 = Sig_t1_f;
if t == 10
    curr_X = [curr_X Xt0];
end
end