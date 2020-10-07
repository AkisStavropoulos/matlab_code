function [err,x_sub,y_sub,v,w,v_final,w_final] = gen_sim_traj(x_tar,y_tar,tau,vmax,wmax,dt,trltime,sw,w_gain)
%% Simulate a circular trajectory from the subject to the target
Tsim = ceil(trltime./dt);
Tswitch = ceil(sw./dt);

js_x = [ones(Tswitch,1) ; -ones(Tsim-Tswitch,1)].*w_gain;
js_y = [ones(Tswitch,1) ; -ones(Tsim-Tswitch,1)];

a = exp(-dt./tau);

% zero-pad the end of the trial
Tpad = ceil(2/dt);
js_x = [js_x ; zeros(Tpad,1)];
js_y = [js_y ; zeros(Tpad,1)];

Ttrl = Tsim+Tpad;

v = zeros(1,Ttrl);  w = zeros(1,Ttrl);  phi = zeros(1,Ttrl);
for t = 2:Ttrl
    beta = (1 - a);
    v(t) = a*v(t-1) + vmax*beta*js_y(t);
    w(t) = a*w(t-1) + wmax*beta*js_x(t);
end

[x_sub, y_sub, ~, ~, ~, phi] = gen_traj(w, v, 0, 0, dt);

err = sqrt( (x_tar-x_sub(end)).^2 + (y_tar-y_sub(end)).^2 );

v_final = v(end);
w_final = w(end);


% figure;plot((1:Ttrl)*dt,js_y,'linewidth',3); axis([-1 Ttrl*dt+1 -1.5 1.5]); hline(0,'k-')