function [av, aw, v, w, t] = plan_traj(x_tar, y_tar, vmax, dt, nr)
%% Plan a circular trajectory towards a given target

    R = (x_tar^2 + y_tar^2) / (2*abs(x_tar));
    r = sqrt(x_tar^2 + y_tar^2);
    
    theta = atan2d(abs(y_tar), abs(x_tar));
    phi = atan2d(r*sind(theta), R - r*cosd(theta));
    phi = phi*pi/180; % transform to radians
    
    t_d = (R*phi) / vmax;
    v = vmax*ones(1,(round(t_d/dt)));
    v(1:nr) = linspace(0, vmax, nr);
    v = [v  linspace(vmax, 0, nr)];      % m per tau
    w = sign(x_tar) * (v/R);
    w = w*180/pi;               % transform to angles from radians
    av = 10*diff(v)/dt;                        % m per tau^2 scaled by 10 so network output need not be too small
    aw = 10*diff(w)/dt;
    t = linspace(0, (dt*(length(v))), dt) + dt;
    t = t(1:end-1);
    
    

