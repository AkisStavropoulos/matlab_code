function [MSE_sum,pred_sub] = FixedTauModelResErr(tau_est,JS,start_pos,exp_sub)

%% calculate residual error based on Tau estimate

dt = 1/60;
g = 2.1; % gain on wmax
x = 400;
T = 8.5;
ang = 38;

PenaliseAngle = 1;
PenaliseDistance = 1;

for j = 1:length(JS.x)
    
    % simulate velocity profile
    a = exp(-dt/tau_est);
    vmax = findvmax(x,T,tau_est);
    wmax = findvmax(ang,T,tau_est);

    Tsim = length(JS.x{j});
    v = zeros(1,Tsim);  w = zeros(1,Tsim);
    for t = 2:Tsim
        beta = (1 - a);
        v(t) = a*v(t-1) + vmax*beta*JS.y{j}(t);
        
        w(t) = a*w(t-1) + g*wmax*beta*JS.x{j}(t);
    end
    
    % construct trajectory
    xmp = zeros(1,Tsim);
    ymp = zeros(1,Tsim);
    phi = zeros(1,Tsim);
    xmp(1) = start_pos.x(j);
    ymp(1) = start_pos.y(j);
    for t=1:length(v)
        v_x = v(t).*sin(phi(t)); % input rads
        v_y = v(t).*cos(phi(t));
        xmp(t+1) = xmp(t) + v_x*dt;
        ymp(t+1) = ymp(t) + v_y*dt;
        phi(t+1) = phi(t) + (w(t)*pi/180)*dt; % converts degrees to rads
        phi(t+1) = (phi(t+1)>-pi & phi(t+1)<=pi).*phi(t+1) + ...
            (phi(t+1)>pi).*(phi(t+1) - 2*pi) + (phi(t+1)<=-pi).*(phi(t+1) + 2*pi);
    end
    % get subject polar coordiantes at end of trial
    pred_sub.r(j) = sqrt((xmp(end)-xmp(1)).^2 + (ymp(end)-ymp(1)).^2);
    pred_sub.th(j) = atan2d((xmp(end) - xmp(1)),(ymp(end) - ymp(1)));
    % get subject cartesian coordinates at end of trial
    pred_sub.x(j) = xmp(end);
    pred_sub.y(j) = ymp(end);
end

% compute MSE
if 0
    if ~isempty(exp_sub)
        MSE_r = median((pred_sub.r - exp_sub.r).^2); % Median or 95% of distribution of squared error removes outliers
        MSE_th = median((pred_sub.th - exp_sub.th).^2);
        
        MSE_sum = PenaliseDistance*MSE_r + PenaliseAngle*MSE_th;
    else
        MSE_sum = [];
    end
else
     if ~isempty(exp_sub)
         exp_x = exp_sub.r.*sind(exp_sub.th) + start_pos.x;        
         exp_y = exp_sub.r.*cosd(exp_sub.th) + start_pos.y;
         MSE_sum = median(sqrt(PenaliseAngle*(pred_sub.x - exp_x).^2 + PenaliseDistance*(pred_sub.y - exp_y).^2));
     else
         MSE_sum = [];
     end
end