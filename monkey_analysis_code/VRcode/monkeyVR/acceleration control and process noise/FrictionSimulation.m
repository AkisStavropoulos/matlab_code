%% Friction simulation

%% Generate velocity profile
dt = 1/60;
taus = [0.3 1 1.5 2];
vmax = 1;
u = [ones(1,100) zeros(1,400)];
ts = dt*(1:numel(u));

clf;
for tau = taus
    
    a = exp(-dt/tau);
    b = 1-a;
    
    v = zeros(1,numel(ts));
    
    for t = 1:numel(ts)-1
        v(t+1) = a*v(t) + b*vmax*u(t);
    end
    acc = [0 diff(v)/dt];
    
    subplot(numel(taus),1,find(tau==taus)); hold on;
    plot(ts,u,'k'); plot(ts,v,'b'); xlabel('time [s]'); ylabel('velocity');
    title(sprintf('tau = %.1f',tau));
    
    % pick threshold
    thresh = 0.1; % 20% of max velocity
    ind = find(v < thresh*vmax & u==0,1);
    
    % add friction
    tt = ts(ind:end);
    
    fr = 0.2;
    for t = ind:numel(ts)-1
        v(t+1) = sqrt(v(t)^2 - 2*fr*v(t)*dt);
    end
    plot(tt,v(ind:end),'r'); hline(thresh)
end