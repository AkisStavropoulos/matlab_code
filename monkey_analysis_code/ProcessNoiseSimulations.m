%% Process noise simulations

%% On velocity

%% Multiplicative noise
dt = 1/60;
maxlength = 1000;
Nt = 1000; % Ntrials

% specify control parameters
tau = unique([0.5*ones(1,Nt/2) 2*ones(1,Nt/2) 5]); % exp(subject(1).trials(1).prs.mu_phi+subject(1).trials(1).prs.sig_phi*randn(N,1));
x = 400;
T = 8.5;
vmax = (x/T)*( 1 ./ (-1 + 2*(tau./T) .* log((1 + exp(T./tau))/2))) ;
Ntrials = numel(tau);

traveltime = T*ones(1,Ntrials);
s = arrayfun(@(tau,T) switchtime(tau, T), tau, traveltime);
S = ceil(s/dt);
Ts = ceil(traveltime/dt);

JSinput = arrayfun(@(s,t) [ones(1,s) -ones(1,t-s) zeros(1,maxlength-t)], S, Ts,'un',0);

% specify process noise parameters
gain_eta = 0.1*ones(1,Ntrials); % randsample([0.01 0.05 0.1],Nt,1); % 0.1;
tau_eta = 2*ones(1,Ntrials); % randsample([0.5 2],Nt,1); % 2;

% generate noisy velocity profiles
for j = 1:Ntrials
rng(1);
    
    v{j} = zeros(1,maxlength);
    v_clean{j} = zeros(1,maxlength);
    v_est{j} = zeros(1,maxlength);
    eta{j} = zeros(1,maxlength);
    eta_est{j} = zeros(1,maxlength);
    eta_final{j} = zeros(1,maxlength);
    alpha = exp(-dt./tau(j));
    beta = (1 - alpha);
    gamma = exp(-dt./tau_eta(j));
    delta = (1 - gamma);

    for t = 2:maxlength
        eta{j}(t) = gamma*eta{j}(t-1) + delta*(gain_eta(j)*v{j}(t-1)*randn);
        v{j}(t) = alpha*v{j}(t-1) + vmax(j)*beta*JSinput{j}(t) + eta{j}(t);
        
        v_clean{j}(t) = alpha*v_clean{j}(t-1) + vmax(j)*beta*JSinput{j}(t);
        v_est{j}(t) = alpha*v{j}(t-1) + vmax(j)*beta*JSinput{j}(t);

%         eta_final{j}(t) = alpha*sum(eta{j}(1:t-1)) + beta*eta{j}(t);

    end
    
    N = 600;
    p0 = [alpha];
    [p,fval(j)] = fminsearch(@(p0) mean(((v{j}(2:N)-v_clean{j}(2:N)) - (exp(-dt./p0)*cumsum(eta{j}(1:N-1)) + (1-exp(-dt./p0))*eta{j}(2:N))/gain_eta(j)).^2),p0);
    
    eta_final{j}(2:end) = (exp(-dt./p)*cumsum(eta{j}(1:end-1)) + (1-exp(-dt./p))*eta{j}(2:end))/gain_eta(j);
    
    w(j) = p;
    
    eta_est{j} = v{j} - v_est{j};
    
    d_true(j) = sum(v{j})*dt;
    d_blv(j) = sum(v_clean{j})*dt;
end


ts = (1:maxlength)*dt;

% Plot
if Ntrials < 4
    
aaa = figure;
for j = 1:Ntrials
    subplot(Ntrials,1,j); hold on;
    plot(ts,JSinput{j}*vmax(j),'color',[.7 .7 .7]);
    plot(ts,v_clean{j},'k');
    plot(ts,v{j},'b');
    plot(ts,eta_final{j},'r'); % plot(ts,10*(eta{j}),'r');
    plot(ts,cumsum(eta{j}),'r:');
%     plot(ts,v{j}-cumsum(eta{j}),'k--')
    plot(ts,v{j}-v_clean{j},'k:')% plot(ts,10*eta_est{j},'k:'); % plot(ts,v{j}-v_clean{j},'k:')
    title(['\tau_{trial} = ' num2str(tau(j)) ' s']); xlabel('time [s]'); ylabel('velocity [cm/s]'); ylim([-130 130]);
end
suptitle(['MULTIPLICATIVE: G_{\eta} = ' num2str(unique(gain_eta)) ', \tau_{\eta} = ' num2str(unique(tau_eta)) ' s']);
legend({'joystick input','clean velocity','noisy velocity','noise input', 'noise cumsum','v_{true}- v_{clean}'})
aaa.Position = [1458 136 324 813];



else
    
    taugroup = unique(tau);
    tauindx = tau == taugroup(1); % indexes for small taus
    
    colr = brewermap(numel(taugroup),'Set1');
    figure;
    subplot(2,2,1); hold on;
    gains = unique(gain_eta);
    arrayfun(@(i) errorbar(gains(i),mean(d_true(gain_eta == gains(i))),std(d_true(gain_eta == gains(i))),...
        'bo','capsize',0,'markerfacecolor','b'), 1:numel(gains)); % plot(gain_eta,d_true,'b.'); 
    axis([-0.05 0.15 100 700]); hline(x,'k--'); xlabel('noise gain'); ylabel('stopping distance [cm]');
    
    subplot(2,2,2); hold on;
    noisetaus = unique(tau_eta);
    arrayfun(@(i) errorbar(noisetaus(i),mean(d_true(tau_eta == noisetaus(i))),std(d_true(tau_eta == noisetaus(i))),...
        'bo','capsize',0,'markerfacecolor','b'), 1:numel(noisetaus)); % plot(tau_eta,d_true,'b.'); 
    axis([0 2.5 200 600]); hline(x,'k--'); xlabel('noise time constant'); ylabel('stopping distance [cm]');
    
    subplot(2,2,3); hold on;
    arrayfun(@(i) errorbar(gains(i)-0.005,mean(d_true(gain_eta == gains(i) & tauindx)),std(d_true(gain_eta == gains(i) & tauindx)),...
        'o','color',colr(1,:),'capsize',0,'markerfacecolor',colr(1,:),'handlevisibility','off'), 1:numel(gains)); 
    arrayfun(@(i) errorbar(gains(i)+0.005,mean(d_true(gain_eta == gains(i) & ~tauindx)),std(d_true(gain_eta == gains(i) & ~tauindx)),...
        'o','color',colr(2,:),'capsize',0,'markerfacecolor',colr(2,:),'handlevisibility','off'), 1:numel(gains));
    plot(gains(1)-0.005,mean(d_true(gain_eta == gains(1) & tauindx)),'o','color',colr(1,:),'markerfacecolor',colr(1,:))
    plot(gains(1)+0.005,mean(d_true(gain_eta == gains(1) & ~tauindx)),'o','color',colr(2,:),'markerfacecolor',colr(2,:))
%     plot(gain_eta(tauindx),d_true(tauindx),'.','color',colr(1,:)); plot(gain_eta(~tauindx),d_true(~tauindx),'.','color',colr(2,:)); 
    axis([-0.05 0.15 100 700]); hline(x,'k--'); xlabel('noise gain'); ylabel('stopping distance [cm]'); legend({['tau = ' num2str(taugroup(1))],['tau = ' num2str(taugroup(2))]});
    
    subplot(2,2,4); hold on;
    arrayfun(@(i) errorbar(noisetaus(i)-0.1,mean(d_true(tau_eta == noisetaus(i) & tauindx)),std(d_true(tau_eta == noisetaus(i) & tauindx)),...
        'o','color',colr(1,:),'capsize',0,'markerfacecolor',colr(1,:),'handlevisibility','off'), 1:numel(noisetaus)); 
    arrayfun(@(i) errorbar(noisetaus(i)+0.1,mean(d_true(tau_eta == noisetaus(i) & ~tauindx)),std(d_true(tau_eta == noisetaus(i) & ~tauindx)),...
        'o','color',colr(2,:),'capsize',0,'markerfacecolor',colr(2,:),'handlevisibility','off'), 1:numel(noisetaus)); 
    plot(noisetaus(1)-0.1,mean(d_true(tau_eta == noisetaus(1) & tauindx)),'o','color',colr(1,:),'markerfacecolor',colr(1,:))
    plot(noisetaus(1)+0.1,mean(d_true(tau_eta == noisetaus(1) & ~tauindx)),'o','color',colr(2,:),'markerfacecolor',colr(2,:))
%     plot(tau_eta(tauindx),d_true(tauindx),'.','color',colr(1,:)); plot(tau_eta(~tauindx),d_true(~tauindx),'.','color',colr(2,:)); 
    axis([0 2.5 200 600]); hline(x,'k--'); xlabel('noise time constant'); ylabel('stopping distance [cm]'); legend({['tau = ' num2str(taugroup(1))],['tau = ' num2str(taugroup(2))]});
    
end

%% Additive noise

% specify control parameters
velthresh = 5;
Nt = 500;
dt = 1/60;
maxlength = 1000;
x = 400;
T = 8.5;
tau = [0.2 5]; % exp(subject(1).trials(1).prs.mu_phi+subject(1).trials(1).prs.sig_phi*randn(N,1));
vmax = (x/T)*( 1 ./ (-1 + 2*(tau./T) .* log((1 + exp(T./tau))/2))) ;
Ntrials = numel(tau);

traveltime = T*ones(1,Ntrials);
s = arrayfun(@(tau,T) switchtime(tau, T), tau, traveltime);
S = ceil(s/dt);
Ts = ceil(traveltime/dt);

JSinput = arrayfun(@(s,t) [ones(1,s) -ones(1,t-s) zeros(1,maxlength-t)], S, Ts,'un',0);

% specify process noise parameters
gain_eta = 0.1;
tau_eta = 2;
gamma = exp(-dt./tau_eta);
delta = (1 - gamma);

% generate noisy velocity profiles
for j = 1:Ntrials
 rng(1);
    
    v{j} = zeros(1,maxlength);
    v_clean{j} = zeros(1,maxlength);
    eta{j} = zeros(1,maxlength);
    alpha = exp(-dt./tau(j));
    beta = (1 - alpha);
    
    for t = 2:maxlength
        if abs(v{j}(t-1)) > velthresh;   eta{j}(t) = gamma*eta{j}(t-1) + delta*(gain_eta*vmax(j)*randn);
        else;                          eta{j}(t) = gamma*eta{j}(t-1) + delta*(gain_eta*vmax(j)*0);
        end
        v{j}(t) = alpha*v{j}(t-1) + vmax(j)*beta*JSinput{j}(t) + eta{j}(t);
        
        v_clean{j}(t) = alpha*v_clean{j}(t-1) + vmax(j)*beta*JSinput{j}(t);
        
    end
    d_true(j) = sum(v{j})*dt;
    d_blv(j) = sum(v_clean{j})*dt;
end

ts = (1:maxlength)*dt;

% Plot
if Ntrials < 4
    
aaa = figure;
for j = 1:Ntrials
    subplot(Ntrials,1,j); hold on;
    plot(ts,JSinput{j}*vmax(j),'color',[.7 .7 .7]);
    plot(ts,v_clean{j},'k');
    plot(ts,v{j},'b');
    plot(ts,10*(eta{j}),'r');
    plot(ts,cumsum(eta{j}),'r:');
%     plot(ts,v{j}-cumsum(eta{j}),'k--')
    plot(ts,v{j}-v_clean{j},'k:')
    title(['\tau_{trial} = ' num2str(tau(j)) ' s']); xlabel('time [s]'); ylabel('velocity [cm/s]');
end
suptitle(['ADDITIVE: G_{\eta} = ' num2str(gain_eta) ', \tau_{\eta} = ' num2str(tau_eta) ' s']);
legend({'joystick input','clean velocity','noisy velocity','noise input', 'noise cumsum','v_{true}- v_{clean}'})
aaa.Position = [1458 136 324 813];

end


%% On control input

%% Multiplicative noise
dt = 1/60;
maxlength = 1000;
Nt = 1000; % Ntrials

% specify control parameters
tau = unique([0.2*ones(1,Nt/2) 5*ones(1,Nt/2)]); % exp(subject(1).trials(1).prs.mu_phi+subject(1).trials(1).prs.sig_phi*randn(N,1));
x = 400;
T = 8.5;
vmax = (x/T)*( 1 ./ (-1 + 2*(tau./T) .* log((1 + exp(T./tau))/2))) ;
Ntrials = numel(tau);

traveltime = T*ones(1,Ntrials);
s = arrayfun(@(tau,T) switchtime(tau, T), tau, traveltime);
S = ceil(s/dt);
Ts = ceil(traveltime/dt);

JSinput = arrayfun(@(s,t) [ones(1,s) -ones(1,t-s) zeros(1,maxlength-t)], S, Ts,'un',0);

% specify process noise parameters
gain_eta = 0.1*ones(1,Ntrials); % randsample([0.01 0.05 0.1],Nt,1); % 0.1;
tau_eta = 2*ones(1,Ntrials); % randsample([0.5 2],Nt,1); % 2;

% generate noisy velocity profiles
for j = 1:Ntrials
rng(1);
    
    v{j} = zeros(1,maxlength);
    v_clean{j} = zeros(1,maxlength);
    v_est{j} = zeros(1,maxlength);
    eta{j} = zeros(1,maxlength);
    eta_est{j} = zeros(1,maxlength);
    alpha = exp(-dt./tau(j));
    beta = (1 - alpha);
    gamma = exp(-dt./tau_eta(j));
    delta = (1 - gamma);

    for t = 2:maxlength
        eta{j}(t) = gamma*eta{j}(t-1) + delta*(gain_eta(j)*vmax(j)*JSinput{j}(t-1)*randn);
        v{j}(t) = alpha*v{j}(t-1) + vmax(j)*beta*JSinput{j}(t) + eta{j}(t);
        
        v_clean{j}(t) = alpha*v_clean{j}(t-1) + vmax(j)*beta*JSinput{j}(t);
        v_est{j}(t) = alpha*v{j}(t-1) + vmax(j)*beta*JSinput{j}(t);

    end
    eta_est{j} = v{j} - v_est{j};
    
    d_true(j) = sum(v{j})*dt;
    d_blv(j) = sum(v_clean{j})*dt;
end

ts = (1:maxlength)*dt;

% Plot
if Ntrials < 4
    
aaa = figure;
for j = 1:Ntrials
    subplot(Ntrials,1,j); hold on;
    plot(ts,JSinput{j}*vmax(j),'color',[.7 .7 .7]);
    plot(ts,v_clean{j},'k');
    plot(ts,v{j},'b');
    plot(ts,10*(eta{j}),'r');
    plot(ts,cumsum(eta{j}),'r:');
%     plot(ts,v{j}-cumsum(eta{j}),'k--')
    plot(ts,v{j}-v_clean{j},'k:')
    title(['\tau_{trial} = ' num2str(tau(j)) ' s']); xlabel('time [s]'); ylabel('velocity [cm/s]');
end
suptitle(['MULTIPLICATIVE in JS: G_{\eta} = ' num2str(unique(gain_eta)) ', \tau_{\eta} = ' num2str(unique(tau_eta)) ' s']);
legend({'joystick input','clean velocity','noisy velocity','noise input', 'noise cumsum','v_{true}- v_{clean}'})
aaa.Position = [1458 136 324 813];

else
    
    taugroup = unique(tau);
    tauindx = tau == taugroup(1); % indexes for small taus
    
    colr = brewermap(numel(taugroup),'Set1');
    figure;
    subplot(2,2,1); hold on;
    gains = unique(gain_eta);
    arrayfun(@(i) errorbar(gains(i),mean(d_true(gain_eta == gains(i))),std(d_true(gain_eta == gains(i))),...
        'bo','capsize',0,'markerfacecolor','b'), 1:numel(gains)); % plot(gain_eta,d_true,'b.'); 
    axis([-0.05 0.15 100 700]); hline(x,'k--'); xlabel('noise gain'); ylabel('stopping distance [cm]');
    
    subplot(2,2,2); hold on;
    noisetaus = unique(tau_eta);
    arrayfun(@(i) errorbar(noisetaus(i),mean(d_true(tau_eta == noisetaus(i))),std(d_true(tau_eta == noisetaus(i))),...
        'bo','capsize',0,'markerfacecolor','b'), 1:numel(noisetaus)); % plot(tau_eta,d_true,'b.'); 
    axis([0 2.5 200 600]); hline(x,'k--'); xlabel('noise time constant'); ylabel('stopping distance [cm]');
    
    subplot(2,2,3); hold on;
    arrayfun(@(i) errorbar(gains(i)-0.005,mean(d_true(gain_eta == gains(i) & tauindx)),std(d_true(gain_eta == gains(i) & tauindx)),...
        'o','color',colr(1,:),'capsize',0,'markerfacecolor',colr(1,:),'handlevisibility','off'), 1:numel(gains)); 
    arrayfun(@(i) errorbar(gains(i)+0.005,mean(d_true(gain_eta == gains(i) & ~tauindx)),std(d_true(gain_eta == gains(i) & ~tauindx)),...
        'o','color',colr(2,:),'capsize',0,'markerfacecolor',colr(2,:),'handlevisibility','off'), 1:numel(gains));
    plot(gains(1)-0.005,mean(d_true(gain_eta == gains(1) & tauindx)),'o','color',colr(1,:),'markerfacecolor',colr(1,:))
    plot(gains(1)+0.005,mean(d_true(gain_eta == gains(1) & ~tauindx)),'o','color',colr(2,:),'markerfacecolor',colr(2,:))
%     plot(gain_eta(tauindx),d_true(tauindx),'.','color',colr(1,:)); plot(gain_eta(~tauindx),d_true(~tauindx),'.','color',colr(2,:)); 
    axis([-0.05 0.15 100 700]); hline(x,'k--'); xlabel('noise gain'); ylabel('stopping distance [cm]'); legend({['tau = ' num2str(taugroup(1))],['tau = ' num2str(taugroup(2))]});
    
    subplot(2,2,4); hold on;
    arrayfun(@(i) errorbar(noisetaus(i)-0.1,mean(d_true(tau_eta == noisetaus(i) & tauindx)),std(d_true(tau_eta == noisetaus(i) & tauindx)),...
        'o','color',colr(1,:),'capsize',0,'markerfacecolor',colr(1,:),'handlevisibility','off'), 1:numel(noisetaus)); 
    arrayfun(@(i) errorbar(noisetaus(i)+0.1,mean(d_true(tau_eta == noisetaus(i) & ~tauindx)),std(d_true(tau_eta == noisetaus(i) & ~tauindx)),...
        'o','color',colr(2,:),'capsize',0,'markerfacecolor',colr(2,:),'handlevisibility','off'), 1:numel(noisetaus)); 
    plot(noisetaus(1)-0.1,mean(d_true(tau_eta == noisetaus(1) & tauindx)),'o','color',colr(1,:),'markerfacecolor',colr(1,:))
    plot(noisetaus(1)+0.1,mean(d_true(tau_eta == noisetaus(1) & ~tauindx)),'o','color',colr(2,:),'markerfacecolor',colr(2,:))
%     plot(tau_eta(tauindx),d_true(tauindx),'.','color',colr(1,:)); plot(tau_eta(~tauindx),d_true(~tauindx),'.','color',colr(2,:)); 
    axis([0 2.5 200 600]); hline(x,'k--'); xlabel('noise time constant'); ylabel('stopping distance [cm]'); legend({['tau = ' num2str(taugroup(1))],['tau = ' num2str(taugroup(2))]});
    
end

%% On position

%% Multiplicative noise on velocity added to position
dt = 1/60;
maxlength = 1000;
Nt = 1000; % Ntrials

% specify control parameters
tau = unique([0.2*ones(1,Nt/2) 5*ones(1,Nt/2)]); % exp(subject(1).trials(1).prs.mu_phi+subject(1).trials(1).prs.sig_phi*randn(N,1));
x = 400;
T = 8.5;
vmax = (x/T)*( 1 ./ (-1 + 2*(tau./T) .* log((1 + exp(T./tau))/2))) ;
Ntrials = numel(tau);

traveltime = T*ones(1,Ntrials);
s = arrayfun(@(tau,T) switchtime(tau, T), tau, traveltime);
S = ceil(s/dt);
Ts = ceil(traveltime/dt);

JSinput = arrayfun(@(s,t) [ones(1,s) -ones(1,t-s) zeros(1,maxlength-t)], S, Ts,'un',0);

% specify process noise parameters
gain_eta = 0.1*ones(1,Ntrials); % randsample([0.01 0.05 0.1],Nt,1); % 0.1;
tau_eta = 2*ones(1,Ntrials); % randsample([0.5 2],Nt,1); % 2;

% generate noisy velocity profiles
for j = 1:Ntrials
rng(1);
    
    v_clean{j} = zeros(1,maxlength);
    v_est{j} = zeros(1,maxlength);
    v_inst{j} = zeros(1,maxlength);
    v{j} = zeros(1,maxlength);
    y{j} = zeros(1,maxlength);
    eta{j} = zeros(1,maxlength);
    ksi{j} = zeros(1,maxlength);
    eta_est{j} = zeros(1,maxlength);
    alpha = exp(-dt./tau(j));
    beta = (1 - alpha);
    gamma = exp(-dt./tau_eta(j));
    delta = (1 - gamma);
    epsilon = exp(-dt./tau_eta(j));
    zeta = (1 - epsilon);

    for t = 2:maxlength
        ksi{j}(t) = epsilon*ksi{j}(t-1) + zeta*randn;
        
        eta{j}(t) = gamma*eta{j}(t-1) + delta*(gain_eta(j)*v{j}(t-1)*ksi{j}(t));
        v_clean{j}(t) = alpha*v_clean{j}(t-1) + vmax(j)*beta*JSinput{j}(t);
        
        y{j}(t) = y{j}(t-1) + v_clean{j}(t)*dt + eta{j}(t);
        
        v{j}(t) = (y{j}(t) - y{j}(t-1))/dt;
        v_est{j}(t) = v_clean{j}(t) + eta{j}(t)/dt;
    end
    eta_est{j} = v{j} - v_clean{j};
    
    d_true(j) = y{j}(end);
    d_blv(j) = sum(v_clean{j})*dt;
end

ts = (1:maxlength)*dt;

% Plot
if Ntrials < 4
    
aaa = figure;
for j = 1:Ntrials
    subplot(Ntrials,1,j); hold on;
    plot(ts,JSinput{j}*vmax(j),'color',[.7 .7 .7]);
    plot(ts,v_clean{j},'k');
    plot(ts,v{j},'b'); plot(ts,v_est{j},'g:','handlevisibility','off');
    plot(ts,eta_est{j},'r'); % plot(ts,10*[0 diff(eta{j})]/dt,'r'); % plot(ts,10*(eta{j}),'r');
    plot(ts,cumsum(eta_est{j})*dt,'r:');
%     plot(ts,v{j}-cumsum(eta{j}),'k--')
    plot(ts,v{j}-v_clean{j},'k:')
    title(['\tau_{trial} = ' num2str(tau(j)) ' s']); xlabel('time [s]'); ylabel('velocity [cm/s]');
end
suptitle(['ADDITIVE in POS: G_{\eta} = ' num2str(unique(gain_eta)) ', \tau_{\eta} = ' num2str(unique(tau_eta)) ' s']);
legend({'joystick input','clean velocity','noisy velocity','noise input', 'noise cumsum','v_{true}- v_{clean}'})
aaa.Position = [1458 136 324 813];

else
    
    taugroup = unique(tau);
    tauindx = tau == taugroup(1); % indexes for small taus
    
    colr = brewermap(numel(taugroup),'Set1');
    figure;
    subplot(2,2,1); hold on;
    gains = unique(gain_eta);
    arrayfun(@(i) errorbar(gains(i),mean(d_true(gain_eta == gains(i))),std(d_true(gain_eta == gains(i))),...
        'bo','capsize',0,'markerfacecolor','b'), 1:numel(gains)); % plot(gain_eta,d_true,'b.'); 
    axis([-0.05 0.15 100 700]); hline(x,'k--'); xlabel('noise gain'); ylabel('stopping distance [cm]');
    
    subplot(2,2,2); hold on;
    noisetaus = unique(tau_eta);
    arrayfun(@(i) errorbar(noisetaus(i),mean(d_true(tau_eta == noisetaus(i))),std(d_true(tau_eta == noisetaus(i))),...
        'bo','capsize',0,'markerfacecolor','b'), 1:numel(noisetaus)); % plot(tau_eta,d_true,'b.'); 
    axis([0 2.5 200 600]); hline(x,'k--'); xlabel('noise time constant'); ylabel('stopping distance [cm]');
    
    subplot(2,2,3); hold on;
    arrayfun(@(i) errorbar(gains(i)-0.005,mean(d_true(gain_eta == gains(i) & tauindx)),std(d_true(gain_eta == gains(i) & tauindx)),...
        'o','color',colr(1,:),'capsize',0,'markerfacecolor',colr(1,:),'handlevisibility','off'), 1:numel(gains)); 
    arrayfun(@(i) errorbar(gains(i)+0.005,mean(d_true(gain_eta == gains(i) & ~tauindx)),std(d_true(gain_eta == gains(i) & ~tauindx)),...
        'o','color',colr(2,:),'capsize',0,'markerfacecolor',colr(2,:),'handlevisibility','off'), 1:numel(gains));
    plot(gains(1)-0.005,mean(d_true(gain_eta == gains(1) & tauindx)),'o','color',colr(1,:),'markerfacecolor',colr(1,:))
    plot(gains(1)+0.005,mean(d_true(gain_eta == gains(1) & ~tauindx)),'o','color',colr(2,:),'markerfacecolor',colr(2,:))
%     plot(gain_eta(tauindx),d_true(tauindx),'.','color',colr(1,:)); plot(gain_eta(~tauindx),d_true(~tauindx),'.','color',colr(2,:)); 
    axis([-0.05 0.15 100 700]); hline(x,'k--'); xlabel('noise gain'); ylabel('stopping distance [cm]'); legend({['tau = ' num2str(taugroup(1))],['tau = ' num2str(taugroup(2))]});
    
    subplot(2,2,4); hold on;
    arrayfun(@(i) errorbar(noisetaus(i)-0.1,mean(d_true(tau_eta == noisetaus(i) & tauindx)),std(d_true(tau_eta == noisetaus(i) & tauindx)),...
        'o','color',colr(1,:),'capsize',0,'markerfacecolor',colr(1,:),'handlevisibility','off'), 1:numel(noisetaus)); 
    arrayfun(@(i) errorbar(noisetaus(i)+0.1,mean(d_true(tau_eta == noisetaus(i) & ~tauindx)),std(d_true(tau_eta == noisetaus(i) & ~tauindx)),...
        'o','color',colr(2,:),'capsize',0,'markerfacecolor',colr(2,:),'handlevisibility','off'), 1:numel(noisetaus)); 
    plot(noisetaus(1)-0.1,mean(d_true(tau_eta == noisetaus(1) & tauindx)),'o','color',colr(1,:),'markerfacecolor',colr(1,:))
    plot(noisetaus(1)+0.1,mean(d_true(tau_eta == noisetaus(1) & ~tauindx)),'o','color',colr(2,:),'markerfacecolor',colr(2,:))
%     plot(tau_eta(tauindx),d_true(tauindx),'.','color',colr(1,:)); plot(tau_eta(~tauindx),d_true(~tauindx),'.','color',colr(2,:)); 
    axis([0 2.5 200 600]); hline(x,'k--'); xlabel('noise time constant'); ylabel('stopping distance [cm]'); legend({['tau = ' num2str(taugroup(1))],['tau = ' num2str(taugroup(2))]});
    
end

%% Alternative Additive noise added on position
dt = 1/60;
maxlength = 1000;
Nt = 1000; % Ntrials

% specify control parameters
tau = unique([0.2*ones(1,Nt/2) 2*ones(1,Nt/2) 5]); % exp(subject(1).trials(1).prs.mu_phi+subject(1).trials(1).prs.sig_phi*randn(N,1));
x = 400;
T = 8.5;
vmax = (x/T)*( 1 ./ (-1 + 2*(tau./T) .* log((1 + exp(T./tau))/2))) ;
Ntrials = numel(tau);

traveltime = T*ones(1,Ntrials);
s = arrayfun(@(tau,T) switchtime(tau, T), tau, traveltime);
S = ceil(s/dt);
Ts = ceil(traveltime/dt);

JSinput = arrayfun(@(s,t) [ones(1,s) -ones(1,t-s) zeros(1,maxlength-t)], S, Ts,'un',0);

% specify process noise parameters
gain_eta = 400*ones(1,Ntrials); % randsample([0.01 0.05 0.1],Nt,1); % 0.1;
tau_eta = 2*ones(1,Ntrials); % randsample([0.5 2],Nt,1); % 2;

% generate noisy velocity profiles
for j = 1:Ntrials
rng(1);
    
    y{j} = zeros(1,maxlength);
    v{j} = zeros(1,maxlength);
    v_clean{j} = zeros(1,maxlength);
    v_est{j} = zeros(1,maxlength);
    ksi{j} = zeros(1,maxlength);
    eta{j} = zeros(1,maxlength);
    eta_est{j} = zeros(1,maxlength);
    alpha = exp(-dt./tau(j));
    beta = (1 - alpha);
    gamma = exp(-dt./tau_eta(j));
    delta = (1 - gamma);
    epsilon = exp(-dt./tau_eta(j));
    zeta = (1 - epsilon);

    for t = 2:maxlength
        ksi{j}(t) = epsilon*ksi{j}(t-1) + zeta*randn;
        
        eta{j}(t) = gamma*eta{j}(t-1) + delta*gain_eta(j)*ksi{j}(t);
        v_clean{j}(t) = alpha*v_clean{j}(t-1) + vmax(j)*beta*JSinput{j}(t);
        
        y{j}(t) = y{j}(t-1) + v_clean{j}(t)*dt + eta{j}(t)*dt;
        
        v{j}(t) = (y{j}(t) - y{j}(t-1))/dt;
        v_est{j}(t) = v_clean{j}(t) + eta{j}(t);

    end
        
    eta_est{j} = v{j} - v_est{j};
    
    d_true(j) = sum(v{j})*dt;
    d_blv(j) = sum(v_clean{j})*dt;
end


ts = (1:maxlength)*dt;

% Plot
if Ntrials < 4
    
aaa = figure;
for j = 1:Ntrials
    subplot(Ntrials,1,j); hold on;
    plot(ts,JSinput{j}*vmax(j),'color',[.7 .7 .7]);
    plot(ts,v_clean{j},'k');
    plot(ts,v{j},'b');
    plot(ts,eta{j},'r'); % plot(ts,10*(eta{j}),'r');
%     plot(ts,cumsum(eta{j})*dt,'r:');
%     plot(ts,v{j}-cumsum(eta{j}),'k--')
    plot(ts,v{j}-v_clean{j},'k:')% plot(ts,10*eta_est{j},'k:'); % plot(ts,v{j}-v_clean{j},'k:')
    title(['\tau_{trial} = ' num2str(tau(j)) ' s']); xlabel('time [s]'); ylabel('velocity [cm/s]'); ylim([-130 130]);
end
suptitle(['x_t=x_{t-1}+(v_t+\eta_t)\Deltat :  G_{\eta} = ' num2str(unique(gain_eta)) ', \tau_{\eta} = ' num2str(unique(tau_eta)) ' s']);
legend({'joystick input','clean velocity','noisy velocity','noise input','v_{true}- v_{clean}'})
aaa.Position = [1458 136 324 813];

figure; hold on; cellfun(@(x,xx) plot(x-xx), v, v_clean); xlabel('time [s]'); ylabel('velocity modulation [cm/s]');
title('noise across control dynamics'); legend(arrayfun(@(x) ['tau = ' num2str(x) 's'], tau,'un',0))
suptitle('x_t=x_{t-1}+(v_t+\eta_t)\Deltat')
else
    
    taugroup = unique(tau);
    tauindx = tau == taugroup(1); % indexes for small taus
    
    colr = brewermap(numel(taugroup),'Set1');
    figure;
    subplot(2,2,1); hold on;
    gains = unique(gain_eta);
    arrayfun(@(i) errorbar(gains(i),mean(d_true(gain_eta == gains(i))),std(d_true(gain_eta == gains(i))),...
        'bo','capsize',0,'markerfacecolor','b'), 1:numel(gains)); % plot(gain_eta,d_true,'b.'); 
    axis([-0.05 0.15 100 700]); hline(x,'k--'); xlabel('noise gain'); ylabel('stopping distance [cm]');
    
    subplot(2,2,2); hold on;
    noisetaus = unique(tau_eta);
    arrayfun(@(i) errorbar(noisetaus(i),mean(d_true(tau_eta == noisetaus(i))),std(d_true(tau_eta == noisetaus(i))),...
        'bo','capsize',0,'markerfacecolor','b'), 1:numel(noisetaus)); % plot(tau_eta,d_true,'b.'); 
    axis([0 2.5 200 600]); hline(x,'k--'); xlabel('noise time constant'); ylabel('stopping distance [cm]');
    
    subplot(2,2,3); hold on;
    arrayfun(@(i) errorbar(gains(i)-0.005,mean(d_true(gain_eta == gains(i) & tauindx)),std(d_true(gain_eta == gains(i) & tauindx)),...
        'o','color',colr(1,:),'capsize',0,'markerfacecolor',colr(1,:),'handlevisibility','off'), 1:numel(gains)); 
    arrayfun(@(i) errorbar(gains(i)+0.005,mean(d_true(gain_eta == gains(i) & ~tauindx)),std(d_true(gain_eta == gains(i) & ~tauindx)),...
        'o','color',colr(2,:),'capsize',0,'markerfacecolor',colr(2,:),'handlevisibility','off'), 1:numel(gains));
    plot(gains(1)-0.005,mean(d_true(gain_eta == gains(1) & tauindx)),'o','color',colr(1,:),'markerfacecolor',colr(1,:))
    plot(gains(1)+0.005,mean(d_true(gain_eta == gains(1) & ~tauindx)),'o','color',colr(2,:),'markerfacecolor',colr(2,:))
%     plot(gain_eta(tauindx),d_true(tauindx),'.','color',colr(1,:)); plot(gain_eta(~tauindx),d_true(~tauindx),'.','color',colr(2,:)); 
    axis([-0.05 0.15 100 700]); hline(x,'k--'); xlabel('noise gain'); ylabel('stopping distance [cm]'); legend({['tau = ' num2str(taugroup(1))],['tau = ' num2str(taugroup(2))]});
    
    subplot(2,2,4); hold on;
    arrayfun(@(i) errorbar(noisetaus(i)-0.1,mean(d_true(tau_eta == noisetaus(i) & tauindx)),std(d_true(tau_eta == noisetaus(i) & tauindx)),...
        'o','color',colr(1,:),'capsize',0,'markerfacecolor',colr(1,:),'handlevisibility','off'), 1:numel(noisetaus)); 
    arrayfun(@(i) errorbar(noisetaus(i)+0.1,mean(d_true(tau_eta == noisetaus(i) & ~tauindx)),std(d_true(tau_eta == noisetaus(i) & ~tauindx)),...
        'o','color',colr(2,:),'capsize',0,'markerfacecolor',colr(2,:),'handlevisibility','off'), 1:numel(noisetaus)); 
    plot(noisetaus(1)-0.1,mean(d_true(tau_eta == noisetaus(1) & tauindx)),'o','color',colr(1,:),'markerfacecolor',colr(1,:))
    plot(noisetaus(1)+0.1,mean(d_true(tau_eta == noisetaus(1) & ~tauindx)),'o','color',colr(2,:),'markerfacecolor',colr(2,:))
%     plot(tau_eta(tauindx),d_true(tauindx),'.','color',colr(1,:)); plot(tau_eta(~tauindx),d_true(~tauindx),'.','color',colr(2,:)); 
    axis([0 2.5 200 600]); hline(x,'k--'); xlabel('noise time constant'); ylabel('stopping distance [cm]'); legend({['tau = ' num2str(taugroup(1))],['tau = ' num2str(taugroup(2))]});
    
end

%% Alternative Multiplicative noise added on position
dt = 1/60;
maxlength = 1000;
Nt = 1000; % Ntrials

% specify control parameters
tau = unique([0.2*ones(1,Nt/2) 2*ones(1,Nt/2) 5]); % exp(subject(1).trials(1).prs.mu_phi+subject(1).trials(1).prs.sig_phi*randn(N,1));
x = 400;
T = 8.5;
vmax = (x/T)*( 1 ./ (-1 + 2*(tau./T) .* log((1 + exp(T./tau))/2))) ;
Ntrials = numel(tau);

traveltime = T*ones(1,Ntrials);
s = arrayfun(@(tau,T) switchtime(tau, T), tau, traveltime);
S = ceil(s/dt);
Ts = ceil(traveltime/dt);

JSinput = arrayfun(@(s,t) [ones(1,s) -ones(1,t-s) zeros(1,maxlength-t)], S, Ts,'un',0);

% specify process noise parameters
gain_eta = 10*ones(1,Ntrials); % randsample([0.01 0.05 0.1],Nt,1); % 0.1;
tau_eta = 2*ones(1,Ntrials); % randsample([0.5 2],Nt,1); % 2;

% generate noisy velocity profiles
for j = 1:Ntrials
rng(1);
    
    y{j} = zeros(1,maxlength);
    v{j} = zeros(1,maxlength);
    v_clean{j} = zeros(1,maxlength);
    v_est{j} = zeros(1,maxlength);
    ksi{j} = zeros(1,maxlength);
    eta{j} = zeros(1,maxlength);
    eta_est{j} = zeros(1,maxlength);
    alpha = exp(-dt./tau(j));
    beta = (1 - alpha);
    gamma = exp(-dt./tau_eta(j));
    delta = (1 - gamma);
    epsilon = exp(-dt./tau_eta(j));
    zeta = (1 - epsilon);

    for t = 2:maxlength
        ksi{j}(t) = epsilon*ksi{j}(t-1) + zeta*randn;
        
        eta{j}(t) = gamma*eta{j}(t-1) + delta*gain_eta(j)*ksi{j}(t);
        v_clean{j}(t) = alpha*v_clean{j}(t-1) + vmax(j)*beta*JSinput{j}(t);
        
        y{j}(t) = y{j}(t-1) + (1+eta{j}(t))*v_clean{j}(t)*dt;

        v{j}(t) = (y{j}(t) - y{j}(t-1))/dt;
        v_est{j}(t) = v_clean{j}(t) + eta{j}(t);

    end
        
    eta_est{j} = v{j} - v_est{j};
    
    d_true(j) = sum(v{j})*dt;
    d_blv(j) = sum(v_clean{j})*dt;
end


ts = (1:maxlength)*dt;

% Plot
if Ntrials < 4
    
aaa = figure;
for j = 1:Ntrials
    subplot(Ntrials,1,j); hold on;
    plot(ts,JSinput{j}*vmax(j),'color',[.7 .7 .7]);
    plot(ts,v_clean{j},'k');
    plot(ts,v{j},'b');
    plot(ts,eta{j}.*v_clean{j},'r'); % plot(ts,10*(eta{j}),'r');
%     plot(ts,cumsum(eta{j})*dt,'r:');
%     plot(ts,v{j}-cumsum(eta{j}),'k--')
    plot(ts,v{j}-v_clean{j},'k:')% plot(ts,10*eta_est{j},'k:'); % plot(ts,v{j}-v_clean{j},'k:')
    title(['\tau_{trial} = ' num2str(tau(j)) ' s']); xlabel('time [s]'); ylabel('velocity [cm/s]'); ylim([-130 130]);
end
suptitle(['x_t=x_{t-1}+(1+\eta_t)v_t\Deltat :  G_{\eta} = ' num2str(unique(gain_eta)) ', \tau_{\eta} = ' num2str(unique(tau_eta)) ' s']);
legend({'joystick input','clean velocity','noisy velocity','noise input','v_{true}- v_{clean}'})
aaa.Position = [1458 136 324 813];

figure; hold on; cellfun(@(x,xx) plot(x-xx), v, v_clean); xlabel('time [s]'); ylabel('velocity modulation [cm/s]');
title('noise across control dynamics'); legend(arrayfun(@(x) ['tau = ' num2str(x) 's'], tau,'un',0))
suptitle('x_t=x_{t-1}+(1+\eta_t)v_t\Deltat')

else
    
    taugroup = unique(tau);
    tauindx = tau == taugroup(1); % indexes for small taus
    
    colr = brewermap(numel(taugroup),'Set1');
    figure;
    subplot(2,2,1); hold on;
    gains = unique(gain_eta);
    arrayfun(@(i) errorbar(gains(i),mean(d_true(gain_eta == gains(i))),std(d_true(gain_eta == gains(i))),...
        'bo','capsize',0,'markerfacecolor','b'), 1:numel(gains)); % plot(gain_eta,d_true,'b.'); 
    axis([-0.05 0.15 100 700]); hline(x,'k--'); xlabel('noise gain'); ylabel('stopping distance [cm]');
    
    subplot(2,2,2); hold on;
    noisetaus = unique(tau_eta);
    arrayfun(@(i) errorbar(noisetaus(i),mean(d_true(tau_eta == noisetaus(i))),std(d_true(tau_eta == noisetaus(i))),...
        'bo','capsize',0,'markerfacecolor','b'), 1:numel(noisetaus)); % plot(tau_eta,d_true,'b.'); 
    axis([0 2.5 200 600]); hline(x,'k--'); xlabel('noise time constant'); ylabel('stopping distance [cm]');
    
    subplot(2,2,3); hold on;
    arrayfun(@(i) errorbar(gains(i)-0.005,mean(d_true(gain_eta == gains(i) & tauindx)),std(d_true(gain_eta == gains(i) & tauindx)),...
        'o','color',colr(1,:),'capsize',0,'markerfacecolor',colr(1,:),'handlevisibility','off'), 1:numel(gains)); 
    arrayfun(@(i) errorbar(gains(i)+0.005,mean(d_true(gain_eta == gains(i) & ~tauindx)),std(d_true(gain_eta == gains(i) & ~tauindx)),...
        'o','color',colr(2,:),'capsize',0,'markerfacecolor',colr(2,:),'handlevisibility','off'), 1:numel(gains));
    plot(gains(1)-0.005,mean(d_true(gain_eta == gains(1) & tauindx)),'o','color',colr(1,:),'markerfacecolor',colr(1,:))
    plot(gains(1)+0.005,mean(d_true(gain_eta == gains(1) & ~tauindx)),'o','color',colr(2,:),'markerfacecolor',colr(2,:))
%     plot(gain_eta(tauindx),d_true(tauindx),'.','color',colr(1,:)); plot(gain_eta(~tauindx),d_true(~tauindx),'.','color',colr(2,:)); 
    axis([-0.05 0.15 100 700]); hline(x,'k--'); xlabel('noise gain'); ylabel('stopping distance [cm]'); legend({['tau = ' num2str(taugroup(1))],['tau = ' num2str(taugroup(2))]});
    
    subplot(2,2,4); hold on;
    arrayfun(@(i) errorbar(noisetaus(i)-0.1,mean(d_true(tau_eta == noisetaus(i) & tauindx)),std(d_true(tau_eta == noisetaus(i) & tauindx)),...
        'o','color',colr(1,:),'capsize',0,'markerfacecolor',colr(1,:),'handlevisibility','off'), 1:numel(noisetaus)); 
    arrayfun(@(i) errorbar(noisetaus(i)+0.1,mean(d_true(tau_eta == noisetaus(i) & ~tauindx)),std(d_true(tau_eta == noisetaus(i) & ~tauindx)),...
        'o','color',colr(2,:),'capsize',0,'markerfacecolor',colr(2,:),'handlevisibility','off'), 1:numel(noisetaus)); 
    plot(noisetaus(1)-0.1,mean(d_true(tau_eta == noisetaus(1) & tauindx)),'o','color',colr(1,:),'markerfacecolor',colr(1,:))
    plot(noisetaus(1)+0.1,mean(d_true(tau_eta == noisetaus(1) & ~tauindx)),'o','color',colr(2,:),'markerfacecolor',colr(2,:))
%     plot(tau_eta(tauindx),d_true(tauindx),'.','color',colr(1,:)); plot(tau_eta(~tauindx),d_true(~tauindx),'.','color',colr(2,:)); 
    axis([0 2.5 200 600]); hline(x,'k--'); xlabel('noise time constant'); ylabel('stopping distance [cm]'); legend({['tau = ' num2str(taugroup(1))],['tau = ' num2str(taugroup(2))]});
    
end

