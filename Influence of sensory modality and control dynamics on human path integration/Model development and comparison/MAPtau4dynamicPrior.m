function Tau_estim = MAPtau4dynamicPrior(tau,params,sampling_distribution,log_space)
% log_space: determines whether you formulate the distributions in log-space
%            (1) or normal-space (0). If 1, use phi = log(tau), instead of
%            just tau.

%%
if log_space
        
phi = log(tau);
x = -5:.01:4; % range of phi = log(tau) values
p = zeros(1,length(tau)); k = zeros(1,length(tau)); % initialize kalman parameters

% assign parameter values
a = params(1); c = params(3); d = params(2); % [sig_L_tau   sig_prior   mu_prior]

% initial prior
Prior_tau = (1/(d*sqrt(2*pi))).*exp(-((x - c).^2)./(2*d^2));
prior_mu = c;

% estimated Prior initial variance
p(1) = d^2;

for j = 1:length(tau)
    
    % Likelihood function
    measurement = phi(j);
    L_tau = (1/(a*sqrt(2*pi))).*exp(-((measurement - x).^2)./(2*a^2)); % normally distributed phi = log(tau)
    % Posterior
    Pos_tau = Prior_tau.*L_tau;
    % MAP estimate
    [~,indx] = max(Pos_tau);
    phi_estim(j) = x(indx);
    Tau_estim(j) = exp(phi_estim(j));

    % Update Prior with Kalman gain
    Prior_old = Prior_tau;
    k(j+1) = p(j) ./ (p(j) + a^2);
    p(j+1) = p(j); % p(j+1) = k(j+1)*a^2;
    prior_mu = (1-k(j+1)).*prior_mu + k(j+1).*measurement;  % Prior_tau = (1-k(j+1)).*Prior_old + k(j+1).*L_tau;
    Prior_tau = (1/(d*sqrt(2*pi))).*exp(-((x - prior_mu).^2)./(2*d^2));

    % sanity check
    if 0
    figure;
    
    plot(x,L_tau/max(L_tau),'k','linewidth',3); vline(measurement,'k'); hold on; 
    plot(x,Prior_old/max(Prior_old),'b','linewidth',3); % vline(sampling_distribution.mu,'k');
    plot(x,Pos_tau/max(Pos_tau),'r'); vline(x(find(cumsum(Pos_tau)/sum(Pos_tau) >= 0.5, 1)),'b'); 
    vline(phi_estim(j),'r'); % vline(Tau_estim(j),'r--');
    plot(x,Prior_tau/max(Prior_tau),'k--'); xlim([-5 3]); xlabel('\phi = log(\tau)');
    hold off;
    end
end
%%  
else
    
x = .01:.01:10; % range of tau values
p = zeros(1,length(tau)); k = zeros(1,length(tau)); % initialize kalman parameters

% assign parameter values
a = params(1); c = params(3); d = params(2); % [sig_L_tau   sig_prior   mu_prior]

% initial prior
Prior_tau = (1./x).*(1/(d*sqrt(2*pi))).*exp(-((log(x) - c).^2)./(2*d^2));
prior_mu = c;

% estimated Prior initial variance
p(1) = d^2;

for j = 1:length(tau)
    %  sig_log = sqrt(log(1 + ((sig_tau/measurement)^2)));
    %  mu_log = log(measurement/sqrt(1 + ((sig_tau/measurement)^2)));
    
    % Likelihood function
    measurement = tau(j);
    L_tau = (1./measurement).*(1/(a*sqrt(2*pi))).*exp(-((log(measurement) - (log(x) - a.^2/2)).^2)./(2*a^2)); % log-normal
    % Posterior
    Pos_tau = Prior_tau.*L_tau;
    % MAP estimate
    [~,indx] = max(Pos_tau);
    Tau_estim(j) = x(indx);
    
    % Update Prior with Kalman gain
    Prior_old = Prior_tau;
    k(j+1) = p(j) ./ (p(j) + a^2);
    p(j+1) = p(j); % p(j+1) = k(j+1)*a^2;
    prior_mu = (1-k(j+1)).*prior_mu + k(j+1).*measurement;  % Prior_tau = (1-k(j+1)).*Prior_old + k(j+1).*L_tau;
    Prior_tau = (1./x).*(1/(d*sqrt(2*pi))).*exp(-((log(x) - prior_mu).^2)./(2*d^2));
    
    % sanity check
    if 0
    figure;
    
    plot(x,L_tau/max(L_tau),'k','linewidth',3); vline(measurement,'k'); hold on; 
    plot(x,Prior_old/max(Prior_old),'b','linewidth',3); vline(exp(sampling_distribution.mu + sampling_distribution.sig^2/2),'k');
    plot(x,Pos_tau/max(Pos_tau),'r'); vline(x(find(cumsum(Pos_tau)/sum(Pos_tau) >= 0.5, 1)),'b'); vline(Tau_estim(j),'r');
    plot(x,Prior_tau/max(Prior_tau),'k--'); xlim([0 6]);
    hold off;
    end
    
end

end
