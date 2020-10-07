function Tau_estim = MAPtau4staticPrior(tau,params,sampling_distribution,log_space)
% log_space: determines whether you formulate the distributions in log-space
%            (1) or normal-space (0). If 1, use phi = log(tau), instead of
%            just tau.

%%
if log_space
    
    phi = log(tau);
    x = -5:.01:3; % range of phi = log(tau) values
    
    % assign parameter values
    a = params(1); b = params(2); c = params(3); % [sig_L_tau  sig_prior  mu_prior]
    
    for j = 1:length(tau)
        
        measurement = phi(j);
        if 0
            % Likelihood function
            L_tau = (1/(a*sqrt(2*pi))).*exp(-((measurement - x).^2)./(2*a^2)); % normally distributed phi = log(tau)
            % Prior
            Prior_tau = (1./(b*sqrt(2*pi))).*exp(-((x - c).^2)./(2*b^2));
            % Posterior
            Pos_tau = Prior_tau.*L_tau;
            % MAP estimate
            [~,indx] = max(Pos_tau);
            phi_estim(j) = x(indx);
            Tau_estim(j) = exp(phi_estim(j));
        end
        % analytical posterior
        Pos_tau = (1/(a*b*2*pi)).*exp(-((measurement-x).^2)./(2*a^2) - ((x-c).^2)./(2*b^2));
        Pos_tau_dot = ((((measurement-x)./a^2)+((c-x)./b^2))./(a*b*2*pi)).*...
            exp(-((measurement-x).^2)./(2*a^2) - ((x-c).^2)./(2*b^2));
        phi_estim(j) = (b^2*measurement + a^2*c)./(a^2 + b^2);
        Tau_estim(j) = exp(phi_estim(j));
        
        
        % sanity check
        if 0
            figure;
            
            plot(x,L_tau/max(L_tau),'k','linewidth',3); vline(measurement,'k'); hold on;
            plot(x,Prior_tau/max(Prior_tau),'b','linewidth',3); vline(sampling_distribution.mu,'k');
            hold on; plot(x,Pos_tau/max(Pos_tau),'r'); vline(x(find(cumsum(Pos_tau)/sum(Pos_tau) >= 0.5, 1)),'b');
            vline(phi_estim(j),'r'); vline(Tau_estim(j),'r--');
            hold off;
        end
    end

    %%
else
    
    x = .01:.01:10; % range of tau values
    
    % assign parameter values
    a = params(1); b = params(2); c = params(3); % [sig_L_tau  sig_prior  mu_prior]
    
    for j = 1:length(tau)
        
        measurement = tau(j);
        if 0
            % Likelihood function
            L_tau = (1./measurement).*(1/(a*sqrt(2*pi))).*exp(-((log(measurement) - (log(x) - a.^2/2)).^2)./(2*a^2)); % log-normal
            % Prior
            Prior_tau = (1./x).*(1/(b*sqrt(2*pi))).*exp(-((log(x) - c).^2)./(2*b^2));
            % Posterior
            Pos_tau = Prior_tau.*L_tau;
            % MAP estimate
            [~,indx] = max(Pos_tau);
            Tau_estim(j) = x(indx);
        end
        % analytical posterior
        Pos_tau = (1./(2*pi*measurement*a*b*x)).*exp(-(log(x)-c).^2./(2*b^2) - (-log(x)+(a^2/2)+log(measurement)).^2./(2*a^2));
        Pos_tau_dot = (-1./4*pi*measurement*b^3*a^3.*x.^2).*...
            exp(-(log(x)-c).^2./(2*b^2) - (-log(x)+(a^2/2)+log(measurement)).^2./(2*a^2)).*...
            ((2*a^2+2*b^2).*log(x) + (b^2-2*c)*a^2 - 2*log(measurement)*b^2);
        Tau_estim(j) = exp((-b^2*a^2 + 2*c*a^2 + 2*log(measurement)*b^2)./(2*b^2+2*a^2));
        
        
        % sanity check
        if 0
            figure;
            
            plot(x,L_tau/max(L_tau),'k','linewidth',3); vline(measurement,'k'); hold on;
            plot(x,Prior_tau/max(Prior_tau),'b','linewidth',3); vline(sampling_distribution.mu,'k');
            hold on; plot(x,Pos_tau/max(Pos_tau),'r'); vline(x(find(cumsum(Pos_tau)/sum(Pos_tau) >= 0.5, 1)),'b');
            vline(phi_estim(j),'r'); vline(Tau_estim(j),'r--');
            hold off;
        end
    end
    
end
