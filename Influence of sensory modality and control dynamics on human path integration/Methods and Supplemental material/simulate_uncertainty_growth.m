function [b2,rho,pval] = simulate_uncertainty_growth(subject, power_par, interc, plt)
%% Simulate increase of uncertainty over trajectories
% Check for sub-linear, linear, supra-linear scaling of variance
if 0
interc = 10;
power_par = 2;
plt = 0;
end

a = power_par;
b = interc;

dt = 1/60;
params = 'stimtype';
[poolindx,legend_input] = get_poolindx(subject,params);

U = []; tau = []; r_sub = [];
for i = 1:length(subject)
    for s = 1:size(poolindx,2)
        indx = poolindx{i,s};
        ntrls = length(indx);
        
        for j = 1:ntrls
            r_sub{i,s}(j) = subject(i).trials(indx(j)).prs.r_sub;
            tau{i,s}(j) = subject(i).trials(indx(j)).prs.tau;
            vel = subject(i).trials(indx(j)).continuous.v; vel(vel<0) = nan;
            U{i,s}(j) = (nansum(vel.^a) + b).*dt; % units are (cm/s)^2 or cm^2 ??
            
        end
        b1{i,s} = regress(r_sub{i,s}(:),[ones(ntrls,1) tau{i,s}(:)]);
        b2{i,s} = regress(U{i,s}(:),[ones(ntrls,1) tau{i,s}(:)]);
        if any(~isreal(U{i,s}))
            rho(i,s) = nan; pval(i,s) = nan;
        else
        [rho(i,s),pval(i,s)] = corr(U{i,s}(:),tau{i,s}(:));
        end
    end
end

if plt
colr = brewermap(size(poolindx,2),'Dark2');
for i = 1:length(subject)
    figure;
    for s = 1:size(poolindx,2)
        
        subplot(2,size(poolindx,2),s); hold on;
        plot(tau{i,s},r_sub{i,s},'.','color',colr(s,:)); 
        plot(sort(tau{i,s}),[ones(numel(tau{i,s}),1) sort(tau{i,s}(:))]*b1{i,s},'k'); xlabel('tau [s]'); ylabel('travel distance [cm]');
        
        subplot(2,size(poolindx,2),s+size(poolindx,2)); hold on;
        plot(tau{i,s},U{i,s},'.','color',colr(s,:));
        plot(sort(tau{i,s}),[ones(numel(tau{i,s}),1) sort(tau{i,s}(:))]*b2{i,s},'k'); xlabel('tau [s]'); ylabel('accumulated uncertainty [cm^2/s]');
        
    end
    suptitle(subject(i).name)
end
end

b2 = cellfun(@(x) x(2),b2);