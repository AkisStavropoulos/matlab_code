function [rho,pval] = trl_duration_vs_tau(subject,params)

%% Plot trial duration and average velocity vs tau
[poolindx,legend_input] = get_poolindx(subject,params);
colr = brewermap(size(poolindx,2),'Dark2');

trldur = []; tau = []; avg_vel = [];
for i = 1:length(subject)
    for s = 1:size(poolindx,2)
        indx = poolindx{i,s};
        for j = 1:length(indx)
            if any(~isnan(subject(i).trials(indx(j)).mc.flag2))
                
%             dt = subject(i).trials(indx(j)).continuous.ts(2) - subject(i).trials(indx(j)).continuous.ts(1);
%             dist = cumsum(subject(i).trials(indx(j)).continuous.v)*dt;
%             stopindx = find(dist > 0.8*dist(end), 1);
%             stopindx = find(subject(i).trials(indx(j)).continuous.v > 10, 1, 'last');
            stopindx = find(subject(i).trials(indx(j)).mc.flag2, 1)-0;
            if ~isempty(stopindx)
            trldur{i,s}(j) = subject(i).trials(indx(j)).continuous.ts(stopindx);  
            tau{i,s}(j) = subject(i).trials(indx(j)).prs.tau;
            avg_vel{i,s}(j) = mean(subject(i).trials(indx(j)).continuous.v(60:stopindx));
            else
            trldur{i,s}(j) = nan;
            tau{i,s}(j) = nan;
            avg_vel{i,s}(j) = nan;
            end
            
            else
            trldur{i,s}(j) = nan;
            tau{i,s}(j) = nan;  
            avg_vel{i,s}(j) = nan;
            end
        end
        if sum(isnan(tau{i,s})) < length(indx)
        [rho.t(i,s),pval.t(i,s)] = nancorr(tau{i,s}(:),trldur{i,s}(:));
        [rho.v(i,s),pval.v(i,s)] = nancorr(tau{i,s}(:),avg_vel{i,s}(:));
        else
            rho.t(i,s) = nan;
            pval.t(i,s) = nan;
            rho.v(i,s) = nan;
            pval.v(i,s) = nan;
        end
    end
end

% vertical scatter
figure;
subplot(1,2,1); hold on;
shake = 0.4;
c = [];
for s = 1:size(poolindx,2)
    % significant
    statindx = pval.t(:,s) < 0.05;
    c{s} = rho.t(statindx,s);
    plot(s*ones(length(c{s}),1) + shake*rand(length(c{s}),1), c{s}(:),'o','markerEdgeColor',colr(s,:),'markerfacecolor',colr(s,:));
    % not significant
    d{s} = rho.t(~statindx,s);
    plot(s*ones(length(d{s}),1) + shake*rand(length(d{s}),1), d{s}(:),'o','markerEdgeColor',colr(s,:));
end
title('tau and travel duration correlation'); axis([0 4 -1 1]); hline(0,'k');
xticks([1 2 3]);   xticklabels(legend_input);    ylabel('corr. coefficients');    

subplot(1,2,2); hold on;
c = [];
for s = 1:size(poolindx,2)
    % significant
    statindx = pval.v(:,s) < 0.05;
    c{s} = rho.v(statindx,s);
    plot(s*ones(length(c{s}),1) + shake*rand(length(c{s}),1), c{s}(:),'o','markerEdgeColor',colr(s,:),'markerfacecolor',colr(s,:));
    % not significant
    d{s} = rho.v(~statindx,s);
    plot(s*ones(length(d{s}),1) + shake*rand(length(d{s}),1), d{s}(:),'o','markerEdgeColor',colr(s,:));
end
title('tau and average velocity correlation'); axis([0 4 -1 1]); hline(0,'k');
xticks([1 2 3]);   xticklabels(legend_input);    ylabel('corr. coefficients');    


