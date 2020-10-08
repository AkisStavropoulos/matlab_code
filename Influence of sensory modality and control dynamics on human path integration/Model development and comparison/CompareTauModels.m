function [p] = CompareTauModels(rho1, rho2, model_labels)
%% Compare correlation of Residual Error and Tau before and after applying model
% rho: struct with fields "(modelnames{1})" and "(modelnames{2})"
% each field has subfields "r" and "th"

rot = 30;
colr = brewermap(size(rho1.r,2),'Dark2');
%% Correlation subject average
figure;
for s = 1:size(rho1.r,2)
    subplot(2,2,1); hold on;
    c = rho1.r(:,s); ylim([-.4 .6]);
    bar(s,mean(c),'FaceColor',colr(s,:),'edgecolor','none');
    errorbar(s,mean(c),std(c)/sqrt(size(rho1.r,1)),'k','capsize',0,'linestyle','none'); 
    xticks([1 2 3]);   xticklabels({'vestibular','visual','combined'}); % xtickangle(rot);
    title(['Distance - ' model_labels{1}])

    subplot(2,2,2); hold on;
    g = rho2.r(:,s); ylim([-.4 .6]);
    bar(s,mean(g),'FaceColor',colr(s,:)*0.5,'edgecolor','none');
    errorbar(s,mean(g),std(g)/sqrt(size(rho1.r,1)),'k','capsize',0,'linestyle','none');
    xticks([1 2 3]);   xticklabels({'vestibular','visual','combined'});
    title(['Distance - '  model_labels{2}])

    subplot(2,2,3); hold on;
    c = rho1.th(:,s); ylim([-.4 .6]);
    bar(s,mean(c),'FaceColor',colr(s,:),'edgecolor','none');
    errorbar(s,mean(c),std(c)/sqrt(size(rho1.r,1)),'k','capsize',0,'linestyle','none');
    xticks([1 2 3]);   xticklabels({'vestibular','visual','combined'});
    title(['Angle - ' model_labels{1}])

    subplot(2,2,4); hold on;
    g = rho2.th(:,s); ylim([-.4 .6]);
    bar(s,mean(g),'FaceColor',colr(s,:)*0.5,'edgecolor','none');
    errorbar(s,mean(g),std(g)/sqrt(size(rho1.r,1)),'k','capsize',0,'linestyle','none');
    xticks([1 2 3]);   xticklabels({'vestibular','visual','combined'});
    title(['Angle - '  model_labels{2}])
end
suptitle('Model corr(\tau,\epsilon_{res}) comparison')

%% Paired test of correlations
if 0
% sign rank
for s = 1:size(rho1.r,2)
    [p.corr.r(s),h.corr.r(s),stats.corr.r(s)] = signrank(rho1.r(:,s)-rho2.r(:,s));
    [p.corr.th(s),h.corr.th(s),stats.corr.th(s)] = signrank(rho1.th(:,s)-rho2.th(:,s));
end

for s = 1:size(rho1.r,2)
    [signrank(rho1.r(:,s))  -  signrank(rho2.r(:,s))];
    [signrank(rho1.th(:,s))  -  signrank(rho2.th(:,s))];
end
p.corr.test = 'signed-rank';
else
% t-test
for s = 1:size(rho1.r,2)
    [h.corr.r(s),p.corr.r(s),~,stats.corr.r(s)] = ttest(rho1.r(:,s)-rho2.r(:,s));
    [h.corr.th(s),p.corr.th(s),~,stats.corr.th(s)] = ttest(rho1.th(:,s)-rho2.th(:,s));
end

for s = 1:size(rho1.r,2)
    [signrank(rho1.r(:,s))  -  signrank(rho2.r(:,s))];
    [signrank(rho1.th(:,s))  -  signrank(rho2.th(:,s))];
end
p.corr.test = 't-test';
end

