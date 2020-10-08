function CompareTauCorr(rho1,rho2)
%% Compare correlation of Residual Error and Tau before and after applying model
% rho1: data correlation
% rho2: must be the model produced correlation
% each field has subfields "r" and "th"

colr = brewermap(size(rho1.r,2),'Dark2');

%% Connecting-Lines plots
figure; hold on; f = []; g = [];
for s = 1:size(rho1.r,2)
    
    subplot(2,size(rho1.r,2),s); hold on;
    for i = 1:size(rho1.r,1)
        plot([rho1.r(i,s) rho2.r(i,s)],'-o','color',colr(s,:),'markerfacecolor',colr(s,:),'markeredgecolor','none');
    end
    xticks([1 2]);   xticklabels({'actual','subj.'});
    ylabel('correlation \rho'); axis equal; xlim([.5 2.5]); ylim([-1 1]);
    hline(0,'k');
    title('Compare corr (actual/subjective)');
    
    subplot(2,size(rho1.r,2),s+size(rho1.r,2)); hold on;
    for i = 1:size(rho1.r,1)
        plot([rho1.th(i,s) rho2.th(i,s)],'-o','color',colr(s,:),'markerfacecolor',colr(s,:),'markeredgecolor','none');
    end
    xticks([1 2]);   xticklabels({'actual','subj.'});
    ylabel('correlation \rho'); axis equal; xlim([.5 2.5]); ylim([-1 1]);
    hline(0,'k');
    title('Compare corr (actual/subjective)');
end

%% Unity line plots
if 0
    figure; hold on; f = []; g = [];
    for s = 1:size(rho1.r,2)
        
        subplot(2,size(rho1.r,2),s); hold on;
        plot(rho1.r(:,s),rho2.r(:,s),'o','markerfacecolor',colr(s,:),'markeredgecolor','k');
        plot(-1:1,-1:1,'k--','handlevisibility','off'); hline(0,'k');
        xlabel('actual'); ylabel('subjective'); axis equal; xlim([-1 1]); ylim([-1 1]);
        title('Compare corr (after/before)');

        subplot(2,size(rho1.r,2),s+size(rho1.r,2)); hold on;
        plot(rho1.th(:,s),rho2.th(:,s),'o','markerfacecolor',colr(s,:),'markeredgecolor','k');
        plot(-1:1,-1:1,'k--','handlevisibility','off'); hline(0,'k');
        xlabel('actual'); ylabel('subjective'); axis equal; xlim([-1 1]); ylim([-1 1]);
        title('Compare corr (actual/subjective)');
    end
end

%% subject average
figure;
for s = 1:size(rho1.r,2)
    subplot(2,2,1); hold on;
    c{i} = rho1.r(:,s); ylim([-.4 .6]);
    g = bar(s,mean(c{i}),'FaceColor',colr(s,:),'edgecolor','none');
    errorbar(s,mean(c{i}),std(c{i})/sqrt(size(rho1.r,1)),'k','capsize',0,'linestyle','none'); 
    xticks([1 2 3]);   xticklabels({'vestibular','visual','combined'});
    title('Distance - Actual corr')

    subplot(2,2,2); hold on;
    c{i} = rho2.r(:,s); ylim([-.4 .6]);
    g = bar(s,mean(c{i}),'FaceColor',colr(s,:),'edgecolor','none');
    errorbar(s,mean(c{i}),std(c{i})/sqrt(size(rho1.r,1)),'k','capsize',0,'linestyle','none');
    xticks([1 2 3]);   xticklabels({'vestibular','visual','combined'});
    title('Distance - Subjective corr')

    subplot(2,2,3); hold on;
    c{i} = rho1.th(:,s); ylim([-.4 .6]);
    g = bar(s,mean(c{i}),'FaceColor',colr(s,:),'edgecolor','none');
    errorbar(s,mean(c{i}),std(c{i})/sqrt(size(rho1.r,1)),'k','capsize',0,'linestyle','none');
    xticks([1 2 3]);   xticklabels({'vestibular','visual','combined'});
    title('Angle - Actual corr')

    subplot(2,2,4); hold on;
    c{i} = rho2.th(:,s); ylim([-.4 .6]);
    g = bar(s,mean(c{i}),'FaceColor',colr(s,:),'edgecolor','none');
    errorbar(s,mean(c{i}),std(c{i})/sqrt(size(rho1.r,1)),'k','capsize',0,'linestyle','none');
    xticks([1 2 3]);   xticklabels({'vestibular','visual','combined'});
    title('Angle - Subjective corr')
end
