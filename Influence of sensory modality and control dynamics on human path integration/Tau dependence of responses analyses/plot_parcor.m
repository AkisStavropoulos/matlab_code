function plot_parcor(parcor,full_model,legend_input)
%% plot partial correlation coefficients of multilinear model of tau and target location

Nsubjects = size(parcor.r.rho,1);
Nstim = size(parcor.r.rho,2);

colr = brewermap(Nstim,'Dark2');

figure; fonts = 8;
c = []; d = [];
for s = 1:Nstim
    subplot(2,Nstim,s); hold on;
    for i = 1:Nsubjects
        for n = 1:length(parcor.r.rho{i,s}.Row)-1
        c{s}(n,i) = parcor.r.rho{i,s}{n+1,1};
        d{s}(n,i) = parcor.th.rho{i,s}{n+1,1};
        end
    end
    h = bar(c{s},'FaceColor',colr(s,:));
    if full_model
   xticks([1 2 3 4 ]);   xticklabels({'c','r_{tar}','\tau','\tau*r_{tar}'}); 
    else
        xticks([1 2 3]);   xticklabels({'c','r_{tar}','\tau*r_{tar}'}); 
    end
    temp_mean = nanmean(c{s}(2:end,:),2);    temp_sd = nanstd(c{s}(2:end,:),[],2);  T = []; 
    for t = 1:length(temp_mean)
    T{t} = [num2str(temp_mean(t),'%.2f') '\pm' num2str(temp_sd(t),'%.2f')];
    end
    text(2:size(c{s},1),[-.5 -.75 -.5],T,'FontSize',fonts,'HorizontalAlignment','center')
    title(legend_input{s});    ylim([-1 1]);  ylabel('\rho');

    subplot(2,Nstim,s+Nstim); hold on;
    g = bar(d{s},'FaceColor',colr(s,:));
    if full_model
   xticks([1 2 3 4 ]);   xticklabels({'c','\theta_{tar}','\tau','\tau*\theta_{tar}'}); 
    else
        xticks([1 2 3]);   xticklabels({'c','\theta_{tar}','\tau*\theta_{tar}'}); 
    end
    temp_mean = nanmean(d{s}(2:end,:),2);    temp_sd = nanstd(d{s}(2:end,:),[],2);  T = [];
    for t = 1:length(temp_mean)
        T{t} = [num2str(temp_mean(t),'%.2f') '\pm' num2str(temp_sd(t),'%.2f')];
    end
    text(2:size(d{s},1),[-.5 -.75 -.5],T,'FontSize',fonts,'HorizontalAlignment','center')

   title(legend_input{s});    ylim([-1 1]);  ylabel('\rho');
end
if full_model
    sgtitle(['All subjects - \bfPartial correlations: \rm\alphar_{tar} + \beta\tau + \gamma\taur_{tar} + c' ]);
else
    sgtitle(['All subjects - \bfreduced model: \rm\alphar_{tar} + \gamma\taur_{tar}']);
end
