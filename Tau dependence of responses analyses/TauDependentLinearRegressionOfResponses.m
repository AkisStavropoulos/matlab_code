function [b,h,p] = TauDependentLinearRegressionOfResponses(subject)

%% Linear regression with tau-dependent bias
% First (1) coefficient is for distance, second (2) is for distance*tau
params = 'stimtype'; 
[poolindx,legend_input] = get_poolindx(subject,params); 
colr = brewermap(size(poolindx,2),'Dark2');
b = [];
for i = 1:length(subject)
   for s = 1:size(poolindx,2) 
    indx = poolindx{i,s};
    r_tar = cell2mat(arrayfun(@(x) x.prs.r_tar,subject(i).trials(indx),'uniformoutput',false));
    r_sub = cell2mat(arrayfun(@(x) x.prs.r_sub,subject(i).trials(indx),'uniformoutput',false));
    th_tar = cell2mat(arrayfun(@(x) x.prs.th_tar,subject(i).trials(indx),'uniformoutput',false));
    th_sub = cell2mat(arrayfun(@(x) x.prs.th_sub,subject(i).trials(indx),'uniformoutput',false));
    tau = cell2mat(arrayfun(@(x) x.prs.tau,subject(i).trials(indx),'uniformoutput',false));
    % standardize quantities
    y_r = r_sub(:)./std(r_sub) ;
    y_th = th_sub(:)./std(th_sub) ;
    X_r = [r_tar(:)./std(r_tar)      (r_tar(:).*tau(:)) ./ std(r_tar(:).*tau(:))];
    X_th = [th_tar(:)./std(th_tar)      (th_tar(:).*tau(:)) ./ std(th_tar(:).*tau(:))];
    % regression
    b.r{i,s} = regress(y_r, X_r);
    b.th{i,s} = regress(y_th, X_th);
   end
end
b.r = cellfun(@(x) x',b.r,'uniformoutput',false);   b.th = cellfun(@(x) x',b.th,'uniformoutput',false);
bias = [];
nx = -.1:.01:1;
for s = 1:size(poolindx,2)
    % average
   bias.r.mu{s} = mean( cell2mat(b.r(:,s)) );     bias.r.se{s} = std( cell2mat(b.r(:,s)) )./sqrt(length(subject));
   bias.th.mu{s} = mean( cell2mat(b.th(:,s)) );     bias.th.se{s} = std( cell2mat(b.th(:,s)) )./sqrt(length(subject));

end
% test significance
[h.r.vis,p.r.vis] = ttest( cell2mat(b.r(:,1)), cell2mat(b.r(:,2)) );    
[h.r.comb,p.r.comb] = ttest( cell2mat(b.r(:,1)), cell2mat(b.r(:,3)) );    
[h.th.vis,p.th.vis] = ttest( cell2mat(b.th(:,1)), cell2mat(b.th(:,2)) );    
[h.th.comb,p.th.comb] = ttest( cell2mat(b.th(:,1)), cell2mat(b.th(:,3)) );    

% bar graph
figure; 
for s = 1:size(poolindx,2)
    subplot(1,2,1); hold on;
    bar(1+s,bias.r.mu{s}(1),'facecolor',colr(s,:),'edgecolor','none'); errorbar(1+s,bias.r.mu{s}(1),bias.r.se{s}(1),'k','capsize',0);
    bar(6+s,bias.r.mu{s}(2),'facecolor',colr(s,:),'edgecolor','none'); errorbar(6+s,bias.r.mu{s}(2),bias.r.se{s}(2),'k','capsize',0);
    xlabel('predictors'); ylabel('Regression coefficients'); title('Radial distance'); ylim([0 1]);
    xticks([3 8]);   xticklabels({'r_{tar}','\tau r_{tar}'}); 
    
    subplot(1,2,2); hold on; 
    bar(1+s,bias.th.mu{s}(1),'facecolor',colr(s,:),'edgecolor','none'); errorbar(1+s,bias.th.mu{s}(1),bias.th.se{s}(1),'k','capsize',0);
    bar(6+s,bias.th.mu{s}(2),'facecolor',colr(s,:),'edgecolor','none'); errorbar(6+s,bias.th.mu{s}(2),bias.th.se{s}(2),'k','capsize',0);
    xlabel('predictors'); ylabel('Regression coefficients'); title('Angular eccentricity'); ylim([0 1]);
    xticks([3 8]);   xticklabels({'\theta_{tar}','\tau \theta_{tar}'});
end
