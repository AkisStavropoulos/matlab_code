function [bias_mu,bias_sem,bias] = bias_taugroups(subject,params,button_push)

%% Radial vs angular bias for different TAU BINS for each modality
[poolindx,legend_input] = get_poolindx(subject,params);

[bias,skipsub] = get_bias(subject,params);

subject(skipsub) = [];
bias.r(skipsub,:) = [];             bias.th(skipsub,:) = [];

count = 1;
% plot biases
colr = brewermap(3,'Dark2'); colr = repmat(colr(:),1,size(poolindx,2)/3)'; colr = reshape(colr,[],3);
if size(poolindx,2)/3 == 3
    taushape = {'o';'o';'o'}; markersize = [2;10;20];
elseif size(poolindx,2)/3 == 4
    taushape = {'o';'o';'o';'o'};  markersize = [2;5;10;20];
end
taushape = repmat(taushape,size(poolindx,2)/3,1);
markersize = repmat(markersize,size(poolindx,2)/3,1);
count1 = 1; count2 = 1;
bias_mu = []; bias_sem = [];
figure; hold on;    axis equal;
for s = 1:size(bias.r,2)
    if size(bias.r,2) > 3;  cnd = count1;  taugrp = count2;
    else;                cnd = s;       taugrp = count2;    end

    bias_sem.r(cnd,taugrp) = std(bias.r(:,s))/sqrt(length(subject) - length(skipsub));
    bias_sem.th(cnd,taugrp) = std(bias.th(:,s))/sqrt(length(subject) - length(skipsub));
    bias_mu.r(cnd,taugrp) = mean(bias.r(:,s));
    bias_mu.th(cnd,taugrp) = mean(bias.th(:,s));
    errorbar(bias_mu.r(cnd,taugrp),bias_mu.th(cnd,taugrp),bias_sem.r(cnd,taugrp),bias_sem.r(cnd,taugrp),bias_sem.th(cnd,taugrp),bias_sem.th(cnd,taugrp),...
        taushape{s},'Color',colr(s,:),'MarkerSize',markersize(s),'markerfacecolor',colr(s,:));
    
    xlim([0 1]);ylim([0 1]);    hline(1,'k--'); vline(1,'k--');xlabel('radial bias');ylabel('angular bias');
    plot(0:2,0:2,'k','HandleVisibility','off'); axis([0.5 1.2 0.5 1.2]);
    if size(bias.r,2) > 3
        if mod(s,3) == 0
            count1 = count1 + 1;
            count2 = 1;
        else
            count2 = count2 + 1;
        end
    end
end
if ~(size(bias.r,2)>3); [~,legend_input] = get_poolindx(subject,[]); end
legend(legend_input{:},'Location','southeast');
xlabel('radial bias');ylabel('angular bias'); 
if button_push ; suptitle('1st button push'); end
