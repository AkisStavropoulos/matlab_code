function [p,h] = CompareFittedPrs(C)

%% Fitted parameters differences across modalities
% Order of comparison: 
% (1): vestibular vs visual
% (2): vestibular vs combined
% (3): visual vs combined

Nsubs = size(C.sd_tau,1);
Nstim = size(C.sd_tau,2);
Nboots = 100000;

disp('************COMPARE FITTED PARAMETERS ACROSS CONDITIONS************');
disp('LIKELIHOOD WIDTH:')
disp('SIGN RANK')
s = 1;
[p.sd_tau.signrank(1),h.sd_tau.signrank(1)] = signrank(C.sd_tau(:,s), C.sd_tau(:,s+1));
[p.sd_tau.signrank(2),h.sd_tau.signrank(2)] = signrank(C.sd_tau(:,s), C.sd_tau(:,s+2));
s = 2;
[p.sd_tau.signrank(3),h.sd_tau.signrank(3)] = signrank(C.sd_tau(:,s), C.sd_tau(:,s+1));
p.sd_tau.signrank

disp('T-TEST')
s = 1;
[h.sd_tau.ttest(1),p.sd_tau.ttest(1)] = ttest(C.sd_tau(:,s),C.sd_tau(:,s+1));
[h.sd_tau.ttest(2),p.sd_tau.ttest(2)] = ttest(C.sd_tau(:,s),C.sd_tau(:,s+2));
s = 2;
[h.sd_tau.ttest(3),p.sd_tau.ttest(3)] = ttest(C.sd_tau(:,s),C.sd_tau(:,s+1));
p.sd_tau.ttest

disp('BOOTSTRAP')
tmp = mat2cell(C.sd_tau, Nsubs, ones(Nstim,1));
bstrpmat = randi(Nsubs,Nsubs,Nboots);
tmp = cellfun(@(x) mean(x(bstrpmat)), tmp,'uniformoutput',false);
pval = sum(tmp{1} - tmp{2} < 0)/Nboots;     p.sd_tau.bootstrap(1) = min([pval 1-pval]);
pval = sum(tmp{1} - tmp{3} < 0)/Nboots;     p.sd_tau.bootstrap(2) = min([pval 1-pval]);
pval = sum(tmp{2} - tmp{3} < 0)/Nboots;     p.sd_tau.bootstrap(3) = min([pval 1-pval]);
p.sd_tau.bootstrap

disp(' ')
disp('PRIOR WIDTH:')
disp('SIGN RANK')
s = 1;
[p.sd_prior.signrank(1),h.sd_prior.signrank(1)] = signrank(C.sd_prior(:,s), C.sd_prior(:,s+1));
[p.sd_prior.signrank(2),h.sd_prior.signrank(2)] = signrank(C.sd_prior(:,s), C.sd_prior(:,s+2));
s = 2;
[p.sd_prior.signrank(3),h.sd_prior.signrank(3)] = signrank(C.sd_prior(:,s), C.sd_prior(:,s+1));
p.sd_prior.signrank

disp('T-TEST')
s = 1;
[h.sd_prior.ttest(1),p.sd_prior.ttest(1)] = ttest(C.sd_prior(:,s), C.sd_prior(:,s+1));
[h.sd_prior.ttest(2),p.sd_prior.ttest(2)] = ttest(C.sd_prior(:,s), C.sd_prior(:,s+2));
s = 2;
[h.sd_prior.ttest(3),p.sd_prior.ttest(3)] = ttest(C.sd_prior(:,s), C.sd_prior(:,s+1));
p.sd_prior.ttest

disp('BOOTSTRAP')
tmp = mat2cell(C.sd_prior, Nsubs, ones(Nstim,1));
bstrpmat = randi(Nsubs,Nsubs,Nboots);
tmp = cellfun(@(x) mean(x(bstrpmat)), tmp,'uniformoutput',false);
pval = sum(tmp{1} - tmp{2} < 0)/Nboots;     p.sd_prior.bootstrap(1) = min([pval 1-pval]);
pval = sum(tmp{1} - tmp{3} < 0)/Nboots;     p.sd_prior.bootstrap(2) = min([pval 1-pval]);
pval = sum(tmp{2} - tmp{3} < 0)/Nboots;     p.sd_prior.bootstrap(3) = min([pval 1-pval]);
p.sd_prior.bootstrap

disp(' ')
disp('PRIOR MEAN:')
disp('SIGN RANK')
s = 1;
[p.mu_prior.signrank(1),h.mu_prior.signrank(1)] = signrank(C.mu_prior(:,s), C.mu_prior(:,s+1));
[p.mu_prior.signrank(2),h.mu_prior.signrank(2)] = signrank(C.mu_prior(:,s), C.mu_prior(:,s+2));
s = 2;
[p.mu_prior.signrank(3),h.mu_prior.signrank(3)] = signrank(C.mu_prior(:,s), C.mu_prior(:,s+1));
p.mu_prior.signrank

disp('T-TEST')
s = 1;
[h.mu_prior.ttest(1),p.mu_prior.ttest(1)] = ttest(C.mu_prior(:,s), C.mu_prior(:,s+1));
[h.mu_prior.ttest(2),p.mu_prior.ttest(2)] = ttest(C.mu_prior(:,s), C.mu_prior(:,s+2));
s = 2;
[h.mu_prior.ttest(3),p.mu_prior.ttest(3)] = ttest(C.mu_prior(:,s), C.mu_prior(:,s+1));
p.mu_prior.ttest

disp('BOOTSTRAP')
tmp = mat2cell(C.mu_prior, Nsubs, ones(Nstim,1));
bstrpmat = randi(Nsubs,Nsubs,Nboots);
tmp = cellfun(@(x) mean(x(bstrpmat)), tmp,'uniformoutput',false);
pval = sum(tmp{1} - tmp{2} < 0)/Nboots;     p.mu_prior.bootstrap(1) = min([pval 1-pval]);
pval = sum(tmp{1} - tmp{3} < 0)/Nboots;     p.mu_prior.bootstrap(2) = min([pval 1-pval]);
pval = sum(tmp{2} - tmp{3} < 0)/Nboots;     p.mu_prior.bootstrap(3) = min([pval 1-pval]);
p.mu_prior.bootstrap


