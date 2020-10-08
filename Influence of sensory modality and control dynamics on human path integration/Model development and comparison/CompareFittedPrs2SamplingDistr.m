function [p,h] = CompareFittedPrs2SamplingDistr(C,RealPrior)

%% Fitted parameters differences across modalities

Nsubs = size(C.sd_tau,1);
Nstim = size(C.sd_tau,2);
Nboots = 10000;

D.mu_prior = C.mu_prior - RealPrior.mu;
D.sd_prior = C.sd_prior - RealPrior.sig;

a = 1-(1-0.05)^(1/Nstim); % adjustment for multiple comparisons
a = 0.05;

disp('************COMPARE FITTED PARAMETERS WITH SAMPLING DISTRIBUTION************');
disp('PRIOR WIDTH:')
for s = 1:Nstim
    [p.sd_prior.signrank(s),h.sd_prior.signrank(s)] = signrank(D.sd_prior(:,s));
end

for s = 1:Nstim
    [h.sd_prior.ttest(s),p.sd_prior.ttest(s)] = ttest(D.sd_prior(:,s));
end

tmp = mat2cell(D.sd_prior, Nsubs, ones(Nstim,1));
bstrpmat = randi(Nsubs,Nsubs,Nboots);
tmp = cellfun(@(x) mean(x(bstrpmat)), tmp,'uniformoutput',false);
pval = cellfun(@(x) sum(x < 0)/Nboots, tmp);     p.sd_prior.bootstrap = min([pval ; 1-pval]);
% print whether statistical significant
structfun(@(x) x < a, p.sd_prior,'uniformoutput',false)

disp(' ')
disp('PRIOR MEAN:')
for s = 1:Nstim
    [p.mu_prior.signrank(s),h.mu_prior.signrank(s)] = signrank(D.mu_prior(:,s));
end

for s = 1:Nstim
    [h.mu_prior.ttest(s),p.mu_prior.ttest(s)] = ttest(D.mu_prior(:,s));
end

tmp = mat2cell(D.mu_prior, Nsubs, ones(Nstim,1));
bstrpmat = randi(Nsubs,Nsubs,Nboots);
tmp = cellfun(@(x) mean(x(bstrpmat)), tmp,'uniformoutput',false);
pval = cellfun(@(x) sum(x < 0)/Nboots, tmp);     p.mu_prior.bootstrap = min([pval ; 1-pval]);
% print whether statistical significant
structfun(@(x) x < a, p.mu_prior,'uniformoutput',false)



