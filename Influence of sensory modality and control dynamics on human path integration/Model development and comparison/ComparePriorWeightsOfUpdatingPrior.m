function p_weight = ComparePriorWeightsOfUpdatingPrior(kalman_prior)
%% Compare differences between prior weights across conditions, for the updating prior model
% kalman_prior = importdata('kalmanprior_prs028.mat');

Nsubs = size(kalman_prior.prs.sd_tau,1);
Nstim = size(kalman_prior.prs.sd_tau,2);
Nboots = 100000;

for i = 1:Nsubs
    for s = 1:Nstim
        a = kalman_prior.prs.sd_tau(i,s);
        b = 0;
        p = kalman_prior.prs.sd_prior(i,s);
        
        k(i,s) = (p + b^2) ./ (p + b^2 + a^2);
    end
end

w_prior = 1-k;

% check differences
disp('***********DIFFERENCES OF PRIOR WEIGHTS***********')
tmp = mat2cell(w_prior, Nsubs, ones(Nstim,1));
bstrpmat = randi(Nsubs,Nsubs,Nboots);
tmp = cellfun(@(x) mean(x(bstrpmat)), tmp,'uniformoutput',false);
pval = sum(tmp{1} - tmp{2} < 0)/Nboots;     p_weight.bootstrap(1) = min([pval 1-pval]);
pval = sum(tmp{1} - tmp{3} < 0)/Nboots;     p_weight.bootstrap(2) = min([pval 1-pval]);
pval = sum(tmp{2} - tmp{3} < 0)/Nboots;     p_weight.bootstrap(3) = min([pval 1-pval]);
p_weight
