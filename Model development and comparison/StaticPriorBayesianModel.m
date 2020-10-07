function [FitPrs,MSE] = StaticPriorBayesianModel(tau,JS,start_pos,exp_sub,prior,realindx,log_space,badsfit)
% First (1) parameter: sigma of likelihood
% Second (2) parameter: sigma of prior
% Third (3) parameter: mu of prior
%
% log_space = 1;
% badsfit = 1; 

%% Static prior Model fit
disp('++++++++++++++++++++ FITTING BAYESIAN MODEL WITH STATIC PRIOR ++++++++++++++++++++');

% Initial Parameters
ed_guess = importdata('fitprs188.mat');

B = []; Fval = []; MSE = []; sig_Likelihood = []; sig_Prior = []; mu_Prior = [];
for i = 1:size(realindx,1)
    for s = 1:size(realindx,2)
        
        ntrls = length(realindx{i,s});
        indx = 1:ntrls;
        
        tautemp = tau{i,s}(indx); JStemp = structfun(@(x) x(indx), JS{i,s},'uniformoutput',false);
        start_postemp = structfun(@(x) x(indx), start_pos{i,s},'uniformoutput',false);
        exp_subtemp = structfun(@(x) x(indx),exp_sub{i,s},'uniformoutput',false);
        
        % initialize parameters and bounds (add noise to initial guess)
        init_prs = [ed_guess.prs.sd_tau(i,s) ed_guess.prs.sd_prior(i,s) ed_guess.prs.mu_prior(i,s)]; %  [sig_L_tau, sig_prior, mu_prior]
        params0 = init_prs + 0.1*rand(1,3) - 0.05; % params(1): sig_L_tau, params(2): sig_prior, params(3): mu_prior
        LB = [0.0001 0.0001 -2]; UB = [1 1 1.8];
        LPB = [.001 .001 -1];  UPB = [.5 .5 1.5];
        
        % check bounds
        if any(params0 < LB) || any(params0 > UB)
            ind = find(params0 < LB); params0(ind) = LB(ind);
            ind = find(params0 > UB); params0(ind) = UB(ind);
        end
        % fit model
        if ~badsfit
            [BB,fval] = fminsearch(@(params) argminResErr(params,tautemp,prior,JStemp,start_postemp,exp_subtemp,'static',log_space),params0);
        else
            [BB,fval] = bads(@(params) argminResErr(params,tautemp,prior,JStemp,start_postemp,exp_subtemp,'static',log_space),params0,LB,UB,LPB,UPB);
        end
        
        sig_Likelihood(i,s) = BB(1); sig_Prior(i,s) = BB(2); mu_Prior(i,s) = BB(3);
    end
    Fval(i,s) = fval;
    disp(['.......Subject = ' num2str(i)])
   
end
B.sd_tau = sig_Likelihood;   B.sd_prior = sig_Prior;   B.mu_prior = mu_Prior;
MSE.val = Fval;
FitPrs.init_prs = init_prs;
FitPrs.prs = B;
FitPrs.mse = MSE;
FitPrs.logspace = log_space;
