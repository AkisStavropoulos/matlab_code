function [FitPrs,MSE] = DynamicPriorBayesianModel(tau,JS,start_pos,exp_sub,sampling_distribution,realindx,log_space,badsfit)

%% Kalman-like update of prior in the next trial
% Investigate whether prior should be dynamic
% First (1) parameter: sigma of likelihood
% Second (2) parameter: sigma of prior
% Third (3) parameter: mu of prior

init_prs = [.5 .2  1]; % [sig_L_tau   sig_prior   mu_prior]

B = []; Fval = []; MSE = []; sig_Likelihood = []; sig_Prior = []; mu_Prior = [];
for i = 1:size(realindx,1)
    for s = 1:size(realindx,2)
        ntrls = length(realindx{i,s});
        indx = 1:ntrls;
        
        tautemp = tau{i,s}(indx); JStemp = structfun(@(x) x(indx), JS{i,s},'uniformoutput',false);
        start_postemp = structfun(@(x) x(indx), start_pos{i,s},'uniformoutput',false);
        exp_subtemp = structfun(@(x) x(indx),exp_sub{i,s},'uniformoutput',false);
        
        % fit model
        LB = [.0001 .0001 -2];   UB = [50 10 2.5];
        LPB = [.05 .05 -1.5];  UPB = [40 8 2];
        params0 = init_prs;
        if ~badsfit
            [BB,fval] = fminsearch(@(params) argminResErr(params,tautemp,sampling_distribution,JStemp,start_postemp,exp_subtemp,'dynamic',log_space),params0);
        else
            % check bounds
            if any(params0 < LB) || any(params0 > UB)
                ind = find(params0 < LB); params0(ind) = LB(ind);
                ind = find(params0 > UB); params0(ind) = UB(ind);
            end
            [BB,fval] = bads(@(params) argminResErr(params,tautemp,sampling_distribution,JStemp,start_postemp,exp_subtemp,'dynamic',log_space),params0,LB,UB,LPB,UPB);
        end
        sig_Likelihood(i,s) = BB(1); sig_Prior(i,s) = BB(2); mu_Prior(i,s) = BB(3);
        Fval(i,s) = fval;
        
    end
    disp(['.......Subject = ' num2str(i)])
end

B.sd_tau = sig_Likelihood;   B.sd_prior = sig_Prior;   B.mu_prior = mu_Prior;

MSE.val = Fval;
FitPrs.init_prs = init_prs;
FitPrs.prs = B;
FitPrs.mse = MSE;
FitPrs.logspace = log_space;

