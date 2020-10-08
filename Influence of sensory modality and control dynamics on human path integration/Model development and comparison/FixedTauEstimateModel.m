function [fixedtau,MSE_fixedtau] = FixedTauEstimateModel(tau,JS,start_pos,exp_sub,realindx)

%% Fixed tau model

disp('++++++++++++++++++++ FITTING FIXED TAU ESTIMATE MODEL ++++++++++++++++++++');

% initial parameter
tau_est0 = 2;

% fit model
for i = 1:size(realindx,1)
    for s = 1:size(realindx,2)
        
        [fixedtau(i,s),MSE_fixedtau(i,s)] = fminsearch(@(tau_est) argminResErr(tau_est,tau{i,s},[],JS{i,s},start_pos{i,s},exp_sub{i,s},'fixed',[]), tau_est0);
        
    end
    disp(['.......Subject = ' num2str(i)])
end

