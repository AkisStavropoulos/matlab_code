function [MSE_sum,pred_sub] = argminResErr(params,tau,sampling_distribution,JS,start_pos,exp_sub,modelname,log_space)
%% Compute tau estimate based on chosen model

if strcmp(modelname,'static')
    MAP_tau = MAPtau4staticPrior(tau,params,sampling_distribution,log_space);
    
elseif strcmp(modelname,'fixed')
    MAP_tau = params*ones(1,length(tau));
        
elseif strcmp(modelname,'dynamic')
    MAP_tau = MAPtau4dynamicPrior(tau,params,sampling_distribution,log_space);
    
end

% check real and estimated taus
if 0
    figure;plot(tau,MAP_tau,'x');hold on; plot(0:8,0:8,'k--')
end

% Calculate MSE of residual errors
[MSE_sum,pred_sub] = MSE_ResErr(MAP_tau,JS,start_pos,exp_sub);