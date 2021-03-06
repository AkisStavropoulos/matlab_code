%% GenerateSourceData

%% Load Data
data_folder = 'C:\Users\ges6\Documents\MATLAB\Data\Single firefly with motion cuing\Subjects\New\';
cd(data_folder);
folders = dir(data_folder);   folders = folders([folders.isdir]);

indx = cell2mat(arrayfun(@(x) ~strcmp(x.name(1),'.'),folders,'uniformoutput',false));
subject_name = {folders(indx).name};

[subject_backup,~] = LoadData(data_folder,subject_name);

%% All responses (Figure 2B)
clear responses_table
params = 'stimtype';
[poolindx,legend_input] = get_poolindx(subject,params); 
Nstim = size(poolindx,2);

[r,th,tau] = GetResponses(subject,params);

for s = 1:Nstim
    responses_table.distance.response(:,s) = table(r.sub(:,s));
    responses_table.angle.response(:,s) = table(th.sub(:,s));
    responses_table.distance.target(:,s) = table(r.tar(:,s));
    responses_table.angle.target(:,s) = table(th.tar(:,s));
    responses_table.tau(:,s) = table(tau(:,s));    
end
responses_table.distance.response.Properties.VariableNames = legend_input;
responses_table.angle.response.Properties.VariableNames = legend_input;
responses_table.distance.target.Properties.VariableNames = legend_input;
responses_table.angle.target.Properties.VariableNames = legend_input;
responses_table.tau.Properties.VariableNames = legend_input;
responses_table.readme = {'Figure 2B';'Each row contains data from one subject'};

save('responses_table','responses_table');

%% Response gain figure data (Figure 2B,2C,2D,3A)
clear response_gain_table

% Sensory conditions
polyorder = 1;  plt = 0;    intercept = 0;
params = 'stimtype';
[poolindx,legend_input] = get_poolindx(subject,params); 
[bias.r,bias.th] = ErrScatterFit(subject_backup,params,polyorder,intercept,plt);
Nstim = size(poolindx,2);

for s = 1:Nstim
    response_gain_table.distance.conditions(:,s) = table(bias.r(:,s));
    response_gain_table.angle.conditions(:,s) = table(bias.th(:,s));    
end
response_gain_table.distance.conditions.Properties.VariableNames = legend_input;
response_gain_table.angle.conditions.Properties.VariableNames = legend_input;

% Tau groups
polyorder = 1;  plt = 0;    intercept = 0;
params = 'stimtau';
[poolindx,legend_input] = get_poolindx(subject,params); 
[bias.r,bias.th] = ErrScatterFit(subject_backup,params,polyorder,intercept,plt);
Nstim = size(poolindx,2);
sensorygroups = {'vestibular_','visual_','combined_'};
taugroups = {'small_tau','medium_tau','large_tau'};
cnt = 1;
for k = 1:length(sensorygroups)
    for l = 1:length(taugroups)
        condgrps{cnt} = [sensorygroups{k} taugroups{l}];
        cnt = cnt+1;
    end
end

for s = 1:Nstim
    response_gain_table.distance.taugroups(:,s) = table(bias.r(:,s));
    response_gain_table.angle.taugroups(:,s) = table(bias.th(:,s));
end
response_gain_table.distance.taugroups.Properties.VariableNames = condgrps;
response_gain_table.angle.taugroups.Properties.VariableNames = condgrps;
response_gain_table.readme = {'Figure 2B,2C,2D,3A'};

save('response_gain_table','response_gain_table');
%% Residual errors data (Figure 3B,3C)
clear residual_errors_table
params = 'stimtype';
[poolindx,legend_input] = get_poolindx(subject,params); 
Nsubs = length(subject);
Nstim = size(poolindx,2);

[r_res,th_res,tau] = ComputeResidualErrors(subject,params);
for s = 1:Nstim
    residual_errors_table.distance.val(:,s) = table(r_res(:,s));
    residual_errors_table.angle.val(:,s) = table(th_res(:,s));
    residual_errors_table.tau.val(:,s) = table(tau(:,s));
end
residual_errors_table.distance.val.Properties.VariableNames = legend_input;
residual_errors_table.angle.val.Properties.VariableNames = legend_input;
residual_errors_table.tau.val.Properties.VariableNames = legend_input;

rho = []; pval = [];
[rho.r,rho.th,pval.r,pval.th] = ResErr_TauCorr(subject,params,0,0);

for s = 1:Nstim
    residual_errors_table.distance.rho(:,s) = table(rho.r(:,s));
    residual_errors_table.angle.rho(:,s) = table(rho.th(:,s));
    residual_errors_table.distance.pval(:,s) = table(pval.r(:,s));
    residual_errors_table.angle.pval(:,s) = table(pval.th(:,s));
end
residual_errors_table.distance.rho.Properties.VariableNames = legend_input;
residual_errors_table.angle.rho.Properties.VariableNames = legend_input;
residual_errors_table.distance.pval.Properties.VariableNames = legend_input;
residual_errors_table.angle.pval.Properties.VariableNames = legend_input;
residual_errors_table.readme = {'Figure 3B,3C';'In the last field "val", each row contains data from one subject'};

save('residual_errors_table','residual_errors_table');

%% Linear regression (standardized) of responses with tau-dependent terms (Figure 3D)
clear regr_coef_table
regr_coef = TauDependentLinearRegressionOfResponses(subject);

for s = 1:Nstim
    regr_coef_table.distance(:,s) = table(regr_coef.r(:,s));
    regr_coef_table.angle(:,s) = table(regr_coef.th(:,s));
end
regr_coef_table.distance.Properties.VariableNames = legend_input;
regr_coef_table.angle.Properties.VariableNames = legend_input;
regr_coef_table.readme = {'Figure 3D';'1st element corresponds to target position coefficient ';'2nd element corresponds to interaction term '};

save('regr_coef_table','regr_coef_table');

%% Residual errors correlation pre- and post-model (Figure 4B,6)
clear model_based_corr_table

log_space = 1;
rho = []; pval = []; MSE = []; pred_sub = [];
% Data
[tau,JS,sub,start_pos,exp_sub,~,targ] = extract_data_for_model(subject,'stimtype');
[exp_sub,start_pos,sub] = struct_within_cell_fix(exp_sub,start_pos,sub);

exp_x = cellfun(@(x,y,z) x.*sind(y) + z, exp_sub.r,exp_sub.th,start_pos.x,'uniformoutput',false); % exp_sub.r.*sind(exp_sub.th) + start_pos.x;
exp_y = cellfun(@(x,y,z) x.*cosd(y) + z, exp_sub.r,exp_sub.th,start_pos.y,'uniformoutput',false); % exp_sub.r.*cosd(exp_sub.th) + start_pos.y;
MSE.data = cellfun(@(x1,x2,y1,y2) median(sqrt((x1 - x2).^2 + (y1 - y2).^2)), sub.x, exp_x, sub.y, exp_y); % median(sqrt((sub_x - exp_x).^2 + (sub_y - exp_y).^2));
r_res = cellfun(@(x,y) x - y, sub.r, exp_sub.r,'uniformoutput',false);  
th_res = cellfun(@(x,y) sign(y).*(x - y), sub.th, exp_sub.th,'uniformoutput',false);
[rho.data.r, pval.data.r] = cellfun(@(x,y) nancorr(x(:),y(:)),tau,r_res); % make sure they're columns
[rho.data.th, pval.data.th] = cellfun(@(x,y) nancorr(x(:),y(:)),tau,th_res);

% Static Prior model
finalprs = importdata('fitprs188.mat');

[tau,JS,sub,start_pos,exp_sub,realindx] = extract_data_for_model(subject,'stimtype');
sampling_distribution.sig = (subject(1).trials(1).prs.sig_phi);     sampling_distribution.mu = (subject(1).trials(1).prs.mu_phi);
pred_sub.static = [];
for i = 1:length(subject)
    for s = 1:size(realindx,2)
        pars = [finalprs.prs.sd_tau(i,s) finalprs.prs.sd_prior(i,s) finalprs.prs.mu_prior(i,s)];
        [MSE.static(i,s),pred_sub.static{i,s}] = argminResErr(pars,tau{i,s},sampling_distribution,JS{i,s},start_pos{i,s},exp_sub{i,s},'static',log_space);
    end
end
[pred_sub.static,exp_sub,start_pos,sub] = struct_within_cell_fix(pred_sub.static,exp_sub,start_pos,sub);

r_res = cellfun(@(x,y) x - y, pred_sub.static.r, exp_sub.r,'uniformoutput',false); 
th_res = cellfun(@(x,y) sign(y).*(x - y), pred_sub.static.th, exp_sub.th,'uniformoutput',false);
[rho.static.r, pval.static.r] = cellfun(@(x,y) corr(x(:),y(:)),tau,r_res); % make sure they're columns
[rho.static.th, pval.static.th] = cellfun(@(x,y) corr(x(:),y(:)),tau,th_res);

% Dynamic Prior model 
dynamic_prior = importdata('kalmanprior_prs028.mat');

[tau,JS,sub,start_pos,exp_sub,realindx] = extract_data_for_model(subject,'stimtype');
sampling_distribution.sig = (subject(1).trials(1).prs.sig_phi);     sampling_distribution.mu = (subject(1).trials(1).prs.mu_phi);
pred_sub.dynamic = [];
for i = 1:length(subject)
    for s = 1:size(realindx,2)
        pars = [dynamic_prior.prs.sd_tau(i,s)   dynamic_prior.prs.sd_prior(i,s)   dynamic_prior.prs.mu_prior(i,s)];
        [MSE.dynamic(i,s),pred_sub.dynamic{i,s}] = argminResErr(pars,tau{i,s},sampling_distribution,JS{i,s},start_pos{i,s},exp_sub{i,s},'dynamic',log_space);
    end
end
[pred_sub.dynamic,exp_sub,start_pos,sub] = struct_within_cell_fix(pred_sub.dynamic,exp_sub,start_pos,sub);

r_res = cellfun(@(x,y) x - y, pred_sub.dynamic.r, exp_sub.r,'uniformoutput',false); 
th_res = cellfun(@(x,y) sign(y).*(x - y), pred_sub.dynamic.th, exp_sub.th,'uniformoutput',false);
[rho.dynamic.r, pval.dynamic.r] = cellfun(@(x,y) corr(x(:),y(:)),tau,r_res); % make sure they're columns
[rho.dynamic.th, pval.dynamic.th] = cellfun(@(x,y) corr(x(:),y(:)),tau,th_res);

% Fixed tau model 
varname = 'fixedtaus.mat';
fixedtau = importdata(varname);

[tau,JS,sub,start_pos,exp_sub,realindx] = extract_data_for_model(subject,'stimtype');
pred_sub.fixedtau = [];
for i = 1:size(realindx,1)
    for s = 1:size(realindx,2)
        [MSE.fixedtau(i,s),pred_sub.fixedtau{i,s}] = FixedTauModelResErr(fixedtau(i,s),JS{i,s},start_pos{i,s},exp_sub{i,s});
    end
end
[pred_sub.fixedtau,exp_sub,start_pos,sub] = struct_within_cell_fix(pred_sub.fixedtau,exp_sub,start_pos,sub);

r_res = cellfun(@(x,y) x - y, pred_sub.fixedtau.r, exp_sub.r,'uniformoutput',false); 
th_res = cellfun(@(x,y) sign(y).*(x - y), pred_sub.fixedtau.th, exp_sub.th,'uniformoutput',false);
[rho.fixedtau.r, pval.fixedtau.r] = cellfun(@(x,y) corr(x(:),y(:)),tau,r_res); % make sure they're columns
[rho.fixedtau.th, pval.fixedtau.th] = cellfun(@(x,y) corr(x(:),y(:)),tau,th_res);


% make table
for s = 1:Nstim
    model_based_corr_table.distance.actual(:,s) = table(rho.data.r(:,s));
    model_based_corr_table.angle.actual(:,s) = table(rho.data.th(:,s));
    model_based_corr_table.distance.static(:,s) = table(rho.static.r(:,s));
    model_based_corr_table.angle.static(:,s) = table(rho.static.th(:,s));
    model_based_corr_table.distance.dynamic(:,s) = table(rho.dynamic.r(:,s));
    model_based_corr_table.angle.dynamic(:,s) = table(rho.dynamic.th(:,s));
    model_based_corr_table.distance.fixedtau(:,s) = table(rho.fixedtau.r(:,s));
    model_based_corr_table.angle.fixedtau(:,s) = table(rho.fixedtau.th(:,s));
end
model_based_corr_table.distance.actual.Properties.VariableNames = legend_input;
model_based_corr_table.angle.actual.Properties.VariableNames = legend_input;
model_based_corr_table.distance.static.Properties.VariableNames = legend_input;
model_based_corr_table.angle.static.Properties.VariableNames = legend_input;
model_based_corr_table.distance.dynamic.Properties.VariableNames = legend_input;
model_based_corr_table.angle.dynamic.Properties.VariableNames = legend_input;
model_based_corr_table.distance.fixedtau.Properties.VariableNames = legend_input;
model_based_corr_table.angle.fixedtau.Properties.VariableNames = legend_input;
model_based_corr_table.readme = {'Figure 4B,6','Residual errors vs tau Correlations from actual data, Static Bayesian, Dynamic Bayesian, Fixed tau estimate models'};

save('model_based_corr_table','model_based_corr_table');

%% Real vs Model estimated tau (Figure 5A)
clear real_vs_estimated_tau_table

varname = 'fitprs188.mat'; % varname = FitPrs;
finalprs = importdata(varname);
[tau,JS,sub,start_pos,exp_sub,realindx] = extract_data_for_model(subject,'stimtype');
sampling_distribution.sig = (subject(1).trials(1).prs.sig_phi);     sampling_distribution.mu = (subject(1).trials(1).prs.mu_phi);

for i = 1:size(realindx,1)
    for s = 1:size(realindx,2)
        pars = [finalprs.prs.sd_tau(i,s) finalprs.prs.sd_prior(i,s) finalprs.prs.mu_prior(i,s)];
        tau_estimate{i,s} = MAPtau4staticPrior(tau{i,s},pars,sampling_distribution,log_space);
    end
end

for s = 1:Nstim
    real_vs_estimated_tau_table.real(:,s) = table(tau(:,s));
    real_vs_estimated_tau_table.estimated(:,s) = table(tau_estimate(:,s));
end
real_vs_estimated_tau_table.real.Properties.VariableNames = legend_input;
real_vs_estimated_tau_table.estimated.Properties.VariableNames = legend_input;
real_vs_estimated_tau_table.readme = {'Figure 5A'};

save('real_vs_estimated_tau_table','real_vs_estimated_tau_table');

%% Model Parameters (Figure 5B)
clear model_parameters_table

varname = 'fitprs188.mat'; % varname = FitPrs;
finalprs = importdata(varname);
sampling_distribution.sig = (subject(1).trials(1).prs.sig_phi);     sampling_distribution.mu = (subject(1).trials(1).prs.mu_phi);

model_parameters_table.sampling_distribution_mu = sampling_distribution.mu;
model_parameters_table.sampling_distribution_sig = sampling_distribution.sig;

for s = 1:Nstim
    model_parameters_table.likelihood_sig(:,s) = table(finalprs.prs.sd_tau(:,s));
    model_parameters_table.prior_sig(:,s) = table(finalprs.prs.sd_prior(:,s));
    model_parameters_table.prior_mu(:,s) = table(finalprs.prs.mu_prior(:,s));
end
model_parameters_table.likelihood_sig.Properties.VariableNames = legend_input;
model_parameters_table.prior_sig.Properties.VariableNames = legend_input;
model_parameters_table.prior_mu.Properties.VariableNames = legend_input;
model_parameters_table.readme = {'Figure 5B',};

save('model_parameters_table','model_parameters_table');
