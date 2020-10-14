%% Effect of Sensory Modality and Control Dynamics on Human Path Intergration (Stavropoulos et al. 2020)
%% ++++++++++++++++++++++++++++++++ANALYSIS++++++++++++++++++++++++++++++++
%% Choose files
clearall = input('Clear? ');
if clearall
    clear
end

data_folder = 'Z:\Data\Human Data\Acceleration Control Humans\';
cd(data_folder);
folders = dir(data_folder);   folders = folders([folders.isdir]);

indx = cell2mat(arrayfun(@(x) ~strcmp(x.name(1),'.'),folders,'uniformoutput',false));
subject_name = {folders(indx).name};
%% Load data
tic;
[subject_backup,~] = LoadData(data_folder,subject_name);
toc;
%% Methods figures
% Suppl. Fig. S1A,B
bangbang_dynamics;

% Suppl. Fig. S1C
ntrls = 1000; 
tautau = subject(1).trials(1).prs.tau_tau;
taumin = subject(1).trials(1).prs.min_tau;
taumax = subject(1).trials(1).prs.max_tau;
[taulist,Prob,taurange] = GenerateContinuousTaus(ntrls,tautau,taumin,taumax);

% Suppl. Fig. S5
TraveledDistanceForGivenControlAsFunctionOfTau; 

%% R^2 with and without intercept
params = 'stimtype';
[R2, DR2, F] = R2_bias_intercept(subject_backup,params);
% intercept doesn't increase the variance of data explained by the multiplicative model

%% Bias Detection and Multiplicative Model Fit (+ compare 1st-button-push bias with end-of-trial bias)
polyorder = 1;  plt = 0;    intercept = 0;
params = 'stimtype';
[poolindx,legend_input] = get_poolindx(subject,params); 
[bias.r,bias.th] = ErrScatterFit(subject_backup,params,polyorder,intercept,plt);

%% Plot angular vs radial bias
count = 1;
% plot biases
colr = brewermap(size(bias.r,2),'Dark2');

figure; hold on;    axis equal;
xlim([0 1]);ylim([0 1]);    hline(1,'k--'); vline(1,'k--');
plot(0:2,0:2,'k','HandleVisibility','off'); xlim([0 1.5]);  ylim([0 1.5]);

conf_int = .68;
for s = 1:size(bias.r,2)
    [Xellipse,Yellipse,x0,y0] = error_ellipse(bias.r(:,s),bias.th(:,s),conf_int,0);
    
    plot(bias.r(:,s),bias.th(:,s),'x','Color',colr(s,:),'MarkerSize',10);
    plot(Xellipse + x0, Yellipse + y0,'Color',colr(s,:),'HandleVisibility','off');
    plot(x0, y0,'d','Color',colr(s,:),'MarkerFaceColor',colr(s,:),'MarkerSize',15,'HandleVisibility','off');
end
xlabel('radial bias');ylabel('angular bias'); 

% F-test
regress_ftest;

%% Radial vs angular bias for different TAU BINS for each modality
params = 'stimtau';
[mu,sem,~] = bias_taugroups(subject,params,0);

[~,legend_input] = get_poolindx(subject,params); 
BiasTbl = BiasTable(mu,sem,legend_input);

%% residual error (epsilon) vs tau correlation
plt = 1;    enclose_data = 1;
params = 'stimtype';
rho = []; pval = [];
[rho.r,rho.th,pval.r,pval.th,Rsq_r,Rsq_th] = ResErr_TauCorr(subject,params,plt,enclose_data);

[CorrTbl,~] = CorrTable(rho,pval,legend_input);

% Plot corr(r_epsilon,tau) vs corr(th_epsilon,tau)
colr = brewermap(size(rho.r,2),'Dark2');
[~,legend_input] = get_poolindx(subject,params);
figure;    hold on;    axis equal; xlim([0 1]);  ylim([0 1]);

conf_int = .68;
for s = 1:size(rho.r,2)
    [Xellipse,Yellipse,x0,y0] = error_ellipse(rho.r(:,s),rho.th(:,s),conf_int,0);

    plot(0:2,0:2,'k','HandleVisibility','off');
    plot(Xellipse + x0, Yellipse + y0,'Color',colr(s,:),'HandleVisibility','off');
    plot(rho.r(:,s),rho.th(:,s),'o','Color',colr(s,:),'MarkerFaceColor',colr(s,:));
    xlabel('Radial corr(\epsilon_r,\tau)');ylabel('Angular corr(\epsilon_\theta,\tau)');
end
legend(legend_input{:});    

%% Partial correlations
full_model = 1;         int = 0;
params = 'stimtype';    plt = 0;
[poolindx,legend_input] = get_poolindx(subject,params); 
[b_r,stats_r,b_th,stats_th,parcor] = MultiRegrFit(subject_backup,params,full_model,int,plt);

plot_parcor(parcor,full_model,legend_input);

%% Linear regression of responses with tau-dependent terms

[b,h,p] = TauDependentLinearRegressionOfResponses(subject);

%% MODELLING TAU EFFECT ON SUBJECTS RESPONSES
%% TAU ESTIMATION MODEL - Bayesian approach
%% Static prior
params = 'stimtype';
% format data appropriately
[tau,JS,~,start_pos,exp_sub,realindx] = extract_data_for_model(subject,params);
% get Tau Sampling Distribution
sampling_distribution.sig = (subject(1).trials(1).prs.sig_phi);
sampling_distribution.mu = (subject(1).trials(1).prs.mu_phi);

log_space = 1; % fit in logspace
badsfit = 1; % fit with BADS (by Luigi Acerbi)

[FitPrs,MSE] = StaticPriorBayesianModel(tau,JS,start_pos,exp_sub,sampling_distribution,realindx,log_space,badsfit);

%% plot real vs estimated Tau
varname = 'fitprs188.mat'; % varname = FitPrs;
finalprs = importdata(varname);

tau = [];   MAP_tau = [];
for i = 1:size(realindx,1)
    for s = 1:size(realindx,2)
        pars = [finalprs.prs.sd_tau(i,s) finalprs.prs.sd_prior(i,s) finalprs.prs.mu_prior(i,s)];
        tau = .7:.1:6;
        tau_estimate{i,s} = MAPtau4staticPrior(tau,pars,sampling_distribution,log_space);
    end
end
% Plot
plt = 1;
if plt;  figure;  hold on;  axis equal;  xlim([0 6]);  ylim([0 6]); end
RealVsEstimatedTau(tau,tau_estimate,plt);

%% Look at fitted parameters
colr = brewermap(size(realindx,2),'Dark2');
c=[]; d=[]; e=[]; f=[]; C = [];
% get log-space parameters
C.sd_tau = finalprs.prs.sd_tau;         sigL_title = '\sigma likelihood [log sec]';
C.sd_prior = finalprs.prs.sd_prior;     sigP_title = '\sigma prior [log sec]';
C.mu_prior = finalprs.prs.mu_prior;     muP_title = '\mu prior [log sec]';
RealPrior.sig = sampling_distribution.sig;
RealPrior.mu = sampling_distribution.mu;
lb = [0 0 -.5]; ub = [1 1 2];

% average
figure;
for s = 1:size(realindx,2)
    subplot(3,2,1:2); hold on;
    c = C.sd_tau(:,s); ylim([lb(1) ub(1)]);
    h = bar(s,mean(c),'FaceColor',colr(s,:),'edgecolor','none'); ylabel(sigL_title);
    errorbar(s,mean(c),std(c)/sqrt(length(subject)),'k','capsize',0,'linestyle','none'); title(sigL_title);
    xticks([1 2 3]);   xticklabels({'vestibular','visual','combined'});
    
    subplot(3,2,3:4); hold on;
    c = C.sd_prior(:,s);  ylim([lb(2) ub(2)]);
    h = bar(s,mean(c),'FaceColor',colr(s,:),'edgecolor','none'); ylabel(sigP_title);
    errorbar(s,mean(c),std(c)/sqrt(length(subject)),'k','capsize',0,'linestyle','none'); title(sigP_title);
    xticks([1 2 3]);   xticklabels({'vestibular','visual','combined'}); hline(RealPrior.sig,'k--');
    
    subplot(3,2,5:6); hold on;
    c = C.mu_prior(:,s);  ylim([lb(3) ub(3)]);
    h = bar(s,mean(c),'FaceColor',colr(s,:),'edgecolor','none'); ylabel(muP_title);
    errorbar(s,mean(c),std(c)/sqrt(length(subject)),'k','capsize',0,'linestyle','none'); title(muP_title);
    xticks([1 2 3]);   xticklabels({'vestibular','visual','combined'}); hline(RealPrior.mu,'k--');
end
suptitle(varname);

% compare differences in parameters
[p,h] = CompareFittedPrs(C);

% significance of difference of prior parameters from sampling distribution
[p,h] = CompareFittedPrs2SamplingDistr(C,RealPrior);

%% TAU ESTIMATION MODEL - Fixed tau estimate
%% Fixed tau model
params = 'stimtype';

[tau, JS, ~, start_pos, exp_sub, realindx] = extract_data_for_model(subject,params);        

[fixedtau,MSE_fixedtau] = FixedTauEstimateModel(tau,JS,start_pos,exp_sub,realindx);

%% TAU ESTIMATION MODEL - Dynamic Prior model
%% Kalman-like update of prior in the next trial
% Investigate whether prior should be dynamic
log_space = 1;
badsfit = 1;   
params = 'stimtype';
[poolindx,legend_input] = get_poolindx(subject,params);

% Tau Sampling Distribution
sampling_distribution.sig = (subject(1).trials(1).prs.sig_phi);
sampling_distribution.mu = (subject(1).trials(1).prs.mu_phi);

% extract data
[tau,JS,sub_resp,start_pos,exp_sub,realindx] = extract_data_for_model(subject,params);

[FitPrs,MSE] = DynamicPriorBayesianModel(tau,JS,start_pos,exp_sub,sampling_distribution,realindx,log_space,badsfit);

%% Compare correlations
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

CompareTauCorr(rho.data,rho.static); % Actual vs Subjective Correlations

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

% Compare correlations
model_labels = {'Data','Static Prior'};
[p] = CompareTauModels(rho.data, rho.static, model_labels);
suptitle('Static Prior');

model_labels = {'Data','Fixed \tau'};
[p] = CompareTauModels(rho.data, rho.fixedtau, model_labels);
suptitle('Fixed \tau');

model_labels = {'Data','Dynamic prior'};
[p] = CompareTauModels(rho.data, rho.dynamic, model_labels);
suptitle('Dynamic prior');

model_labels = {'Static Prior','Fixed \tau'};
[p] = CompareTauModels(rho.static, rho.fixedtau, model_labels);
suptitle('Static Prior vs Fixed \tau');

model_labels = {'Static Prior','Dynamic prior'};
[p] = CompareTauModels(rho.static, rho.dynamic, model_labels);
suptitle('Static Prior vs Dynamic prior');

% differences of predicted correlations from zero
a = 1-(1-0.05)^(1/Nstim); % adjustment for multiple comparisons
Nsubs = length(subject); Nstim = size(realindx,2);
Nboots = 1000000;
bstrpmat = randi(Nsubs,Nsubs,Nboots);
% distance
tmp = mat2cell(rho.static.r, Nsubs, ones(Nstim,1));
[p,h] = cellfun(@(x) signrank(x), tmp);
p < a
tmp = cellfun(@(x) mean(x(bstrpmat)), tmp,'uniformoutput',false);
pval = cellfun(@(x) sum(x < 0)/Nboots, tmp);     p_r = min([pval ; 1-pval],[],1);
% angle
tmp = mat2cell(rho.static.th, Nsubs, ones(Nstim,1));
[p,h] = cellfun(@(x) signrank(x), tmp);
p < a
tmp = cellfun(@(x) mean(x(bstrpmat)), tmp,'uniformoutput',false);
pval = cellfun(@(x) sum(x < 0)/Nboots, tmp);     p_th = min([pval ; 1-pval],[],1);

% differences of updating prior weights
p_weight = ComparePriorWeightsOfUpdatingPrior(dynamic_prior);

%% Check correlation between tau and travel duration or mean travel velocity (Suppl. Fig. S4A)

[rho,pval] = trl_duration_vs_tau(subject,params);

%% Simulate uncertainty growth for different models of uncertainty scaling (Suppl. Fig. S4B)
interc = 100;
power_par = [0.5 1 2 2.5];
modelname = {'sublinear','linear','quadratic','supraquadratic'};

clear regr r p rho pval slope
for i = 1:numel(power_par)
    [regr.(modelname{i}),r.(modelname{i}),p.(modelname{i})] = ...
        simulate_uncertainty_growth(subject, power_par(i), interc, plt);
end

slope.mu = structfun(@(x) nanmean(real(x)),regr,'uniformoutput',false);
slope.se = structfun(@(x) nanstd(real(x))./sqrt(length(subject)),regr,'uniformoutput',false);
rho.mu = structfun(@(x) nanmean(real(x)),r,'uniformoutput',false);
rho.se = structfun(@(x) nanstd(real(x))./sqrt(length(subject)),r,'uniformoutput',false);
pval = structfun(@(x) sum(x < 0.05),p,'uniformoutput',false);

tempy = cell2mat(struct2cell(rho.mu));
temperr = cell2mat(struct2cell(rho.se));

colr = brewermap(size(poolindx,2),'Dark2');
figure; hold on;
sz = size(poolindx,2);
sp = numel(power_par);
for s = 1:size(poolindx,2)
    
    for i = 1:numel(power_par)
    bar(sz*(i-1)*2+s,tempy(i,s),'facecolor',colr(s,:),'edgecolor','none'); 
    errorbar(sz*(i-1)*2+s,tempy(i,s),temperr(i,s),'k','capsize',0);
    end
end
axis([0 sz*(1+2*(sp-1))+1 0 0.6]); xticks([sz*((1:sp)-1)*2+2]); xticklabels(modelname); ylabel('corr. coefficient');
title('Accumulated uncertainty over tau');
