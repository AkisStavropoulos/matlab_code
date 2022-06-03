%% Quick perturbation analysis
clear;
session_date = {'Feb 17 2022'}; % 'Feb 09 2022','Feb 15 2022',

data_path = 'Z:\Data\Monkey2_newzdrive\Jimmy\Perturbations';
file_path = '\neural data\Pre-processing X E';

trials = []; gdtrls_indx = [];
for n = 1:numel(session_date)
data_folder = fullfile(data_path,session_date{n},file_path);
cd(data_folder);

tmp = dir(data_folder);
filename = tmp(arrayfun(@(x) ~isempty(regexp(x.name,'m\w*s\w*.mat','once')),tmp));

disp(['......... Loading ' filename.name]);
tmp = load(filename.name);
disp('Loaded.')

gdtrls_indx = [gdtrls_indx  tmp.behv_stats.trialtype.all.trlindx];
for i = 1:numel(tmp.trials_behv); tmp.trials_behv(i).prs.session = filename; end
trials = [trials  tmp.trials_behv];
end
disp(['........Total trials = ' num2str(numel(trials))])
clear tmp
%% Remove bad trials
sf = 1/.006;

max_dur = 7;
min_dist = 40;
max_err_ratio = 0.5; % 0.7
[gdtrls, rewtrls, targON_indx] = SortTrials(trials, max_dur, min_dist, max_err_ratio);

ptbindx = arrayfun(@(x) x.logical.ptb,trials);

trlindx_np = gdtrls & ~ptbindx;
trlindx_p = gdtrls & ptbindx;
%% Targets vs Responses

% unperturbed
[tar_tmp, sub_tmp] = TargetsVsResponses(trials(trlindx_np),1);
suptitle('unperturbed');

% perturbed
[tar_tmp, sub_tmp] = TargetsVsResponses(trials(trlindx_p),1);
suptitle('perturbed');

% simulated
x_targ_np = arrayfun(@(x) x.prs.xfp, trials(trlindx_np));
y_targ_np = arrayfun(@(x) x.prs.yfp, trials(trlindx_np));
x_targ_p = arrayfun(@(x) x.prs.xfp, trials(trlindx_p));
y_targ_p = arrayfun(@(x) x.prs.yfp, trials(trlindx_p));
v_ts_np = arrayfun(@(x) x.continuous.v, trials(trlindx_np),'un',0);
w_ts_np = arrayfun(@(x) x.continuous.w, trials(trlindx_np),'un',0);
t_ptbn = arrayfun(@(x) x.events.t_ptb, trials(trlindx_p));
v_ptb = arrayfun(@(x) x.prs.ptb_linear, trials(trlindx_p));
w_ptb = arrayfun(@(x) x.prs.ptb_angular, trials(trlindx_p));
[sim_resp,sim_targ] = gen_sim_ptb(x_targ_np,y_targ_np,v_ts_np,w_ts_np,x_targ_p,y_targ_p,t_ptbn,v_ptb,w_ptb);

figure; 
subplot(1,2,1); hold on;
plot(sim_targ.r,sim_resp.r,'.'); plot(0:500,0:500,'k--'); axis equal; axis([0 500 0 500]); 
xlabel('target distance [cm]'); ylabel('response [cm]');

subplot(1,2,2); hold on;
plot(sim_targ.theta,sim_resp.theta,'.'); plot(-50:50,-50:50,'k--'); vline(0,'k');hline(0,'k'); axis equal; axis([-50 50 -50 50]); 
xlabel('target angle [deg]'); ylabel('response [deg]');
suptitle('simulated');

%% ROC Curve

x_targ_np = arrayfun(@(x) x.prs.xfp, trials(trlindx_np));
y_targ_np = arrayfun(@(x) x.prs.yfp, trials(trlindx_np));
x_targ_p = arrayfun(@(x) x.prs.xfp, trials(trlindx_p));
y_targ_p = arrayfun(@(x) x.prs.yfp, trials(trlindx_p));

v_ts_np = arrayfun(@(x) x.continuous.v, trials(trlindx_np),'un',0);
w_ts_np = arrayfun(@(x) x.continuous.w, trials(trlindx_np),'un',0);

t_ptbn = arrayfun(@(x) x.events.t_ptb, trials(trlindx_p));
v_ptb = arrayfun(@(x) x.prs.ptb_linear, trials(trlindx_p));
w_ptb = arrayfun(@(x) x.prs.ptb_angular, trials(trlindx_p));

[sim_resp,sim_targ] = gen_sim_ptb(x_targ_np,y_targ_np,v_ts_np,w_ts_np,x_targ_p,y_targ_p,t_ptbn,v_ptb,w_ptb);


maxrewardwin = 400;
npermutations = 100;
t_stop = arrayfun(@(x) numel(x.continuous.ts), trials); % arrayfun(@(x) floor(x.events.t_end*sf), trials);

% unperturbed
[r_targ_np,th_targ_np] = cart2polarY(x_targ_np,y_targ_np);
X_fly = [r_targ_np(:) th_targ_np(:)];
x_subj_np = arrayfun(@(x,e) x.continuous.xmp(e-10), trials(trlindx_np), t_stop(trlindx_np));
y_subj_np = arrayfun(@(x,e) x.continuous.ymp(e-10), trials(trlindx_np), t_stop(trlindx_np));
[r_subj_np,th_subj_np] = cart2polarY(x_subj_np,y_subj_np);
X_monk = [r_subj_np(:) th_subj_np(:)]; 
[rewardwin, pCorrect_noptb, pCorrect_shuffled_mu_noptb] = ComputeROCFirefly(X_fly,X_monk,maxrewardwin,npermutations);

% perturbed
[r_targ_p,th_targ_p] = cart2polarY(x_targ_p,y_targ_p);
X_fly = [r_targ_p(:) th_targ_p(:)];
[r_targ_p,th_targ_p] = cart2polarY(x_targ_p,y_targ_p);
x_subj_p = arrayfun(@(x,e) x.continuous.xmp(e), trials(trlindx_p), t_stop(trlindx_p));
y_subj_p = arrayfun(@(x,e) x.continuous.ymp(e), trials(trlindx_p), t_stop(trlindx_p));
[r_subj_p,th_subj_p] = cart2polarY(x_subj_p,y_subj_p);
X_monk = [r_subj_p(:) th_subj_p(:)]; 
[rewardwin, pCorrect_ptb, pCorrect_shuffled_mu_ptb] = ComputeROCFirefly(X_fly,X_monk,maxrewardwin,npermutations);

% simulated
X_fly = [sim_targ.r(:) sim_targ.theta(:)];
X_monk = [sim_resp.r(:) sim_resp.theta(:)];
[rewardwin, pCorrect_sim, pCorrect_shuffled_mu_sim] = ComputeROCFirefly(X_fly,X_monk,maxrewardwin,npermutations);


% plot 
figure;hold on;
plot(pCorrect_shuffled_mu_noptb,pCorrect_noptb,'b');
plot(pCorrect_shuffled_mu_ptb,pCorrect_ptb,'r');
plot(pCorrect_shuffled_mu_sim,pCorrect_sim,'k');
plot(0:1,0:1,'k--')
xlabel('shuffled'); ylabel('true'); legend({'no ptb','ptb','sim'},'location','southeast'); title('ROC curve');

%% Cumulative distribution of error
max_err = 400;
N = 200;
nx = linspace(0,max_err,N+1);

% unperturbed
err_np = sqrt((x_targ_np-x_subj_np).^2 + (y_targ_np-y_subj_np).^2);
[ny_np,~] =  hist(err_np,nx);
cdf_np = cumsum(ny_np(:))./sum(ny_np);

% perturbed
err_p = sqrt((x_targ_p-x_subj_p).^2 + (y_targ_p-y_subj_p).^2);
[ny_p,~] =  hist(err_p,nx);
cdf_p = cumsum(ny_p(:))./sum(ny_p);

% simulated
[x_targ_sim,y_targ_sim] = polar2cartY(sim_targ.r,sim_targ.theta);
[x_subj_sim,y_subj_sim] = polar2cartY(sim_resp.r,sim_resp.theta);
err_sim = sqrt((x_targ_sim-x_subj_sim).^2 + (y_targ_sim-y_subj_sim).^2);
[ny_sim,~] =  hist(err_sim,nx);
cdf_sim = cumsum(ny_sim(:))./sum(ny_sim);

figure; hold on;
plot(nx,cdf_np,'b');
plot(nx,cdf_p,'r');
plot(nx,cdf_sim,'k');
xlabel('euclideian error [cm]'); ylabel('Cum. Prob.');legend({'no ptb','ptb','sim'},'location','southeast');



%% Sanity check: plot trajectories with perturbation
figure;
for i = 1:numel(trials)
    
    if trials(i).events.t_ptb && ~trials(i).logical.ptb
        disp(['Trial ' num2str(i) ' has inconsistent PTB parameters!!!!!'])
    end
    
    if trials(i).logical.ptb && (trials(i).prs.ptb_linear < -30) && gdtrls(i)
    v = trials(i).continuous.v;
    w = trials(i).continuous.w;
    ts = trials(i).continuous.ts;
    t_ptb = trials(i).events.t_ptb;
    
    plot(ts,v); hold on; vline(t_ptb); vline(t_ptb+1); title(['LinPtb = ' num2str(trials(i).prs.ptb_linear)]); hold off;
    
    end
end
