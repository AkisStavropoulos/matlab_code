%% reformat data
clear all
base_file = '/Volumes/WD Edo/firefly_analysis/LFP_band/DATASET/MST/m53s127.mat';
load(base_file)



%% reformat
lfp_beta = struct();
lfp_theta = struct();
lfp_alpha = struct();
lfps_new = struct();
for ch = 1: length(lfps)
    lfp_beta(ch).trials = struct();
    lfp_theta(ch).trials = struct();
    lfp_alpha(ch).trials = struct();
    lfps_new(ch).trials = struct('lfp',[]);
    lfps_new(ch).channel_id = lfps(ch).channel_id;
    lfps_new(ch).electrode_id = lfps(ch).electrode_id;
    lfps_new(ch).brain_area = lfps(ch).brain_area;
    for tr = 1:length(lfps(ch).trials)
        lfp_beta(ch).trials(tr).lfp_beta = lfps(ch).trials(tr).lfp_beta;
        lfp_theta(ch).trials(tr).lfp_beta = lfps(ch).trials(tr).lfp_theta;
        lfp_alpha(ch).trials(tr).lfp_alpha = lfps(ch).trials(tr).lfp_alpha;
        lfps_new(ch).trials(tr) = rmfield(lfps(ch).trials(tr),{'lfp_beta','lfp_alpha','lfp_theta'});
        
    end
    
end
 %% clear lfps and save
clear lfps
lfps = lfps_new;
save('/Volumes/WD Edo/firefly_analysis/LFP_band/DATASET/MST/m53s127_new.mat','behv_stats','lfps','prs','trials_behv','units')
save('/Volumes/WD Edo/firefly_analysis/LFP_band/DATASET/MST/lfp_beta.mat','lfp_beta')
save('/Volumes/WD Edo/firefly_analysis/LFP_band/DATASET/MST/lfp_alpha.mat','lfp_alpha')
save('/Volumes/WD Edo/firefly_analysis/LFP_band/DATASET/MST/lfp_theta.mat','lfp_theta')