%% export tracking index
clear all
clc
addpath(genpath('/Users/edoardo/Work/Code/Firefly_sina/'));

%folder = '/Volumes/WD Edo/firefly_analysis/LFP_band/DATASET/PPC+MST/';
folder = '/Users/jean-paulnoel/Documents/Savin-Angelaki/For_behavior/Data/m44/';


% get session done
listing = dir(folder);
session_done = {};
% for ii = 1:length(listing)
%     if ~endsWith(listing(ii).name,'eyeTrack.mat')
%         continue
%     end
%     splt = strsplit(listing(ii).name,'_eyeTrack.mat');
%     session_done{ii} = splt{1};
% end
% 
% % cycle over not done
exception_save = {};
for ii = 1:length(listing)
    if listing(ii).name(1) ~= 'm'
        continue
    end
    if endsWith(listing(ii).name,'eyeTrack.mat')
        continue
    end
    if ~endsWith(listing(ii).name,'.mat')
        continue
    end
    
    
    sprintf('loading %s...',listing(ii).name)
    splt = strsplit(listing(ii).name,'.');
    session = splt{1};
    
    done_flag = false;
    for jj = 1:length(session_done)
        if strcmp(session, session_done{jj})
            done_flag = true;
        end
    end
    
    if done_flag
        continue
    end
    
    new_name = strcat(splt{1},'_eyeTrack.mat');
    
    load(strcat(folder,listing(ii).name))

    % reproduce the beginning of Analyse Behaviour

    prs.regress_behv = true;
    prs.extractonly = false;
    prs.regress_eye = true;
    prs.compute_linreg = false;

    %cd('/Users/edoardo/Work/Code/Firefly_sina/firefly-monkey/Analysis/methods__analyse')
    cd('/Users/jean-paulnoel/Documents/Savin-Angelaki/Firefly_sina/firefly-monkey/Analysis/methods__analyse')
    
%       stats = align_eye_tracking_with_behavior(trials_behv,prs);
%         behav_stats = stats;
%         save(strcat(folder,new_name),'prs','trials_behv','units','behav_stats')
    try
        stats = align_eye_tracking_with_behavior(trials_behv,prs);
        behav_stats = stats;
        save(strcat(folder,new_name),'prs','trials_behv','units','behav_stats')
    catch ME
        exception_save{ii} = ME;
        session_save{ii} = session;
        save('exception.mat','exception_save','session_save')
    end
    
end

% check_temp_align = true;
% cnt = 1;
% for tr = 1:length(trials_behv)
%     if stats.trialtype.all.trlindx(tr)
%         ts_beh = trials_behv(tr).continuous.ts;
%         ts_eyemov = stats.trialtype.all.eye_movement.ts{cnt};
%         
%         % take the eye true position
%         ver_mean = stats.trialtype.all.eye_movement.eyepos.true.ver_mean.val{cnt};
%         hor_mean = stats.trialtype.all.eye_movement.eyepos.true.hor_mean.val{cnt};
%         ver_diff = stats.trialtype.all.eye_movement.eyepos.true.ver_diff.val{cnt};
%         hor_diff = stats.trialtype.all.eye_movement.eyepos.true.hor_diff.val{cnt};
%         
%         % take the eye position predicted
%         ver_mean_pred = stats.trialtype.all.eye_movement.eyepos.pred.ver_mean.val{cnt};
%         hor_mean_pred = stats.trialtype.all.eye_movement.eyepos.pred.hor_mean.val{cnt};
%         ver_diff_pred = stats.trialtype.all.eye_movement.eyepos.pred.ver_diff.val{cnt};
%         hor_diff_pred = stats.trialtype.all.eye_movement.eyepos.pred.hor_diff.val{cnt};
%         
%         
%         % get indices
%         idx = zeros(length(ts_eyemov),1);
%         ii = 1;
%         for time = ts_eyemov'
%             idx(ii) = find(ts_beh == time);
%             ii = ii + 1;
%         end
%         % check time alignment *optional*
%         if check_temp_align
%             hor_mean_beh = 0.5*(trials_behv(tr).continuous.yle + trials_behv(tr).continuous.yre);
%             non_nan = ~isnan(hor_mean);
%             idx_non_nan = idx(non_nan);
%             assert(all(hor_mean_beh(idx_non_nan) == hor_mean(non_nan)))         
%         end
%         
%         % create the new vectors of the same length of the behav trial
%         ver_mean_new = nan(length(ts_beh),1);
%         hor_mean_new = nan(length(ts_beh),1);
%         ver_diff_new = nan(length(ts_beh),1);
%         hor_diff_new = nan(length(ts_beh),1);
%         
%         ver_mean_pred_new = nan(length(ts_beh),1);
%         hor_mean_pred_new = nan(length(ts_beh),1);
%         ver_diff_pred_new = nan(length(ts_beh),1);
%         hor_diff_pred_new = nan(length(ts_beh),1);
%         
%         % fill the values
%         ver_mean_new(idx) = ver_mean;
%         hor_mean_new(idx) = hor_mean;
%         ver_diff_new(idx) = ver_diff;
%         hor_diff_new(idx) = hor_diff;
%         
%         ver_mean_pred_new(idx) = ver_mean_pred;
%         hor_mean_pred_new(idx) = hor_mean_pred;
%         ver_diff_pred_new(idx) = ver_diff_pred;
%         hor_diff_pred_new(idx) = hor_diff_pred;
%         
%         % replace values
%         stats.trialtype.all.ts{cnt} = ts_beh; 
%         stats.trialtype.all.eye_movement.eyepos.true.ver_mean.val{cnt} = ver_mean_new;
%         stats.trialtype.all.eye_movement.eyepos.true.hor_mean.val{cnt} = hor_mean_new;
%         stats.trialtype.all.eye_movement.eyepos.true.ver_diff.val{cnt} = ver_diff_new;
%         stats.trialtype.all.eye_movement.eyepos.true.hor_diff.val{cnt} = hor_diff_new;
%         
%         stats.trialtype.all.eye_movement.eyepos.pred.ver_mean.val{cnt} = ver_mean_pred_new;
%         stats.trialtype.all.eye_movement.eyepos.pred.hor_mean.val{cnt} = hor_mean_pred_new;
%         stats.trialtype.all.eye_movement.eyepos.pred.ver_diff.val{cnt} = ver_diff_pred_new;
%         stats.trialtype.all.eye_movement.eyepos.pred.hor_diff.val{cnt} = hor_diff_pred_new;
%              
%         cnt = cnt + 1;
%     end
% end
