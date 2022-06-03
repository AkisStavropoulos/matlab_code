%% Script for extracting the data pre-GAM fit
% 
%cd 'C:\Users\eb162\Firefly_sina\firefly-monkey'
cd '/Users/edoardo/Work/Code/Firefly_sina/firefly-monkey'
% Save the for the repo
if contains(computer,'MAC')
    separ = '/';
else
    separ = '\';
end
this_path = pwd;
path_split = split(this_path,separ);
base_fold = join(path_split(1:end-1),separ);
if ~contains(path,base_fold{1})
    addpath(fullfile(base_fold{1},'genpath2'))
    addpath(genpath2(base_fold{1},{'.git','genpath2'}))
    
end

w = warning ('off','all');

%% load and save the monkey info
monkeyInfoFile_joysticktask;
save('monkey_info.mat','monkeyInfo')

%% extract a file

monk_id = 51;
sess_list = [38,40,41, 42, 43, 44,45, 46, 47, 123, 124, 125, 15,...
   16, 17,18,19,20, 21, 22,23 ,24,25,26,27,28,29,30,31,326];%2:26;


not_done = [];
except_struct = struct();
for session_id = sess_list
        experiments = experiment('firefly-monkey');
        
        try
            prs = default_prs(monk_id,session_id);
            prs.compute_coherencyLFP = 1;
            prs.analyse_alpha = 0;
            prs.analyse_beta = 0;
            prs.analyse_theta = 0;
            experiments.AddSessions(monk_id,session_id,{'behv','units','lfps'},prs);
        catch ME
            str_id = sprintf('m%ds%d',monk_id,session_id);
            except_struct.(str_id) = ME;
            not_done = [not_done,session_id];
            save([this_path,separ,'not_done.mat'],'not_done','except_struct')
        end
        disp('...clearing exp after session extraction')
        clear experiments
end
