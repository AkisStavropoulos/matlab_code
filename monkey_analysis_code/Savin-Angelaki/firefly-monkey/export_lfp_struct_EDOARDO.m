%% Script for extracting the data pre-GAM fit
% 
%% load and save the monkey info
cd('Z:\Users\Soren\Savin-Angelaki\firefly-monkey\Analysis') % read monkeyInfo from Soren's folder on the server
monkeyInfoFile_joysticktask; % MAKE SURE YOU EDIT THE VERSION ON THE SERVER!!!!!!
 
cd('G:\My Drive\MATLAB\Code\Savin-Angelaki\firefly-monkey')

save('monkey_info.mat','monkeyInfo')

% Save the for the repo
separ = '\';
this_path = pwd;
path_split = split(this_path,separ);
base_fold = join(path_split(1:end-1),separ);
if ~contains(path,base_fold{1})
    addpath(fullfile(base_fold{1},'genpath2'))
    addpath(genpath2(base_fold{1},{'.git','genpath2'}))
    
end
w = warning ('off','all');


%% extract a file
% monk_id = [71];
% sess_list = [1003 1004 1005]; % [12 13 14 15 16 18 22 23 24 25 26 27]
% monk_id = [71]; sess_list = [480 484 466 488 490]; 
% viktor_unity = [1003 1004 1005 1006 1007 1008 1009 1011 1012 1013]
% viktor_spike2 = [468 480 484 466 488 490]
monk_id = [73]; sess_list = [34 35 36 37 38 39]; 

not_done = [];
except_struct = struct();
for session_id = sess_list
        experiments = experiment('firefly-monkey');
        prs = default_prs(monk_id,session_id);
        try
            experiments.AddSessions_EDOARDO(monk_id,session_id,{'behv','units','lfps'});
%             experiments.AddSessions_EDOARDO(monk_id,session_id,{'behv'});
        catch ME
            str_id = sprintf('m%ds%d',monk_id,session_id);
            except_struct.(str_id) = ME;
            not_done = [not_done,session_id];
            save([this_path,separ,'not_done.mat'],'not_done','except_struct')
        end
        disp('...clearing exp after session extraction')
        clear experiments
end
