%% Script for extracting the data pre-GAM fit
% 
clc; 
cd 'D:\Savin-Angelaki\Firefly_sina\firefly-monkey'
%cd 'C:\Users\eb162\Firefly_sina\firefly-monkey'
%cd '/Users/edoardo/Work/Code/Firefly_sina/neuroGAM'
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
%cd 'D:\Savin-Angelaki\Firefly_sina\firefly-monkey'

monkeyInfoFile_joysticktask;
save('monkey_info.mat','monkeyInfo')

%% extract a file


monk_id = 71;
session = [855];%2:26;
%load('fit_list.mat')

not_done = [];
except_struct = struct();
fld_save = 'Pre-processing X E';
for kk = 1:length(session)
        
        session_id = session(kk);
        %monk_id = monk_id(kk);
        experiments = experiment('firefly-monkey');
        
        try
            prs = default_prs(monk_id,session_id);
           
            experiments.AddSessions(monk_id,session_id,{'behv','units','lfps'});
        catch ME
            str_id = sprintf('m%ds%d',monk_id,session_id);
            except_struct.(str_id) = ME;
            not_done = [not_done,session_id];
            save([this_path,separ,'not_done.mat'],'not_done','except_struct')
        end
        disp('...clearing exp after session extraction')
        clear experiments
end
