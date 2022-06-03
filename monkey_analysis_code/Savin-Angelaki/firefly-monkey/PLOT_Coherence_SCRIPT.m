%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Jean-Paul Noel - Plot Coherence Phase stuff for traveling wave
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Housekeeping
close all; clear all; clc; 

% Saving location
save_folder = 'D:\Savin-Angelaki\Figures Coherence Phase';
data_folder = 'D:\Savin-Angelaki\saved';

%% Session of interest
session = 'm44s187'; 

%%
load(fullfile(data_folder, [session, '.mat'])); % lfps and prs are here
load(fullfile(data_folder, ['LFP_coherence_', session, '.mat'])) % stats LFP is here
pop_lfps.stats = stats_lfp; 

for i = 1:size(lfps, 2)
    if lfps(i).brain_area == 'MST' | lfps(i).brain_area == 'VIP'
        lfps(i).electrode_type = 'linearprobe24'; 
    else
       %lfps(i).electrode_type = 'linearprobe24'; 
       lfps(i).electrode_type = 'utah96';  
    end
end

plotCoherence(lfps,pop_lfps,prs,44,187,save_folder)




