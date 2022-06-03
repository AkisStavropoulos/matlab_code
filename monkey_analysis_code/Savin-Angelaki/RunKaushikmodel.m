%% Run the analysis
clc
clear

for ss=  121 %29 %(Jan 22)
   %DONE [30,31,32,33,34,35]%[83,84,86,87,88,89]%[90:93,95:100]%[104:117,120:134,136]%
   %Trial mismatch [112,118,125,131(no-lfp)]
    clear experiments
    experiments = experiment('firefly-monkey');
    experiments.AddSessions(53,ss,{'behv','units','lfps'});
%     experiments.AddSessions(53,ss,{'behv','units','lfps','pop'});
    
    filname=['m',num2str(experiments.sessions.monk_id),'s',num2str(experiments.sessions.sess_id),'_GAM.mat'];
    save(filname,'experiments','-v7.3');
end


