
%% Receive Unity data
clc
clear

% set initial parameters
sac_thresh = 10;
dt = 1/90;
pre_sac = round(0.2/dt);
post_sac = round(0.2/dt);

% set up plots
figure;
subplot(2,2,1); hold on; xlabel('time since stimulation [s]'); ylabel('hor. eye [deg]'); axis([-2 2 -40 40]); vline(0,'r');
subplot(2,2,2); hold on; xlabel('time since stimulation [s]'); ylabel('ver. eye [deg]'); axis([-2 2 -40 40]); vline(0,'r');
subplot(2,2,3); hold on; xlabel('hor. eye [deg]'); ylabel('ver. eye [deg]'); axis([-40 40 -40 40]); vline(0,'k'); hline(0,'k');
subplot(2,2,4); hold on; xlabel('saccade direction'); ylabel('saccade amplitude [deg]'); axis([0 3 0 40]); xticks(1:2); xticklabels({'HOR','VER'});

% connect to Unity
tcpipServer = tcpip('127.0.0.1',55000,'NetworkRole','Server');
set(tcpipServer, 'InputBufferSize', 30000); 
fopen(tcpipServer);

% receive data
reset = 1;
trl = 0; t_start = 0; t_end = 0; t_stim = 0; t_rew = 0;
cont_data = []; mrk = []; mrkdiff = 0; ts = []; yle = []; zle = []; yre = []; zre = [];
while reset
    fprintf(tcpipServer, 'GET /');
    if (get(tcpipServer, 'BytesAvailable') > 0) 
        tcpipServer.BytesAvailable 
        
        DataReceived = fscanf(tcpipServer);
        tmp1 = split(DataReceived,',');
        tmp = cellfun(@(x) str2double(x),tmp1)';
        cont_data  = [cont_data ; tmp];
        
        % set-up data
        mrk = [mrk ; cont_data(end,end)];
        mrkdiff = [mrkdiff ; mrk(end)-mrk(end-1)];
        ts = [ts ; cont_data(end,2)];
        
        LXnorm = cont_data(end,14);  RXnorm = cont_data(end,17);
        LYnorm = cont_data(end,15);  RYnorm = cont_data(end,18);
        LZnorm = cont_data(end,16);  RZnorm = cont_data(end,19);
        
        yle = [yle ; atan2d(LXnorm,LZnorm)];
        zle = [zle ; atan2d(LYnorm,LZnorm)];
        yre = [yre ; atan2d(RXnorm,RZnorm)];
        zre = [zre ; atan2d(RYnorm,RZnorm)];
        
        % extract trials
        if mrkdiff(end) && mrk(end)==2
            
            plt = 1;
            t_start = ts(end);
            t_starti = size(cont_data,1);
            
        elseif mrkdiff(end) && mrk(end)==3
            
            t_end = ts(end);
            t_endi = size(cont_data,1);
            
        elseif mrkdiff(end) && mrk(end)==5
            
            t_stim = ts(end);
            t_stimi = size(cont_data,1);
            
        elseif mrkdiff(end) && mrk(end)==4
            
            t_rew = ts(end);
            t_rewi = size(cont_data,1);
            
        elseif mrkdiff(end) && mrk(end)==17
            
            reset = 0; % exit loop
            
        end
        
        % plot data
        if (t_stim>t_start && t_stim<t_end) && plt
            plt = 0;
            trl = trl+1;
            
            trlindx = t_starti:t_endi;
            ts1 = ts(trlindx) - t_stim;
            ver_eye = 0.5*(zle(trlindx) + zre(trlindx));
            hor_eye = 0.5*(yle(trlindx) + yre(trlindx));
            
            pre_sac_hor(trl) = hor_eye(t_stimi-pre_sac-t_starti);
            pre_sac_ver(trl) = ver_eye(t_stimi-pre_sac-t_starti);
            post_sac_hor(trl) = hor_eye(t_stimi+post_sac-t_starti);
            post_sac_ver(trl) = ver_eye(t_stimi+post_sac-t_starti);
            sac_amp_hor(trl) = post_sac_hor(trl)-pre_sac_hor(trl);
            sac_amp_ver(trl) = post_sac_ver(trl)-pre_sac_ver(trl);
            
            subplot(2,2,1); plot(ts1,hor_eye,'k');
            subplot(2,2,2); plot(ts1,ver_eye,'k');
            subplot(2,2,3); plot(pre_sac_hor(trl),pre_sac_ver(trl),'bo'); plot(post_sac_hor(trl),post_sac_ver(trl),'ro');
            subplot(2,2,4); plot(1,abs(sac_amp_hor(trl)),'or'); plot(2,abs(sac_amp_ver(trl)),'or');
            
        end

    end 
end

% compute statistics
[~,p_hor] = ttest(pre_sac_hor(:),post_sac_hor(:));
[~,p_ver] = ttest(pre_sac_ver(:),post_sac_ver(:));
valid_tmp = sac_amp_hor > sac_thresh | sac_amp_ver > sac_thresh;
if valid_tmp; valid_saccade = 'YES'; else; valid_saccade = 'NO!'; end

fprintf('++++++++++++ Stimulation Protocol ++++++++++++ \n\n');
fprintf('HOR:\n');
fprintf('Saccade magnitude: %0.3f \nt-test: %0.4f \n\n', mean(sac_amp_hor), p_hor);
fprintf('VER:\n');
fprintf('Saccade magnitude: %0.3f \nt-test: %0.4f \n\n', mean(sac_amp_ver), p_ver);
fprintf('Valid saccade (>%ddeg): %s \n', sac_thresh, valid_saccade)


% sanity checks
dt = (diff(cont_data(:,2)));
ts = cont_data(:,2);
mrk = cont_data(:,end);

%% Plot superimposed eye positions

cont_data = importdata('C:\Users\ges6\Downloads\eye_data_Test_05262022_339.txt');
hdr = cont_data.textdata(1,:);
expflag = cont_data.textdata(2:end,1);
expflag = cellfun(@(x) contains(x,'Stim'), expflag);
cont_data = cellfun(@(x) str2double(x), cont_data.textdata(2:end,:));

mrk = cont_data(:,end);
mrkdiff = [0 ; diff(mrk)];
ts = cont_data(:,2);

LXnorm = cont_data(:,14);  RXnorm = cont_data(:,17);
LYnorm = cont_data(:,15);  RYnorm = cont_data(:,18);
LZnorm = cont_data(:,16);  RZnorm = cont_data(:,19);

yle = atan2d(LXnorm,LZnorm);  % double check
zle = atan2d(LYnorm,LZnorm);
yre = atan2d(RXnorm,RZnorm);
zre = atan2d(RYnorm,RZnorm);

figure;plot(yle,zle);hold on;plot(yre,zre)

%%
figure; 
subplot(2,2,1); hold on; xlabel('time since stimulation [s]'); ylabel('hor. eye [deg]'); axis([-2 2 -40 40]); vline(0,'r');
subplot(2,2,2); hold on; xlabel('time since stimulation [s]'); ylabel('ver. eye [deg]'); axis([-2 2 -40 40]); vline(0,'r');
subplot(2,2,3); hold on; xlabel('hor. eye [deg]'); ylabel('ver. eye [deg]'); axis([-40 40 -40 40]); vline(0,'k'); hline(0,'k');
subplot(2,2,4); hold on; xlabel('saccade direction'); ylabel('saccade amplitude [deg]'); axis([0 3 0 40]); xticks(1:2); xticklabels({'HOR','VER'});

sac_thresh = 10;
dt = 1/90;
pre_sac = round(0.2/dt);
post_sac = round(0.2/dt);
reset = 1;
t = 0;
trl = 0;
t_start = 0; t_end = 0; t_stim = 0; t_rew = 0;
while reset    
    
    t = t+1;
    if mrkdiff(t) & mrk(t)==2
        
        plt = 1;
        t_start = ts(t);
        t_starti = t;
        
    elseif mrkdiff(t) & mrk(t)==3
        
        t_end = ts(t);
        t_endi = t;
        
    elseif mrkdiff(t) & mrk(t)==5
        
        t_stim = ts(t);
        t_stimi = t;
        
    elseif mrkdiff(t) & mrk(t)==4
        
        t_rew = ts(t);
        t_rewi = t;
    
    elseif mrkdiff(t) & mrk(t)==17

        reset = 0;
               
    end
    
    if (t_stim>t_start & t_stim<t_end) & plt
        plt = 0;
        trl = trl+1;
        
        trlindx = t_starti:t_endi;
        ts1 = ts(trlindx) - t_stim;
        ver_eye = 0.5*(zle(trlindx) + zre(trlindx));
        hor_eye = 0.5*(yle(trlindx) + yre(trlindx));
        
        pre_sac_hor(trl) = hor_eye(t_stimi-pre_sac-t_starti);
        pre_sac_ver(trl) = ver_eye(t_stimi-pre_sac-t_starti);
        post_sac_hor(trl) = hor_eye(t_stimi+post_sac-t_starti);
        post_sac_ver(trl) = ver_eye(t_stimi+post_sac-t_starti);
        sac_amp_hor(trl) = post_sac_hor(trl)-pre_sac_hor(trl);
        sac_amp_ver(trl) = post_sac_ver(trl)-pre_sac_ver(trl);
        
        subplot(2,2,1); plot(ts1,hor_eye,'k'); 
        subplot(2,2,2); plot(ts1,ver_eye,'k');
        subplot(2,2,3); plot(pre_sac_hor(trl),pre_sac_ver(trl),'bo'); plot(post_sac_hor(trl),post_sac_ver(trl),'ro');
        subplot(2,2,4); plot(1,abs(sac_amp_hor(trl)),'or'); plot(2,abs(sac_amp_ver(trl)),'or');
        
    end
    
end

[~,p_hor] = ttest(pre_sac_hor(:),post_sac_hor(:));
[~,p_ver] = ttest(pre_sac_ver(:),post_sac_ver(:));
valid_tmp = sac_amp_hor > sac_thresh | sac_amp_ver > sac_thresh;
if valid_tmp; valid_saccade = 'YES'; else; valid_saccade = 'NO!'; end

clc
fprintf('++++++++++++ Stimulation Protocol ++++++++++++ \n\n');
fprintf('HOR:\n');
fprintf('Saccade magnitude: %0.3f \nt-test: %0.4f \n\n', mean(sac_amp_hor), p_hor);
fprintf('VER:\n');
fprintf('Saccade magnitude: %0.3f \nt-test: %0.4f \n\n', mean(sac_amp_ver), p_ver);
fprintf('Valid saccade (>%ddeg): %s \n', sac_thresh, valid_saccade)
