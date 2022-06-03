%% Process noise testing

%% trajectories

figure;
for i = 13:numel(trials) % randi([1 numel(trials)],1,50)
    
    plot(trials(i).continuous.ts, trials(i).continuous.xmp,'b'); hold on;
    plot(trials(i).continuous.ts, trials(i).continuous.ymp,'k'); 
    plot(trials(i).continuous.ts, trials(i).continuous.xfp,'b--'); 
    plot(trials(i).continuous.ts, trials(i).continuous.yfp,'--k'); 
    vline(trials(i).continuous.ts(trials(i).continuous.ts==0),'g'); vline(trials(i).events.t_end); hold off;
    
end

%% velocity

figure;
for i = 1:numel(trials) % randi([1 numel(trials)],1,50)
    
    plot(trials(i).continuous.ts, trials(i).continuous.w,'b:'); hold on;
    plot(trials(i).continuous.ts, trials(i).continuous.v,'k:'); 
    hline(trials(i).prs.wmax,'b--'); 
    hline(trials(i).prs.vmax,'--k'); 
%     plot(trials(i).continuous.ts, trials(i).continuous.w_clean,'b:'); 
%     plot(trials(i).continuous.ts, trials(i).continuous.v_clean,'k:'); 
    plot(trials(i).continuous.ts, trials(i).continuous.w_eta*10,'b--'); 
    plot(trials(i).continuous.ts, trials(i).continuous.v_eta*10,'k--'); 
%     plot(trials(i).continuous.ts, trials(i).continuous.w_ksi*10,'b'); 
%     plot(trials(i).continuous.ts, trials(i).continuous.v_ksi*10,'k'); 
    
    plot(trials(i).continuous.ts, trials(i).continuous.xjs*trials(i).prs.wmax,'b'); 
    plot(trials(i).continuous.ts, trials(i).continuous.yjs*trials(i).prs.vmax,'k'); 
    vline(trials(i).continuous.ts(trials(i).continuous.ts==0),'g'); vline(trials(i).events.t_end); hold off;
    
end


%% joystick

figure;
for i = 1:numel(trials) % randi([1 numel(trials)],1,50)
     
    plot(trials(i).continuous.ts, trials(i).continuous.xjs,'r'); hold on;
    plot(trials(i).continuous.ts, trials(i).continuous.yjs,'k'); hold off;
end

%% noise

figure;
for i = 1:numel(trials) % randi([1 numel(trials)],1,50)
    
    plot(trials(i).continuous.ts, trials(i).continuous.w_eta,'b:'); hold on;
    plot(trials(i).continuous.ts, trials(i).continuous.v_eta,'k:'); 
    plot(trials(i).continuous.ts, trials(i).continuous.w_ksi,'b'); 
    plot(trials(i).continuous.ts, trials(i).continuous.v_ksi,'k'); hold off;
end

%% eye movements

figure;
for i = 13:numel(trials) % randi([1 numel(trials)],1,50)
    
    plot(trials(i).continuous.ts, trials(i).continuous.yle,'b'); hold on;
    plot(trials(i).continuous.ts, trials(i).continuous.zle,'k'); 
    ylim([-60 60]); vline(trials(i).continuous.ts(trials(i).continuous.ts==0),'g'); vline(trials(i).events.t_end); hold off;
end

%% Parameters

% tau
tau = arrayfun(@(x) x.prs.tau, trials);
figure;hist(tau)


