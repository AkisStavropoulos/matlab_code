function [trl,mc] = AddMCData(data,t)
% organize the Motion Cueing variables that are saved from moogdots
% sampled at 60 Hz !!!
% Note that VR position is NOT subject position. It's the position of the SCREEN
% For SUBJECT's position, use SMR data
% If only one input provided , generates only MC variables
%% Put everything in a struct and extract variable names
if ~isempty(data)
    mc.timestamp = data(:,1);
    mc.timestamp = (mc.timestamp - mc.timestamp(1))*.001;
    mc.flag = data(:,2);
    % non-effective button push recordings
    flag2 = mc.flag;
    flag2(flag2 == 1) = 0; % = 2 as long as the button is pressed down
    flag2(flag2 == 2) = 1; % set = 1 to find the timepoints
    
    if isempty(find(flag2))
        flag2 = nan;
    end
    
    mc.flag2 = flag2;
    if ~isempty(mc.flag2)
        mc.flag2 = mc.flag2.*mc.timestamp;
    end
    % actual button push
    mc.flag(mc.flag == 2) = 0;
    mc.flag = mc.flag.*mc.timestamp;
    
    mc.JS_X_Raw = data(:,3); % forward JS input scaled by vmax, looks converted
    mc.JS_Yaw_Raw = data(:,4); % rotational JS input scaled by wmax, looks converted
%     mc.JS_X_final = -data(:,5); % WTF is this?? subject velocity? exactly the same as (col 7 * 100)
%     mc.JS_Yaw_final = data(:,6); % WTF is this??? looks like angular velocity, exactly the same as col 9,10
    
    % mc.filteredBall_X_Vel = data(:,7); % forward
    % mc.filteredBall_Y_Vel = data(:,8); % probably 0 because there is no sideways translation
    % mc.filteredBall_Yaw_Vel = data(:,9); % angular velocity (w)?
    % check 7,8,9
  %%  
%     mc.VR_Yaw_Vel = data(:,10); % angular velocity (w)? same as filteredBall_Yaw_Vel = data(:,5)??
%     mc.VR_Yaw_Pos = -data(:,11);
%     
%     mc.VR_X_Vel = data(:,12); % not final, forward
%     mc.VR_Y_Vel = data(:,13); % not final
%     
%     mc.VR_X_Acc = data(:,14); % forward
%     mc.VR_Y_Acc = data(:,15);
%     
%     mc.Rat_X_Acc = data(:,16); % what is this??
%     mc.Rat_Y_Acc = data(:,17); % what is this??
%     
%     mc.desiredMoog_X_Acc = data(:,18); % forward
%     mc.desiredMoog_Y_Acc = data(:,19);
   %%
    % mc.MC_InternalVariable0_X = data(:,20);
    % mc.MC_InternalVariable1_X = data(:,21);
    % mc.MC_InternalVariable2_X = data(:,22);
    % mc.MC_InternalVariable0_Y = data(:,23);
    % mc.MC_InternalVariable1_Y = data(:,24);
    % mc.MC_InternalVariable2_Y = data(:,25);
    %%
%     mc.MC_X_Acc = data(:,26);
%     mc.MC_Y_Acc = data(:,27);
    %%
    % mc.MC_TiltX_Pos = data(:,28);
    % mc.MC_TiltY_Pos = data(:,29);
    % mc.MC_TiltX_Vel = data(:,30);
    % mc.MC_TiltY_Vel = data(:,31);
    % mc.MC_TiltX_Acc = data(:,32);
    % mc.MC_TiltY_Acc = data(:,33);
    %%
%     mc.Moog_X_Pos = data(:,34);
%     mc.Moog_Y_Pos = data(:,35);
%%
    % mc.Moog_X_Vel = data(:,36);
    % mc.Moog_Y_Vel = data(:,37);
    % mc.Moog_X_Acc = data(:,38);
    % mc.Moog_Y_Acc = data(:,39);
    %
%     mc.Moog_TiltX_Pos = data(:,40);
%     mc.Moog_TiltY_Pos = data(:,41);
    mc.Moog_TiltX_Vel = data(:,42);
    mc.Moog_TiltY_Vel = data(:,43);
    % mc.Moog_TiltX_Acc = data(:,44);
    % mc.Moog_TiltY_Acc = data(:,45);
    %
%     mc.finalMoog_X_Vel = data(:,46);
%     mc.finalMoog_Y_Vel = data(:,47);
%     mc.finalMoog_X_Acc = data(:,48);
%     mc.finalMoog_Y_Acc = data(:,49);
    
    mc.finalMoog_X_GIA = data(:,50);
    mc.finalMoog_Y_GIA = data(:,51);
    
    mc.Rat_X_GIA = data(:,52);
    mc.Rat_Y_GIA = data(:,53);
    
    mc.Rat_X_GIAerror = data(:,54); % what is this??
    mc.Rat_Y_GIAerror = data(:,55);
    
%     mc.VR_X_GIAerror = data(:,56); % what's the relationship with the above error variables?
%     mc.VR_Y_GIAerror = data(:,57);
    
%     mc.finalVR_X_Vel = data(:,58); % forward velocity
%     mc.finalVR_Y_Vel = data(:,59);% sideways velocity
%     
%     mc.VR_X_Pos = 100*data(:,60); % forward position
%     mc.VR_Y_Pos = -100*data(:,61); % sideways position
%     mc.Moog_Yaw_Pos = data(:,62);
    
    % finalVR_X_Vel_Corrected = data(:,59); % corrected VR X position(57)?? in cm
    % finalVR_Y_Vel_Corrected = data(:,60); % corrected VR Y position(56)?? in cm
    
    varnames = fieldnames(mc);
    trl = [];
    % find times of effective button push
    flag_time = unique(mc.flag);
    
    % find times of NON-effective button push
    butt_on = find(diff(flag2) == 1) + 1;
    push_on = flag2(butt_on).*mc.timestamp(butt_on);
    
    butt_off = find(diff(flag2) == -1);
    push_off = flag2(butt_off).*mc.timestamp(butt_off);
    
    mc.push_on = push_on;
    mc.push_off = push_off;
    
    % find when target disappears
    indx_off = find(mc.flag);
    indx_off = indx_off(1:end-1)+1;
    off_time = mc.timestamp(indx_off);
    
    % interpolate missing timepoints
    ts = 0:1/60:mc.timestamp(end);
    for i=2:length(varnames)
        if strcmp(varnames{i},'flag2') && isempty(mc.(varnames{i}))
            continue;
        else
            mc.(varnames{i}) = interp1(mc.timestamp,mc.(varnames{i}),ts,'previous');
        end
    end
    mc.timestamp = interp1(mc.timestamp,mc.timestamp,ts);
    % hold on;plot(mc.timestamp,mc.JS_X_Raw,'.')
    
    for i=1:length(varnames)
        mc.(varnames{i}) = mc.(varnames{i})';
    end
    
    % 
    
    % fix flag (button push)    
    indx = [];
    for i = 1:length(flag_time)
        temp = find(mc.timestamp >= flag_time(i),1);
        indx = [indx temp];
    end
    flag_temp = zeros(length(mc.timestamp),1);
    flag_temp(indx) = 1;
    flag_temp = flag_temp.*mc.timestamp;
    mc.flag = flag_temp;
    % fix target off
    indx_off = [];
    for i = 1:length(off_time)
        temp = find(mc.timestamp >= off_time(i),1);
        indx_off = [indx_off temp];
    end
    target_temp = zeros(length(mc.timestamp),1);
    target_temp(indx_off) = 1;
    target_temp = target_temp.*mc.timestamp;
    mc.target_off = target_temp;
    % could've avoided all that, by using 'next' in interp1 .....
    % get new variable names
    varnames = fieldnames(mc);
    
    if ~isempty(t)
        %% extract trials 
        % dt = dt*prs.factor_downsample;
        for j=1:length(t.end)
            % define pretrial period
            %     pretrial = prs.pretrial; % extract everything from "movement onset - pretrial" or "target onset - pretrial" - whichever is first
            %     posttrial = prs.posttrial; % extract everything until "t_end + posttrial"
            for i=1:length(varnames)
                if ~any(strcmp(varnames{i},[{'push_on'} {'push_off'}])) && ~isempty(mc.(varnames{i}))
                    trl(j).mc.(varnames{i}) = mc.(varnames{i})(mc.timestamp>t.beg(j) & mc.timestamp<t.end(j));
                %             trl(j).continuous.(chnames{i}) = downsample(trl(j).continuous.(chnames{i}),prs.factor_downsample);
                else
                    trl(j).mc.(varnames{i}) = mc.(varnames{i})(mc.(varnames{i}) > t.beg(j) & mc.(varnames{i}) < t.end(j));
                end
            end
            %     trl(j).events.t_move = t.move(j);
            %     trl(j).events.t_stop = t.stop(j);
        end
        
        emptytrl = [];
        for n = 1:length(trl)
            if ~isempty(trl(n).mc.timestamp)
                % start time at t = 0 for every trial
                trl(n).mc.push_on =  trl(n).mc.push_on - trl(n).mc.timestamp(1);
                trl(n).mc.push_off =  trl(n).mc.push_off - trl(n).mc.timestamp(1);
                trl(n).mc.timestamp =  trl(n).mc.timestamp - trl(n).mc.timestamp(1);
            else
                emptytrl = [emptytrl n];
            end
        end
        if emptytrl
            disp(['empty MC trials from trial ' num2str(min(emptytrl)) ' to trial ' num2str(max(emptytrl)) ' !!'])
        end
    end
    
else
    % if MC_Variables doesn't exist or empty
    mc.timestamp = nan;
    mc.flag = nan;
    mc.flag2 = nan;
    
    mc.JS_X_Raw = nan; % forward JS input scaled by vmax, looks converted
    mc.JS_Yaw_Raw = nan; % rotational JS input scaled by wmax, looks converted
    mc.JS_X_final = nan; % WTF is this?? subject velocity? exactly the same as (col 7 * 100)
    mc.JS_Yaw_final = nan; % WTF is this??? looks like angular velocity, exactly the same as col 9,10
    
    % mc.filteredBall_X_Vel = data(:,7); % forward
    % mc.filteredBall_Y_Vel = data(:,8); % probably 0 because there is no sideways translation
    % mc.filteredBall_Yaw_Vel = data(:,9); % angular velocity (w)?
    % check 7,8,9
    mc.VR_Yaw_Vel = nan; % angular velocity (w)? same as filteredBall_Yaw_Vel = data(:,5)??
    mc.VR_Yaw_Pos = nan;
    
    mc.VR_X_Vel = nan; % not final, forward
    mc.VR_Y_Vel = nan; % not final
    
    mc.VR_X_Acc = nan; % forward
    mc.VR_Y_Acc = nan;
    
    mc.Rat_X_Acc = nan; % what is this??
    mc.Rat_Y_Acc = nan; % what is this??
    
    mc.desiredMoog_X_Acc = nan; % forward
    mc.desiredMoog_Y_Acc = nan;
    %
    % mc.MC_InternalVariable0_X = data(:,20);
    % mc.MC_InternalVariable1_X = data(:,21);
    % mc.MC_InternalVariable2_X = data(:,22);
    % mc.MC_InternalVariable0_Y = data(:,23);
    % mc.MC_InternalVariable1_Y = data(:,24);
    % mc.MC_InternalVariable2_Y = data(:,25);
    
    mc.MC_X_Acc = nan;
    mc.MC_Y_Acc = nan;
    %
    % mc.MC_TiltX_Pos = data(:,28);
    % mc.MC_TiltY_Pos = data(:,29);
    % mc.MC_TiltX_Vel = data(:,30);
    % mc.MC_TiltY_Vel = data(:,31);
    % mc.MC_TiltX_Acc = data(:,32);
    % mc.MC_TiltY_Acc = data(:,33);
    %
    mc.Moog_X_Pos = nan;
    mc.Moog_Y_Pos = nan;
    % mc.Moog_X_Vel = data(:,36);
    % mc.Moog_Y_Vel = data(:,37);
    % mc.Moog_X_Acc = data(:,38);
    % mc.Moog_Y_Acc = data(:,39);
    %
    mc.Moog_TiltX_Pos = nan;
    mc.Moog_TiltY_Pos = nan;
    % mc.Moog_TiltX_Vel = data(:,42);
    % mc.Moog_TiltY_Vel = data(:,43);
    % mc.Moog_TiltX_Acc = data(:,44);
    % mc.Moog_TiltY_Acc = data(:,45);
    %
    mc.finalMoog_X_Vel = nan;
    mc.finalMoog_Y_Vel = nan;
    mc.finalMoog_X_Acc = nan;
    mc.finalMoog_Y_Acc = nan;
    
    mc.finalMoog_X_GIA = nan;
    mc.finalMoog_Y_GIA = nan;
    
    mc.Rat_X_GIA = nan;
    mc.Rat_Y_GIA = nan;
    
    mc.Rat_X_GIAerror = nan; % what is this??
    mc.Rat_Y_GIAerror = nan;
    
    mc.VR_X_GIAerror = nan; % what's the relationship with the above error variables?
    mc.VR_Y_GIAerror = nan;
    
    mc.finalVR_X_Vel = nan; % forward velocity
    mc.finalVR_Y_Vel = nan; % sideways velocity
    
    mc.VR_X_Pos = nan; % forward position
    mc.VR_Y_Pos = nan; % sideways position
    mc.Moog_Yaw_Pos = nan;
    
    mc.target_off = nan;
    mc.push_on = nan;
    mc.push_off = nan;
    % finalVR_X_Vel_Corrected = data(:,59); % corrected VR X position(57)?? in cm
    % finalVR_Y_Vel_Corrected = data(:,60); % corrected VR Y position(56)?? in cm
    
    varnames = fieldnames(mc);
    
    for j=1:length(t.end)
        % define pretrial period
        %     pretrial = prs.pretrial; % extract everything from "movement onset - pretrial" or "target onset - pretrial" - whichever is first
        %     posttrial = prs.posttrial; % extract everything until "t_end + posttrial"
        for i=1:length(varnames)
            trl(j).mc.(varnames{i}) = mc.(varnames{i});
            %             trl(j).continuous.(chnames{i}) = downsample(trl(j).continuous.(chnames{i}),prs.factor_downsample);
        end
        %     trl(j).events.t_move = t.move(j);
        %     trl(j).events.t_stop = t.stop(j);
    end
    
    disp('MC file is missing for this block !!')
    
end









