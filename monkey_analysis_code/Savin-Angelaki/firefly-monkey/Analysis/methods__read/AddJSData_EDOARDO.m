function [trl,mc] = AddJSData_EDOARDO(data,t,prs)
% organize the Motion Cueing variables that are saved from moogdots
% sampled at 120 Hz !!! (Polarized projectors setups) 
% Note that VR position is NOT subject position. It's the position of the SCREEN
% For SUBJECT's position, use SMR data
% If only one input provided , generates only MC variables
%% Put everything in a struct and extract variable names
if ~isempty(data)
    dt = 1/prs.fs_mc;

    mc.tstamp = data(:,1);
    mc.tstamp = (mc.tstamp - mc.tstamp(1)) + dt;
    
    mc.jsy = data(:,2); % forward JS input scaled by vmax, looks converted (JS_X_Raw)
    mc.jsx = data(:,3); % rotational JS input scaled by wmax, looks converted (JS_Yaw_Raw)
    varnames = fieldnames(mc);
    
    startend = diff(mc.tstamp) > 0.08;
    endindx = find(diff(startend) == 1);    begindx = find(diff(startend) == -1)+1;
    t_end = mc.tstamp(endindx);     t_beg = [t.beg(1) ; mc.tstamp(begindx)]; % t_beg misses the first trial
    t_beg = t_beg(1:numel(t_end));
    
    if t.ntrls < numel(t_end) % correct markers in case of 2 smr files for one block (until Jing adds markers in MC variables)
        ntrls1 = t.ntrls;        ntrls2 = numel(t_end) - ntrls1;
        t_correction = t_beg(ntrls1+1) - t.beg(ntrls1+1);        
        tfields = fieldnames(t);
        for n = 1:numel(tfields)
            if ~any(strcmp(tfields{n},{'fstart','fstop','eyeblinks','ntrls'}))
                t.(tfields{n})(ntrls1+1:end) = t.(tfields{n})(ntrls1+1:end) + t_correction;
            end
        end
    end
    
    % interpolate missing timepoints (intertrial interval)
    ts = dt:dt:mc.tstamp(end);
    for i=1:length(varnames)
        if ~any(strcmp(varnames{i},'tstamp'))
            mc.(varnames{i}) = interp1(mc.tstamp,mc.(varnames{i}),ts,'previous');
        end
    end
    mc.tstamp = interp1(mc.tstamp,mc.tstamp,ts);
%     mc = structfun(@(x) x', mc, 'un', 0); % hold on;plot(mc.tstamp,mc.JS_X_Raw,'.')
    
    %% Upsample to match .smr frequency (note: .smr has been downsampled)
    dt_original = mc.tstamp(2)-mc.tstamp(1);
    fs_mc = 1/dt_original;
    fs_smr = prs.fs_smr/prs.factor_downsample;
    factor_upsample = fs_mc/fs_smr;
    dt = dt_original*factor_upsample;
    ts_original = ts;
    ts = dt:dt:mc.tstamp(end);
    
    for i=1:length(varnames)
        if ~any(strcmp(varnames{i},'tstamp'))
            mc.(varnames{i}) = interp1(ts_original,mc.(varnames{i}),ts,'pchip');
        end
    end
    mc.tstamp = interp1(ts_original,mc.tstamp,ts,'pchip'); 
    mc = structfun(@(x) x', mc, 'un', 0); 

    %% extract trials
    if ~isempty(t)
        for j=1:length(t.end)
            pretrial = max(t.beg(j) - t.move(j),0) + prs.pretrial; % extract everything from "movement onset - pretrial" or "target onset - pretrial" - whichever is first
            posttrial = prs.posttrial; % extract everything until "t_end + posttrial"
            
            timeindx = mc.tstamp > t.beg(j)-pretrial & mc.tstamp < t.end(j)+posttrial;
            trl(j).mc = structfun(@(x) x(timeindx), mc,'un',0);
        end
        
        % set values prior to target onset to nan
        for i = 1:length(varnames)
            if any(strcmp(varnames{i},{'jsy','jsx'})) 
                trl(j).mc.(varnames{i})(1:floor(pretrial/dt)) = nan;
            end
        end
%         trl(j).mc.tstamp = (dt:dt:length(trl(j).mc.(varnames{1}))*dt)' + ...
%             ((t.beg(j)-pretrial)<0)*(t.beg(j)-pretrial) + ((t.beg(j)-pretrial)>0)*(-pretrial); % because not enough pretrial before 1st trial

        emptytrl = [];
        for n = 1:length(trl)
            if ~isempty(trl(n).mc.tstamp)
                % start time at target onset for each trial
                trl(n).mc.tstamp =  trl(n).mc.tstamp - t.beg(n);
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
    mc.tstamp = nan;
    
    mc.JS_X_Raw = nan; % forward JS input scaled by vmax, looks converted
    mc.JS_Yaw_Raw = nan; % rotational JS input scaled by wmax, looks converted
    
    varnames = fieldnames(mc);
    
    for j=1:length(t.end)
        for i=1:length(varnames)
            trl(j).mc.(varnames{i}) = mc.(varnames{i});
            %             trl(j).continuous.(chnames{i}) = downsample(trl(j).continuous.(chnames{i}),prs.factor_downsample);
        end
    end
    
    disp('MC file is missing for this block !!')
    
end









