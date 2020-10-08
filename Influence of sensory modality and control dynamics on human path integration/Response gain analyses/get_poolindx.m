function [poolindx,legend_input] = get_poolindx(subject,cond,varargin)
%% Get a pool of trial indices of the condition you want for each  subject
% poolindx{i,k}: Rows: subject, Cols: condition bins

if any(strcmp(cond,{'stimtype','tau','stimtau'}))
    %% MODALITIES
    if strcmp(cond,'stimtype')
        
        condition = fieldnames(subject(1).stimtype);
        condition(end) = [];
        condnames = {'vestibular' 'visual', 'combined'};
        
    elseif strcmp(cond,'tau')
        
        condition = fieldnames(subject(1).tau.bins);
        limnames = fieldnames(subject(1).tau.lims);
        
        for i = 1:length(limnames)
            condnames{i} = ['\tau ' condition{i} ' - [' num2str(subject(1).tau.lims.(limnames{i})(1),'%4.2f') ' - ' num2str(subject(1).tau.lims.(limnames{i})(2),'%4.2f') ']'];
        end
        
    elseif strcmp(cond,'stimtau')
        
        condition = fieldnames(subject(1).stimtau.bins);
        limnames = fieldnames(subject(1).stimtau.lims);
        stimnames = {'vestibular' 'visual', 'combined'};
        
        count = 0;
        for s = 1:length(stimnames)
            for i = 1:length(limnames)
                count = count + 1;
                condnames{count} = [stimnames{s} ' and ' '\tau ' condition{count} ' - [' num2str(subject(1).tau.lims.(limnames{i})(1),'%4.2f') ' - ' num2str(subject(1).tau.lims.(limnames{i})(2),'%4.2f') ']'];
            end
        end
    end
    %% Calculate poolindx and legend_input
    for i = 1:length(subject)
        for s = 1:length(condition)
            if strcmp(cond,'stimtype')
                poolindx{i,s} = subject(i).(cond).(condition{s});
                
            else
                poolindx{i,s} = subject(i).(cond).bins.(condition{s});
            end
        end
    end
    
    legend_input = condnames;
elseif strcmp(cond,'angle')
    %% ANGLE
    % 2 extra args: 1) 'abs' value or 'notabs', 2) 'subject' or 'target' angle
    % default selection: 'abs' and 'target'
    
    % choose subject angle or target angle
    if isempty(varargin)
        th_field = 'th_tar';
    elseif length(varargin) == 1
        th_field = 'th_tar';
    else
        if strcmp(varargin{2},'target')
            th_field = 'th_tar';
        elseif strcmp(varargin{2},'subject')
        th_field = 'th_sub';
        end
    end
    
    if strcmp(cond,'angle')
        if isempty(varargin) || strcmp(varargin{1},'abs')
            M = 2; % number of bins
            angbins = linspace(0,40,M+1);
            
            for i = 1:length(subject)
                for s = 1:M
                    count = 0;
                    for j = 1:length(subject(i).trials)
                        if (abs(subject(i).trials(j).prs.(th_field)) >= angbins(s)) && (abs(subject(i).trials(j).prs.(th_field)) < angbins(s+1))
                            count = count + 1;
                            poolindx{i,s}(count) = j;
                            
                        end
                    end
                    condnames{s} = ['\theta_t_a_r bin' num2str(s) ' - [' num2str(angbins(s)) ' - ' num2str(angbins(s+1)) 'deg]'];
                end
            end
            
            legend_input = condnames;
            
        elseif strcmp(varargin{1},'notabs')
            M = 4; % number of bins
            angbins = linspace(-40,40,M+1);
            
            for i = 1:length(subject)
                for s = 1:M
                    count = 0;
                    for j = 1:length(subject(i).trials)
                        if ((subject(i).trials(j).prs.(th_field)) >= angbins(s)) && ((subject(i).trials(j).prs.(th_field)) < angbins(s+1))
                            count = count + 1;
                            poolindx{i,s}(count) = j;
                            
                        end
                    end
                    condnames{s} = ['\theta_t_a_r bin' num2str(s) ' - [' num2str(angbins(s)) ' - ' num2str(angbins(s+1)) 'deg]'];
                end
            end
            legend_input = condnames;
        end
     
    end
elseif strcmp(cond,'distance')
    %% DISTANCE
    if strcmp(cond,'distance')
        M = 2; % number of bins
        distbins = linspace(250,550,M+1);
        
        for i = 1:length(subject)
            for s = 1:M
                count = 0;
                for j = 1:length(subject(i).trials)
                    if (subject(i).trials(j).prs.r_tar >= distbins(s)) && (subject(i).trials(j).prs.r_tar < distbins(s+1))
                        count = count + 1;
                        poolindx{i,s}(count) = j;
                        
                    end
                end
                condnames{s} = ['d_t_a_r bin' num2str(s) ' - [' num2str(distbins(s)) ' - ' num2str(distbins(s+1)) 'm]'];
            end
        end
        
        legend_input = condnames;
    end
elseif strcmp(cond,'SPC')
    %% SPC
        if strcmp(cond,'SPC')
            if ~isempty(varargin)
                % input the dividing edges of the bins
                thresh = varargin{1};
            else
                thresh = [.75];
            end
        M = length(thresh) + 1; % number of bins
        SPCbins = [0 thresh 10];
        
        for i = 1:length(subject)
            for s = 1:M
                count = 0;
                for j = 1:length(subject(i).trials)
                    if (subject(i).trials(j).stats.SPC >= SPCbins(s)) && (subject(i).trials(j).stats.SPC < SPCbins(s+1))
                        count = count + 1;
                        poolindx{i,s}(count) = j;
                        
                    end
                end
                if s == M
                    condnames{s} = ['SPC bin' num2str(s) ' > ' num2str(thresh(s-1))];
                else
                    condnames{s} = ['SPC bin' num2str(s) ' < ' num2str(thresh(s))];
                end
            end
        end
        
        legend_input = condnames;
        end
elseif isempty(cond)
    
    for i = 1:length(subject)
       for j = 1:length(subject(i).trials) 
          poolindx{i,1}(j) = j; 
           
       end
    end
    legend_input = {'all data'};
else
    error('Wrong condition string.')
end