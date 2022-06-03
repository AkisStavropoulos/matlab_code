%% function to add sessions
function AddSessions(this,monk_id,sess_id,content,prs) % e.g. content = {'behv','lfps','units','pop'}
islfps = any(strcmp(content,'lfps')); isunits = any(strcmp(content,'units')); ispop = any(strcmp(content,'pop'));
allsessions = this.sessions; old_instance = find([allsessions.monk_id] == monk_id & [allsessions.sess_id] == sess_id);
if ~isempty(old_instance)
    ovwrt = logical(input('This session was already analysed once. Press 1 to overwrite, 0 to quit \n'));
    if ovwrt, new_instance = old_instance; % overwrite old instance
    else, return;
    end
else
    nsessions = numel(this.sessions);
    new_instance = nsessions + 1; % create new instance
end
if ~exist('prs','var')
    prs = default_prs(monk_id,sess_id);
end
this.sessions(new_instance) = session(monk_id,sess_id,prs.sess_date);
this.sessions(new_instance).AddBehaviours(prs);
if islfps % load and analyse LFPs
    this.sessions(new_instance).AddLfps(prs);
    this.sessions(new_instance).AnalyseLfps(prs);
end
if isunits % load and analyse neurons
    this.sessions(new_instance).AddUnits(prs);
    % This is for export a mat file include all data of a session
    % Need to set prs.extractonly=true in default_prs file, then it ignores all
    if prs.extractonly
        PrepareData(this);
    else
        this.sessions(new_instance).AnalyseUnits(prs);
    end
end
%         return
% >>>>>>> bc519e6e7a17a0fc89e080316e60bd60f48e58a1
%     else
%         nsessions = numel(this.sessions);
%         new_instance = nsessions + 1; % create new instance
%     end
%     prs = default_prs(monk_id,sess_id);
%     this.sessions(new_instance) = session(monk_id,sess_id,prs.sess_date);
%     this.sessions(new_instance).AddBehaviours(prs);
%     if islfps % load and analyse LFPs
%         this.sessions(new_instance).AddLfps(prs);
%         this.sessions(new_instance).AnalyseLfps(prs);
%     end
%     if isunits % load and analyse neurons
%         this.sessions(new_instance).AddUnits(prs);
%         this.sessions(new_instance).AnalyseUnits(prs);
%     end
    if ispop && isunits, this.sessions(new_instance).AddPopulation('units',prs);
    elseif ispop && islfps, this.sessions(new_instance).AddPopulation('lfps',prs);
    end
end