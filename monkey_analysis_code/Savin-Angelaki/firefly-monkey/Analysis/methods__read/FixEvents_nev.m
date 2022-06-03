function events_nev = FixEvents_nev(events_nev,behv_trials)

events_smr = cell2mat({behv_trials.events});
ntrls_smr = length(events_smr);
ntrls_nev = length(events_nev.t_beg);

if (ntrls_smr - ntrls_nev) < 0
    warning('Mismatch in number of trials in NEV and SMR -- could not fix the problem!');

else
    ntrls_missed = ntrls_smr - ntrls_nev;
    % check if neural recording was not started on time
    iti_smr = diff([events_smr.t_beg]); iti_smr = iti_smr((ntrls_missed+1):end);
    iti_nev = diff(events_nev.t_beg);
    badindx = (iti_smr>10 | iti_nev>10);
    iti_smr = iti_smr(~badindx); iti_nev = iti_nev(~badindx);
    if corr(iti_smr(:),iti_nev(:)) > 0.99 % this level of corr. is unlikely by chance
        warning('Mismatch in number of trials in NEV and SMR -- problem most likely fixed!');
        % yes, neural recording was not started on time
        events_nev.t_beg = [nan(1,ntrls_missed) events_nev.t_beg];
        events_nev.t_end = [nan(1,ntrls_missed) events_nev.t_end];
    
    else
        iti_smr = diff([events_smr.t_beg]); iti_smr = iti_smr(1:(end-ntrls_missed));
        iti_nev = diff(events_nev.t_beg);
        badindx = (iti_smr>10 | iti_nev>10);
        iti_smr = iti_smr(~badindx); iti_nev = iti_nev(~badindx);
        trial_start = 1;
        while corr(iti_smr(trial_start:min(trial_start+100, length(iti_smr)))',iti_nev(trial_start: min(trial_start+100, length(iti_smr)))') > 0.95 && trial_start <= (length(iti_nev)-100)
            trial_start = trial_start + 1;
        end
        trial_end = trial_start - 1;
        if corr(iti_smr(1:trial_end)',iti_nev(1:trial_end)') > 0.99 % this level of corr. is unlikely by chance
            warning('Mismatch in number of trials in NEV and SMR -- problem most likely fixed!');
            % yes, neural recording was not started on time
            events_nev.t_beg = [events_nev.t_beg nan(1,ntrls_missed)];
            events_nev.t_end = [events_nev.t_end nan(1,ntrls_missed)];
            events_nev.t_beg(trial_end+1:end) = nan;
            events_nev.t_end(trial_end+1:end) = nan;
        end
        % recheck the correlation
        iti_smr = diff([events_smr.t_beg]);
        iti_nev = diff(events_nev.t_beg);
        badindx = (iti_smr>10 | iti_nev>10 | isnan(iti_nev));
        iti_smr = iti_smr(~badindx); iti_nev = iti_nev(~badindx);
        if corr(iti_smr(:),iti_nev(:)) < 0.99
            ME = MException('TrailMismatch:InconsistentDataType', ...
            'Mismatch in trial number could not be resolved');
            throw(ME)
        end
    end
    
    ntrls_end_nev = length(events_nev.t_end);
    
    if (ntrls_smr == (ntrls_end_nev + 1)) 
        % possibly recording interrupted before in the nev
        iti_nev = diff(events_nev.t_beg);
        iti_smr = diff([events_smr.t_beg]);
        badindx = (iti_smr>10 | iti_nev>10 | isnan(iti_nev));
        iti_smr = iti_smr(~badindx); iti_nev = iti_nev(~badindx);
        if corr(iti_smr(:),iti_nev(:)) > 0.99
            warning('Mismatch in number of trials in NEV and SMR -- problem most likely fixed!');
            % yes, neural recording was not started on time
            events_nev.t_end = [events_nev.t_end nan(1,1)];
            events_nev.t_beg(end) = nan;
        end
    end
 end
end

