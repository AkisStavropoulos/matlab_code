function events_plx = FixEvents_plx(events_plx,behv_trials)

events_smr = cell2mat({behv_trials.events});
ntrls_smr = length(events_smr);
ntrls_plx = length(events_plx.t_beg);

if (ntrls_smr - ntrls_plx) < 0
    warning('Mismatch in number of trials in NEV and SMR -- could not fix the problem!');

else
    ntrls_missed = ntrls_smr - ntrls_plx;
    % check if neural recording was not started on time
    iti_smr = diff([events_smr.t_beg]); iti_smr = iti_smr((ntrls_missed+1):end);
    iti_plx = diff(events_plx.t_beg);
    badindx = (iti_smr>10 | iti_plx>10);
    iti_smr = iti_smr(~badindx); iti_plx = iti_plx(~badindx);
    if corr(iti_smr(:),iti_plx(:)) > 0.99 % this level of corr. is unlikely by chance
        warning('Mismatch in number of trials in NEV and SMR -- problem most likely fixed!');
        % yes, neural recording was not started on time
        events_plx.t_beg = [nan(1,ntrls_missed) events_plx.t_beg];
        events_plx.t_end = [nan(1,ntrls_missed) events_plx.t_end];
        events_plx.t_rew = [nan(1,ntrls_missed) events_plx.t_end];
    end
    ntrls_end_plx = length(events_plx.t_end);
    ntrls_rew_plx = length(events_plx.t_rew);
    if (ntrls_smr == (ntrls_end_plx + 1)) && (ntrls_smr == (ntrls_rew_plx + 1))
        % possibly recording interrupted before in the nev
        iti_plx = diff(events_plx.t_beg);
        iti_smr = diff([events_smr.t_beg]);
        badindx = (iti_smr>10 | iti_plx>10 | isnan(iti_plx));
        iti_smr = iti_smr(~badindx); iti_plx = iti_plx(~badindx);
        if corr(iti_smr(:),iti_plx(:)) > 0.99
            warning('Mismatch in number of trials in NEV and SMR -- problem most likely fixed!');
            % yes, neural recording was not started on time
            events_plx.t_end = [events_plx.t_end nan(1,1)];
            events_plx.t_beg(end) = nan;
            events_plx.t_rew(end) = [events_plx.t_rew nan(1,1)];
            
        end
    elseif (ntrls_smr == (ntrls_end_plx + 1)) && (ntrls_smr == (ntrls_rew_plx))
        % possibly recording interrupted before in the nev
        iti_plx = diff(events_plx.t_beg);
        iti_smr = diff([events_smr.t_beg]);
        badindx = (iti_smr>10 | iti_plx>10 | isnan(iti_plx));
        iti_smr = iti_smr(~badindx); iti_plx = iti_plx(~badindx);
        if corr(iti_smr(:),iti_plx(:)) > 0.99
            warning('Mismatch in number of trials in NEV and SMR -- problem most likely fixed!');
            % yes, neural recording was not started on time
            events_plx.t_end = [events_plx.t_end nan(1,1)];
            events_plx.t_beg(end) = nan;
            events_plx.t_rew(end) = nan;
            
        end
        
    end
 end
end

