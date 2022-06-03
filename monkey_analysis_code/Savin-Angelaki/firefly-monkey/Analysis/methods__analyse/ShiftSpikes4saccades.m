function trials = ShiftSpikes4saccades(trials,eventtimes)
% shifts spike trains by eventtimes - used to align spike trains to events

%% check event times are in a cell array, whose each cell is events form one trial
if length(trials)~=length(eventtimes)
    fprintf('error: event times should be a cell array of same length as trials \n');
    return;
end

%% shift spike train on each trial by the event time on that trial
for i=1:length(trials)
    tspk_temp = trials(i).tspk;
    
    if isempty(tspk_temp)
        trials(i).tspk = nan(1,numel(eventtimes{i}));    
    else
        if numel(eventtimes{i}) > 0
            for n = 1:numel(eventtimes{i})
                trials(i).tspk(:,n) = tspk_temp - eventtimes{i}(n);
            end
        else
            trials(i).tspk = [];
        end
    end
end