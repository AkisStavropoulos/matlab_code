function [out_structure] = addEvents(in_structure)


% phases
% 0 = begin
% 1 = trial
% 2 = stop moving
% 3 = question (for humans)

n_trial = unique(in_structure.trial_num);
idx_trial = find(diff(in_structure.trial_num));
start_trial = [1; idx_trial];
start_trial = start_trial + 1;
end_trial = [start_trial(2:end)-2; size(in_structure.trial_num,1)]; % we are gonna give us a bit of a buffer, 2 samples, so we don't see the jump.
for i = 1:size(start_trial,1)
    trial_phase =  in_structure.phase(start_trial(i):end_trial(i));
    try % in case they never stop...
        stop_trial(i) = start_trial(i)+find(trial_phase == 2, 1)-2; % again for buffer from the jump
    catch
        stop_trial(i) = nan;
    end
end

out_structure = in_structure;
out_structure.start_trial = start_trial(1:end);
out_structure.end_trial = end_trial(1:end);
out_structure.stop_trial = stop_trial(1:end)';
end

