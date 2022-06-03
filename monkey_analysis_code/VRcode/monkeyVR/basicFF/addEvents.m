function [data_out] = addEvents(data_in)


% phases
% 0 = begin
% 1 = trial
% 2 = stop moving
% 3 = question (for humans)
% 4 = ITI?

n_trial = unique(data_in.trial_num);
idx_trial = find(diff(data_in.trial_num(~isnan(data_in.trial_num))));
start_trial = [idx_trial];
start_trial = start_trial + 1;
end_trial = [start_trial(2:end)-2; size(data_in.trial_num,1)]; % we are gonna give us a bit of a buffer, 2 samples, so we don't see the jump.
for i = 1:size(start_trial,1)
    trial_phase =  data_in.phase(start_trial(i):end_trial(i));
    try % in case they never stop...
        stop_trial(i) = start_trial(i)+find(trial_phase == 2, 1)-2; % again for buffer from the jump
    catch
        stop_trial(i) = nan;
    end
end

data_out = data_in;
data_out.start_trial = start_trial(1:end);
data_out.end_trial = end_trial(1:end);
data_out.stop_trial = stop_trial(1:end)';



