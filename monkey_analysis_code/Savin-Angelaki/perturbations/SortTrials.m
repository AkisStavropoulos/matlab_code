function [gdtrls, rewtrls, targON_indx] = SortTrials(trials, max_dur, min_dist, max_err_ratio)

% index where trial starts
startindx = arrayfun(@(x) find(~isnan(x.continuous.xmp),1), trials);

% timed-out
% max_dur = 8;
timeout_indx = arrayfun(@(x) x.continuous.ts(end), trials) >= max_dur;

% not moved
% min_dist = 40;
x_sub = arrayfun(@(x,i) x.continuous.xmp(end) - x.continuous.xmp(i), trials, startindx);
y_sub = arrayfun(@(x,i) x.continuous.ymp(end) - x.continuous.ymp(i), trials, startindx);
nomove_indx = sqrt(x_sub.^2 + y_sub.^2) <= min_dist;

% big errors
% max_err_ratio = 0.7;% trials(1).prs.rew_boundary/100;
x_tar = arrayfun(@(x,i) x.prs.xfp - x.continuous.xmp(i), trials, startindx);
y_tar = arrayfun(@(x,i) x.prs.yfp - x.continuous.ymp(i), trials, startindx);
tar_dist = sqrt(x_tar.^2 + y_tar.^2);
bigerr_indx = sqrt((x_tar-x_sub).^2 + (y_tar-y_sub).^2) >= tar_dist*max_err_ratio;

% bad target positions
max_ang = 35;
min_d = 100;
max_d = 400;
tar_ang = atan2d(x_tar,y_tar);
bad_targ_indx = tar_ang > max_ang | tar_ang < -max_ang | tar_dist > max_d | tar_dist < min_d;

% target stayed ON
targON_indx = arrayfun(@(x) x.logical.firefly_fullON, trials);

% keep good trials
gdtrls = ~(timeout_indx | nomove_indx | bigerr_indx | bad_targ_indx | targON_indx); % ~any([timeout_indx ; nomove_indx ; bigerr_indx ; bad_targ_indx ; targON_indx]);

% rewarded trials
rewtrls = arrayfun(@(x) x.logical.reward > 0, trials);