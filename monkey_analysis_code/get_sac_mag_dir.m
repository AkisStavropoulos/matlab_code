function [trials] = get_sac_mag_dir(trials,presac_t,postsac_t)
%% calculate saccade magnitude and direction
% add them to the trials.events struct

dt = trials(1).continuous.ts(2) - trials(1).continuous.ts(1);

presac = floor(presac_t/dt);
postsac = floor(postsac_t/dt);

for i = 1:numel(trials)
    
    sacindx = arrayfun(@(x) find(trials(i).continuous.ts >= x, 1), trials(i).events.t_sac);
%     sacindx = sacindx(sacindx > presac & sacindx < numel(trials(i).continuous.ts)-postsac);
    
    if ~isempty(sacindx)
        
        trlength = numel(trials(i).continuous.ts);
        
        ye = 0.5*(trials(i).continuous.yle+trials(i).continuous.yre);
        ze = 0.5*(trials(i).continuous.zle+trials(i).continuous.zre);
        
        for n = 1:numel(trials(i).events.t_sac)
            pre_s = sacindx(n) - presac;    if pre_s <= 0; pre_s = 1; end
            post_s = sacindx(n) + postsac;  if post_s > trlength; post_s = trlength; end
            
            y_mag(n) = ye(post_s) - ye(pre_s);
            z_mag(n) = ze(post_s) - ze(pre_s);
        end
        
        trials(i).events.sac_mag = sqrt(y_mag.^2 + z_mag.^2);
        trials(i).events.sac_mag = trials(i).events.sac_mag(:)';
        
        trials(i).events.sac_dir = atan2d(z_mag,y_mag); % angle from the x axis
        trials(i).events.sac_dir = trials(i).events.sac_dir(:)';
        
    else
        
        trials(i).events.sac_mag = [];
        trials(i).events.sac_dir = [];
    
    end
end