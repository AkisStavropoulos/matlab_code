function [rad_acc_smo,ang_acc_smo] = acceleration_estimate(trial_id,continuous,prs)
%     function that estimates acceleration from noisy velocities

   vel = continuous(trial_id).v;
   ang_vel = continuous(trial_id).w;
   
   n_pts = length(vel);
   f = 30./n_pts;
   xx = (1:n_pts) * 0.006;
   
   sm = smooth(xx,vel,f,'loess');
   acc_smo = diff(sm) / prs.dt;
   rad_acc_smo = [acc_smo(1)*0.5 ; acc_smo];
   
   sm = smooth(xx,ang_vel,f,'loess');
   acc_smo = diff(sm) / prs.dt;
   
   
   ang_acc_smo = [acc_smo(1)*0.5 ; acc_smo];
%    v_thresh = prs.v_thresh;
%    w_thresh = prs.w_thresh;
%    v_time2thresh = prs.v_time2thresh;
%    indx = find(vel > v_thresh, 1);
%    orign_indx = 1:length(vel);
%    while ~isempty(indx)
%        % do something
%        indx_end = indx + find(abs(vel(indx:end)) < v_thresh & abs(ang_vel(indx:end)) < w_thresh,1); % first downward threshold-crossing
%        if isempty(indx_end)
%            indx_end = length(vel);
%  
%        end
%        n_pts = indx_end - indx + 1;
%        if n_pts == 1
%            break
%        end
%        f = 0.2;%min([45,n_pts]);
%        xx = (1:n_pts) * 0.006;
%        yy = vel(indx:indx_end);
%        sm = smooth(xx,yy,f,'loess');
%        %XX = zeros(n_pts, 2);
%        
%        %XX = lowess(XX,0.1);
%        acc_smo = diff(sm) / prs.dt;
%        acc_smo = [acc_smo(1)*0.5 ; acc_smo];
%        
%        rad_acc_smo(orign_indx(indx:indx_end)) = acc_smo;
%        
%        yy = ang_vel(indx:indx_end);
%        sm = smooth(xx,yy,f,'loess');
%        
%        acc_smo = diff(sm) / prs.dt;
%        acc_smo = [acc_smo(1)*0.5 ; acc_smo];
%        
%        ang_acc_smo(orign_indx(indx:indx_end)) = acc_smo;
% 
%        vel = vel(indx_end+1:end);
%        ang_vel = ang_vel(indx_end+1:end);
%        orign_indx = orign_indx(indx_end+1:end);
%        indx = find(vel > v_thresh, 1);
       
   end
