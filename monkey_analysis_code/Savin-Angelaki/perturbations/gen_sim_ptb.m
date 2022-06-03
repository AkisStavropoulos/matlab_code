function [sim_resp,sim_targ] = gen_sim_ptb(x_targ_np,y_targ_np,v_ts_np,w_ts_np,x_targ_p,y_targ_p,t_ptbn,v_ptb,w_ptb)

%% Non Parametric ptb generation
%% params
sf = 1/.006;
%% Template ptb wave

sig =.2/.006; %prs.filtwidth; %filter width
sz = 168;%prs.filtsize; %filter size
t2 =linspace(-sz/2, sz/2, sz);
h = exp(-t2.^2/(2*sig^2));
h = h - h(1);

%% load coordinates




for i = 1:numel(x_targ_p)
    
    dist = sqrt((x_targ_np - x_targ_p(i)).^2 + (y_targ_np - y_targ_p(i)).^2);
    
    n_close = 5; %number of closest trajectories picked
    [~,cand_indx] = mink(dist,n_close);
    indx = cand_indx(randi(n_close));
    
    v = v_ts_np{indx};
    w = w_ts_np{indx};
    
    
   %% Find ptb index
    t = (1/sf):(1/sf):numel(v)*(1/sf); t =t(:);
    ptb_indx = find(t >= t_ptbn(i),1);
    
    %case that the ptb happens after the end of the trial
    if isempty(ptb_indx)
        t = (1/sf):(1/sf):t_ptbn(i);
        ptb_indx = numel(t);
        v(end+1:ptb_indx) = 0; % zero pad until the beginning of the ptb
        w(end+1:ptb_indx) = 0;
        
    end
    
    %% generate ptb waves
    v_wave = h*v_ptb(i);
    w_wave = h*w_ptb(i);
    
    %% Add ptb_waves to the trajectories
    
    if (ptb_indx + numel(h) - 1) > numel(t)     % zero pad if there are not enough dat points
        
        pad_indx = ptb_indx + numel(h) - 1 - numel(t);
        
        v((end + 1) : (end+pad_indx)) = 0;
        w((end + 1) : (end+pad_indx)) = 0;
        
        v(ptb_indx:(ptb_indx + numel(h) - 1)) = v(ptb_indx:(ptb_indx + numel(h) - 1)) + v_wave';
        w(ptb_indx:(ptb_indx + numel(h) - 1)) = w(ptb_indx:(ptb_indx + numel(h) - 1)) + w_wave';
        
    else
        
        v(ptb_indx:(ptb_indx + numel(h) - 1)) = v(ptb_indx:(ptb_indx + numel(h) - 1)) + v_wave';
        w(ptb_indx:(ptb_indx + numel(h) - 1)) = w(ptb_indx:(ptb_indx + numel(h) - 1)) + w_wave';
        
    end
    %% generate trajectories here
    [~, ~, x_sim.resp(i), y_sim.resp(i), ~, ~] = gen_traj(w,v, t);
    
    x_sim.targ(i) = x_targ_np(indx);
    y_sim.targ(i) = y_targ_np(indx);
    


end

sim_resp.r = sqrt(x_sim.resp.^2+y_sim.resp.^2) ;       sim_resp.theta = atand(x_sim.resp./y_sim.resp) ;
sim_targ.r = sqrt(x_sim.targ.^2+y_sim.targ.^2) ;       sim_targ.theta = atand(x_sim.targ./y_sim.targ) ;

end