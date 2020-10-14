function [r_res,th_res,tau] = ComputeResidualErrors(subject,params)


[poolindx,legend_input] = get_poolindx(subject,params); 
polyorder = 1;  intercept = 0;    
[bias.r,bias.th] = ErrScatterFit(subject,params,polyorder,intercept,0);

Nsubs = length(subject);
Nstim = size(poolindx,2);

for i = 1:Nsubs
    for s = 1:Nstim
        
        indx = poolindx{i,s};
        tau{i,s} = arrayfun(@(x) x.prs.tau, subject(i).trials(indx));
        r_sub{i,s} = arrayfun(@(x) x.prs.r_sub, subject(i).trials(indx));
        r_tar{i,s} = arrayfun(@(x) x.prs.r_tar, subject(i).trials(indx));
        th_sub{i,s} = arrayfun(@(x) x.prs.th_sub, subject(i).trials(indx));
        th_tar{i,s} = arrayfun(@(x) x.prs.th_tar, subject(i).trials(indx));
        
        % residual error
        r_res{i,s} = r_sub{i,s} - bias.r(i,s)*r_tar{i,s}; 
        anglesign = sign(th_tar{i,s});
        th_res{i,s} = anglesign.*(th_sub{i,s} - bias.th(i,s)*th_tar{i,s});

    end
end