function [r,th,tau] = GetResponses(subject,params)


[poolindx,legend_input] = get_poolindx(subject,params); 

Nsubs = length(subject);
Nstim = size(poolindx,2);

for i = 1:Nsubs
    for s = 1:Nstim
        
        indx = poolindx{i,s};
        tau{i,s} = arrayfun(@(x) x.prs.tau, subject(i).trials(indx));
        r.sub{i,s} = arrayfun(@(x) x.prs.r_sub, subject(i).trials(indx));
        r.tar{i,s} = arrayfun(@(x) x.prs.r_tar, subject(i).trials(indx));
        th.sub{i,s} = arrayfun(@(x) x.prs.th_sub, subject(i).trials(indx));
        th.tar{i,s} = arrayfun(@(x) x.prs.th_tar, subject(i).trials(indx));
        
    end
end