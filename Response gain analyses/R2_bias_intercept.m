function [R2, DR2, F] = R2_bias_intercept(subject,params)

%% R^2 (for polyorder=1) with and without intercept
% params = 'stimtype';
polyorder = 1;
plt = 0;

intercept0 = 0;
[~,~,fval_r_wo,fval_th_wo,~,r_sub,~,th_sub] = ErrScatterFit(subject,params,polyorder,intercept0,plt);
intercept1 = 1;
[~,~,fval_r_w,fval_th_w,~,~,~,~] = ErrScatterFit(subject,params,polyorder,intercept1,plt);

for i = 1:length(subject)
   for s = 1:size(r_sub,2)
      r2_r_wo(i,s) = rsquared(r_sub{i,s}(:),fval_r_wo{i,s}{1}(:));  r2_th_wo(i,s) = rsquared(th_sub{i,s}(:),fval_th_wo{i,s}{1}(:)); 
      
      r2_r_w(i,s) = rsquared(r_sub{i,s}(:),fval_r_w{i,s}{1}(:));    r2_th_w(i,s) = rsquared(th_sub{i,s}(:),fval_th_w{i,s}{1}(:));
   end  
end
F = []; p = []; R2 = []; DR2 = [];
for s = 1:size(r_sub,2)
    R2.r{s} = [r2_r_wo(:,s) r2_r_w(:,s)];
    R2.th{s} = [r2_th_wo(:,s) r2_th_w(:,s)];
    DR2.r{s} = diff(R2.r{s},[],2);
    DR2.th{s} = diff(R2.th{s},[],2);
    for i = 1:length(subject)
        ntrials = length(r_sub{i});
        for j = 2:size(R2.r{s},2)
            F.r{s}(i,j) = (R2.r{s}(i,j-1) - R2.r{s}(i,j))/(R2.r{s}(i,j)/(ntrials-j));       p.r{s}(i,j) = 1-fcdf(R2.r{s}(i,j),1,ntrials);
            F.th{s}(i,j) = (R2.th{s}(i,j-1) - R2.th{s}(i,j))/(R2.th{s}(i,j)/(ntrials-j));   p.th{s}(i,j) = 1-fcdf(R2.th{s}(i,j),1,ntrials);
        end
    end
end
% intercept doesn't increase the variance of data explained by the multiplicative model