function [b_r,stats_r,b_th,stats_th,parcor] = MultiRegrFit(subject,params,full_model,int,plt)
% coefficients are spit out for METERS!!!!! or not...
% int: intercept 1 yes 0 no

% full_model = 0;
% params = 'stimtype';    plt = 0;

[poolindx,legend_input] = get_poolindx(subject,params);
colr = brewermap(size(poolindx,2),'Dark2');
b_r = []; stats_r = []; b_th = []; stats_th = [];
for i = 1:length(subject)
    if plt; figure; hold on; end
   for s = 1:size(poolindx,2) 
       indx = poolindx{i,s};
       r_tar = [];  r_sub = []; th_tar = [];  th_sub = []; tau = [];
       for j = 1:length(indx)
           r_tar(j) = subject(i).trials(indx(j)).prs.r_tar; r_sub(j) = subject(i).trials(indx(j)).prs.r_sub;
           th_tar(j) = subject(i).trials(indx(j)).prs.th_tar; th_sub(j) = subject(i).trials(indx(j)).prs.th_sub;
           tau(j) = subject(i).trials(indx(j)).prs.tau; 
       end
       r_tar = .01*r_tar';  r_sub = .01*r_sub'; tau = tau'; % convert to meters
       r_tar = r_tar(:);  r_sub = r_sub(:);  tau = tau(:); % columnize
       y_r = (r_sub(:));  y_th = (th_sub(:));
       if full_model
           X_r = ([r_tar(:) tau(:) r_tar(:).*tau(:)]);  
           X_th = ([th_tar(:) tau(:) th_tar(:).*tau(:)]);  
           if int
           X_r = [ones(length(r_tar),1) X_r];  X_th = [ones(length(th_tar),1) X_th];    
           else
           X_r = [zeros(length(r_tar),1) X_r];  X_th = [zeros(length(th_tar),1) X_th];
           end
       else
           X_r = ([r_tar(:) r_tar(:).*(tau(:))]);  X_r = [zeros(length(r_tar),1) X_r];
           X_th = ([th_tar(:) th_tar(:).*(tau(:))]);  X_th = [zeros(length(th_tar),1) X_th];
       end
       % r_tar(:) tau(:) r_tar(:).*tau(:) zeros(length(r_tar),1)
       [b_r{i,s},~,~,~,stats_r{i,s} ] = regress(y_r,X_r);
       [b_th{i,s},~,~,~,stats_th{i,s} ] = regress(y_th,X_th);
       
       % Partial correlations
       [parcor.r.rho{i,s},parcor.r.p{i,s}] = partialcorr([y_r X_r],'Rows','all');
       [parcor.th.rho{i,s},parcor.th.p{i,s}] = partialcorr([y_th X_th],'Rows','all');
       if full_model
       parcor.r.rho{i,s} = array2table(parcor.r.rho{i,s}, ...
           'VariableNames',{'r_sub','c','r_tar','tau','r_tarXtau'},...
           'RowNames',{'r_sub','c','r_tar','tau','r_tarXtau'});
       parcor.th.rho{i,s} = array2table(parcor.th.rho{i,s}, ...
           'VariableNames',{'th_sub','c','th_tar','tau','th_tarXtau'},...
           'RowNames',{'th_sub','c','th_tar','tau','th_tarXtau'});
       else
           parcor.r.rho{i,s} = array2table(parcor.r.rho{i,s}, ...
           'VariableNames',{'r_sub','c','r_tar','r_tarXtau'},...
           'RowNames',{'r_sub','c','r_tar','r_tarXtau'});
       parcor.th.rho{i,s} = array2table(parcor.th.rho{i,s}, ...
           'VariableNames',{'th_sub','c','th_tar','th_tarXtau'},...
           'RowNames',{'th_sub','c','th_tar','th_tarXtau'});
       end
       if plt
       % plot multilinear regression
       plot3(r_tar,tau,r_sub,'.','Color',colr(s,:));    view(40,40);
       x1fit = 2.5:.2:5.5;
       x2fit = .5:.5:6;
       [X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
       YFIT = b_r{i,s}(1) + b_r{i,s}(2)*X1FIT + b_r{i,s}(3)*X2FIT + b_r{i,s}(4)*X1FIT.*X2FIT;
       
       mesh(X1FIT,X2FIT,YFIT,'EdgeColor',colr(s,:));
       end
   end
   if plt
       xlabel('r_{tar} [m]'); ylabel('\tau [s]');  zlabel('r_{sub} [m]');
       title([ num2str([b_r{i,:}])]);  suptitle(subject(i).name)
   end
end
