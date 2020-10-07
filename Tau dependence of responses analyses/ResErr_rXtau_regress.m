function [coef_r,coef_th,rho_r,rho_th,pval_r,pval_th,Rsq_r,Rsq_th] = ResErr_rXtau_regress(subject,params,plt)
%% Residual error vs Tau
% instead of scatterplot, show how tau affects accuracy this way 
% plt = 0;
% params = [];
conf_int = .68;

[poolindx,legend_input] = get_poolindx(subject,params); 
[bias,skipsub] = get_bias(subject,params,0);

subject(skipsub) = [];
bias.r(skipsub,:) = [];             bias.th(skipsub,:) = [];

colr = brewermap(size(poolindx,2),'Dark2');

for i = 1:length(subject)
    if plt; figure; end
   for s = 1:size(poolindx,2)
       indx = poolindx{i,s};
      tau = []; r_res = []; th_res = []; r_sub = []; r_tar = []; th_tar = []; sign_tar = [];  th_sub = [];
      for j = 1:length(indx)
          tau(j) = subject(i).trials(indx(j)).prs.tau;
          r_sub(j) = subject(i).trials(indx(j)).prs.r_sub;
          r_tar(j) = subject(i).trials(indx(j)).prs.r_tar;
         r_res(j) = r_sub(j) - bias.r(i,s)*r_tar(j); % residual error
         
         th_tar(j) = subject(i).trials(indx(j)).prs.th_tar;
         sign_tar(j) = sign(th_tar(j));
         th_sub(j) = subject(i).trials(indx(j)).prs.th_sub;
         th_res(j) = sign_tar(j).*(th_sub(j) - bias.th(i,s)*th_tar(j));

      end
      
      % regression
      x = r_tar(:).*tau(:); y = r_res(:);
      X = [x zeros(length(x),1)];
      [b,~,~,~,stats_r]=regress(y,X);  Rsq_r(i,s) = stats_r(1);
      coef_r(i,s) = b(1);
      
      x = abs(th_tar(:)).*tau(:); y = th_res(:);
      X = [x zeros(length(x),1)];
      [g,~,~,~,stats_th]=regress(y,X); Rsq_th(i,s) = stats_th(1);
      coef_th(i,s) = g(1);

      if plt
          subplot(1,2,1);hold on;
          plot(r_tar.*tau,r_res,'.','Color',colr(s,:),'MarkerSize',3);xlabel('r*tau [cm]');ylabel('residual error [cm]');
          xlim([0 2000]);ylim([-300 300]);
          % draw regression
          plot([0:2000],[0:2000]*b(1) + b(2),'Color',colr(s,:),'LineWidth',2,'HandleVisibility','off');
          
          subplot(1,2,2);hold on;
          plot(abs(th_tar).*tau,th_res,'.','Color',colr(s,:),'MarkerSize',3);xlabel('|\theta|*tau [deg]');ylabel('residual error [deg]');
          xlim([0 200]);ylim([-10 10]);
          % draw regression
          plot([0:200],[0:200]*g(1) + g(2),'Color',colr(s,:),'LineWidth',2,'HandleVisibility','off');
          
      end
      % compute correlation
      [rho_r(i,s),pval_r(i,s)] = nancorr(r_tar(:).*tau(:),r_res(:));
      [rho_th(i,s),pval_th(i,s)] = nancorr(abs(th_tar(:)).*tau(:),th_res(:));
      title_r_temp{s} = ['\rho = ' num2str(rho_r(i,s),'%4.3f') ', Pval = ' num2str(pval_r(i,s),'%4.3f')];
      title_th_temp{s} = ['\rho = ' num2str(rho_th(i,s),'%4.3f') ', Pval = ' num2str(pval_th(i,s),'%4.3f')];

   end
   if plt;  subplot(1,2,1);title(title_r_temp(:));  subplot(1,2,2);title(title_th_temp(:));   end
   if plt;  legend(legend_input{:});    suptitle(subject(i).name);  end
end
