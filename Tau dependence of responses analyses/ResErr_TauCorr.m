function [rho_r,rho_th,pval_r,pval_th,Rsq_r,Rsq_th] = ResErr_TauCorr(subject,params,plt,enclose_data)
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
      tau = []; r_res = []; th_res = []; r_sub = []; r_tar = []; th_tar = [];    th_sub = [];
      for j = 1:length(indx)
          tau(j) = subject(i).trials(indx(j)).prs.tau;
          r_sub(j) = subject(i).trials(indx(j)).prs.r_sub;
          r_tar(j) = subject(i).trials(indx(j)).prs.r_tar;
         r_res(j) = r_sub(j) - bias.r(i,s)*r_tar(j); % residual error
         
         th_sub(j) = subject(i).trials(indx(j)).prs.th_sub;
         th_tar(j) = subject(i).trials(indx(j)).prs.th_tar;
         if th_tar(j) >= 0
             th_res(j) = th_sub(j) - bias.th(i,s)*th_tar(j);
         else
             th_res(j) = -(th_sub(j) - bias.th(i,s)*th_tar(j));
         end
      end
      
      % regression
      x = tau(:); y = r_res(:);
      X = [x ones(length(x),1)];
      [b,~,~,~,stats_r]=regress(y,X);  Rsq_r(i,s) = stats_r(1);
      
      x = tau(:); y = th_res(:);
      X = [x ones(length(x),1)];
      [g,~,~,~,stats_th]=regress(y,X); Rsq_th(i,s) = stats_th(1);


      if plt
          subplot(1,2,1);hold on;
          plot(tau,r_res,'.','Color',colr(s,:),'MarkerSize',3);xlabel('tau [s]');ylabel('residual error [cm]');
          xlim([0 5]);ylim([-300 300]);
          if enclose_data
              % draw ellipse
              [Xellipse,Yellipse,x0,y0] = error_ellipse(tau,r_res,conf_int,0);
              plot(Xellipse + x0, Yellipse + y0,'Color',colr(s,:),'HandleVisibility','off');
          end   
              % draw regression
              c = b(2);
              b = b(1);
              plot([-4:8],[-4:8]*b + c,'Color',colr(s,:),'LineWidth',2,'HandleVisibility','off');    hline(0,'k--');
                    
          subplot(1,2,2);hold on;
          plot(tau,th_res,'.','Color',colr(s,:),'MarkerSize',3);xlabel('tau [s]');ylabel('residual error [deg]');
          xlim([0 5]);ylim([-15 15]); 
          if enclose_data
              % draw ellipse
              [Xellipse,Yellipse,x0,y0] = error_ellipse(tau,th_res,conf_int,0);
              plot(Xellipse + x0, Yellipse + y0,'Color',colr(s,:),'HandleVisibility','off');
          end  
              % draw regression
              c = g(2);
              g = g(1);
              plot([-4:8],[-4:8]*g + c,'Color',colr(s,:),'LineWidth',2,'HandleVisibility','off');    hline(0,'k--');
          
      end
      % compute correlation
      [rho_r(i,s),pval_r(i,s)] = nancorr(tau(:),r_res(:));
      [rho_th(i,s),pval_th(i,s)] = nancorr(tau(:),th_res(:));
      title_r_temp{s} = ['\rho = ' num2str(rho_r(i,s),'%4.3f') ', Pval = ' num2str(pval_r(i,s),'%4.3f')];
      title_th_temp{s} = ['\rho = ' num2str(rho_th(i,s),'%4.3f') ', Pval = ' num2str(pval_th(i,s),'%4.3f')];

   end
   if plt;  subplot(1,2,1);title(title_r_temp(:));  subplot(1,2,2);title(title_th_temp(:));   end
   if plt;  legend(legend_input{:},'location','southeast');    suptitle(subject(i).name);  end
end