function [rho_r,rho_th,pval_r,pval_th,Rsq_r,Rsq_th] = ResErr_TauCorr(subject,params,plt,enclose_data)
%% Residual error vs Tau
% instead of scatterplot, show how tau affects accuracy this way 
% plt = 0;
% params = [];
conf_int = .68;

[poolindx,legend_input] = get_poolindx(subject,params); 

[r_res,th_res,tau] = ComputeResidualErrors(subject,params);

colr = brewermap(size(poolindx,2),'Dark2');

for i = 1:length(subject)
    if plt; figure; end
   for s = 1:size(poolindx,2)
      
      % regression
      x = tau{i,s}(:); y = r_res{i,s}(:);
      X = [x ones(length(x),1)];
      [b,~,~,~,stats_r]=regress(y,X);  Rsq_r(i,s) = stats_r(1);
      
      x = tau{i,s}(:); y = th_res{i,s}(:);
      X = [x ones(length(x),1)];
      [g,~,~,~,stats_th]=regress(y,X); Rsq_th(i,s) = stats_th(1);


      if plt
          subplot(1,2,1);hold on;
          plot(tau{i,s},r_res{i,s},'.','Color',colr(s,:),'MarkerSize',3);xlabel('tau [s]');ylabel('residual error [cm]');
          xlim([0 5]);ylim([-300 300]);
          if enclose_data
              % draw ellipse
              [Xellipse,Yellipse,x0,y0] = error_ellipse(tau{i,s},r_res{i,s},conf_int,0);
              plot(Xellipse + x0, Yellipse + y0,'Color',colr(s,:),'HandleVisibility','off');
          end   
              % draw regression
              c = b(2);
              d = b(1);
              plot([-4:8],[-4:8]*d + c,'Color',colr(s,:),'LineWidth',2,'HandleVisibility','off');    hline(0,'k--');
                    
          subplot(1,2,2);hold on;
          plot(tau{i,s},th_res{i,s},'.','Color',colr(s,:),'MarkerSize',3);xlabel('tau [s]');ylabel('residual error [deg]');
          xlim([0 5]);ylim([-15 15]); 
          if enclose_data
              % draw ellipse
              [Xellipse,Yellipse,x0,y0] = error_ellipse(tau{i,s},th_res{i,s},conf_int,0);
              plot(Xellipse + x0, Yellipse + y0,'Color',colr(s,:),'HandleVisibility','off');
          end  
              % draw regression
              c = g(2);
              e = g(1);
              plot([-4:8],[-4:8]*e + c,'Color',colr(s,:),'LineWidth',2,'HandleVisibility','off');    hline(0,'k--');
          
      end
      % compute correlation
      [rho_r(i,s),pval_r(i,s)] = nancorr(tau{i,s}(:),r_res{i,s}(:));
      [rho_th(i,s),pval_th(i,s)] = nancorr(tau{i,s}(:),th_res{i,s}(:));
      title_r_temp{s} = ['\rho = ' num2str(rho_r(i,s),'%4.3f') ', Pval = ' num2str(pval_r(i,s),'%4.3f')];
      title_th_temp{s} = ['\rho = ' num2str(rho_th(i,s),'%4.3f') ', Pval = ' num2str(pval_th(i,s),'%4.3f')];

   end
   if plt;  subplot(1,2,1);title(title_r_temp(:));  subplot(1,2,2);title(title_th_temp(:));   end
   if plt;  legend(legend_input{:},'location','southeast');    suptitle(subject(i).name);  end
end
