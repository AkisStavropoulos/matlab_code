function RealVsEstimatedTau(tau,tau_estimate,plt)
%% Plot Real vs Estimated Tau (model implied)

colr = brewermap(size(tau_estimate,2),'Dark2');
for i = 1:size(tau_estimate,1)
    for s = 1:size(tau_estimate,2) 
       y = tau_estimate{i,s}(:);
       X = [ones(length(tau),1) tau(:) tau(:).^2 tau(:).^3];
       c{i,s} = regress(y,X); 
       if plt
       plot(tau,tau_estimate{i,s},'Color',colr(s,:),'linewidth',1); 
       xlabel('Real \tau');    ylabel('Estimated \tau');
       plot(0:6,0:6,'k--');
       xx = 0:6;    xx = [ones(length(xx),1) xx(:) xx(:).^2 xx(:).^3];
%        plot(xx(:,2),xx*c{i,s},'o-','Color',colr(s,:))
       end
    end
   if 0;  suptitle([subject(i).name ' - Model type = ' num2str(model_type(n))]); else; suptitle('All subjects');  end
end

% All subjects
xx = 0.5:.2:7;    xx = [ones(length(xx),1) xx(:) xx(:).^2 xx(:).^3];
figure; hold on;   axis equal;   xlim([0 6]);  ylim([0 6]);   plot(0:8,0:8,'k--'); 
for s = 1:size(tau_estimate,2)
    PercTau = xx*[c{:,s}];    mu = mean(PercTau,2);    SEM = std(PercTau,[],2)/sqrt(size(tau_estimate,1));
    shadedErrorBar(xx(:,2),mu,SEM,'lineprops',{'Color',colr(s,:)});
end
xlabel('Real \tau');    ylabel('Estimated \tau');   title('All subjects'); 

