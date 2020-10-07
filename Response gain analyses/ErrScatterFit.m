function [b_r,b_th,fval_r,fval_th,r_tar,r_sub,th_tar,th_sub] = ErrScatterFit(subject,params,polyorder,intercept,plt)
%% Fit behavioural bias for radial and angular errors

% params = 'stimtype';
% polyorder = 4; % choose order of polynomials you want to fit
% intercept = 0; % 0 for c=0, 1 for free intercept
% plt = 1;
%
% b_r and b_th: rows: number of subjects
%               cols: number of conditions (poolindx)
%               Each element contains the polynomial coefficients of
%               the corresponding order, e.g. {i,s}{1}: 1st order poly coefs
%                                             {i,s}{2}: 2nd order poly coefs
%                                             {i,s}{n}: n-th order poly coefs

if nargin < 5
    plt = 0;
end

[poolindx,legend_input] = get_poolindx(subject,params); 
colr = brewermap(size(poolindx,2),'Dark2');

r_tar = cell(size(poolindx));   r_sub = cell(size(poolindx));
th_tar = cell(size(poolindx));  th_sub = cell(size(poolindx));
b_r = cell(size(poolindx));     b_th = cell(size(poolindx));
fval_r = cell(size(poolindx)); fval_th = cell(size(poolindx));
for i = 1:length(subject)
    if plt; figure; end
    count = 1;
    for s = 1:size(poolindx,2)
       indx = poolindx{i,s};
       
       for j = 1:length(indx)
           r_tar{i,s}(j) = subject(i).trials(indx(j)).prs.r_tar;
           r_sub{i,s}(j) = subject(i).trials(indx(j)).prs.r_sub;
           th_tar{i,s}(j) = subject(i).trials(indx(j)).prs.th_tar;
           th_sub{i,s}(j) = subject(i).trials(indx(j)).prs.th_sub;
       end
       % Linear regression (or polyfit)
       % radial error
       y_r = r_sub{i,s}(:);
       x_r = r_tar{i,s}(:);
       if intercept
           for n = 1:polyorder
               b_r{i,s}{n} = polyfit(x_r,y_r,n);
               fval_r{i,s}{n} = polyval(b_r{i,s}{n},x_r);
           end
       else
           for n = 1:polyorder
               b_r{i,s}{n} = polyfitZero(x_r,y_r,n);
               fval_r{i,s}{n} = polyval(b_r{i,s}{n},x_r);
           end
           
       end

       % angular error
       y_th = th_sub{i,s}(:);
       x_th = th_tar{i,s}(:);
       if intercept
           for n = 1:polyorder
               b_th{i,s}{n} = polyfit(x_th,y_th,n);
               fval_th{i,s}{n} = polyval(b_th{i,s}{n},x_th);
           end
       else
           for n = 1:polyorder
               b_th{i,s}{n} = polyfitZero(x_th,y_th,n);
               fval_th{i,s}{n} = polyval(b_th{i,s}{n},x_th);
           end
       end
       
       % Plot
       if plt
       % radial error
       subplot(size(poolindx,2),2,count);hold on;
       plot(r_tar{i,s},r_sub{i,s},'.','Color',colr(s,:),'MarkerSize',2);grid on;
       for n = 1:polyorder
           plot(0:800,polyval(b_r{i,s}{n},0:800),'r','LineWidth',2);
       end
       plot(0:1000,0:1000,'k--');
       xlabel('r_t_a_r [cm]');ylabel('r_s_u_b_j [cm]');title(['Bias = ' num2str(b_r{i,s}{1}(1),2)]);
       xlim([0 800]);ylim([0 800]);
       legend(legend_input{s});
       count = count+1;
       
       % angular error
       subplot(size(poolindx,2),2,count);hold on;
       plot(th_tar{i,s},th_sub{i,s},'.','Color',colr(s,:),'MarkerSize',2);grid on;
       for n = 1:polyorder
           plot(-180:180,polyval(b_th{i,s}{n},-180:180),'r','LineWidth',2);
       end
       plot(-180:180,-180:180,'k--');
       xlabel('\theta_t_a_r [cm]');ylabel('\theta_s_u_b_j [cm]');title(['Bias = ' num2str(b_th{i,s}{1}(1),2)]);
       xlim([-50 50]);ylim([-50 50]);
       legend(legend_input{s});
       count = count+1;
       end
    end
    if plt
    suptitle(subject(i).name);
    end
end

%% If only multiplicative bias, convert to matrix
if iscell(b_r)
for i = 1:size(b_r,1)
   for s = 1:size(b_r,2)
       b_r{i,s} = cell2mat(b_r{i,s});   b_r{i,s}(b_r{i,s} == 0) = [];
       b_th{i,s} = cell2mat(b_th{i,s}); b_th{i,s}(b_th{i,s} == 0) = [];
   end
end
b_r = cell2mat(b_r);    b_th = cell2mat(b_th);
end
