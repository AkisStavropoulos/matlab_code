%% compute F-statistic and corresponding p-value for a sequence of regression models

%% gather data from all subjects
for i=1:length(subject)
    for j = 1:length(subject(i).trials)
        r_f{i}(j) = subject(i).trials(j).prs.r_tar;
        r_s{i}(j) = subject(i).trials(j).prs.r_sub;
        theta_f{i}(j) = subject(i).trials(j).prs.th_tar;
        theta_s{i}(j) = subject(i).trials(j).prs.th_sub;
    end
    % force columnization
    r_f{i} = r_f{i}(:);
    r_s{i} = r_s{i}(:);
    theta_f{i} = abs(theta_f{i}(:));
    theta_s{i} = abs(theta_s{i}(:));
end

%% radial distance
for i=1:length(subject)
    for j=0:5 % order of polynomial
        if j==0
            P = [1 mean(r_s{i} - r_f{i})];
        else    
            P = polyfit(r_f{i},r_s{i},j);
        end
        P = fliplr([zeros(1,6-length(P)) P]);
        SSE_radial(i,j+1) = mean((r_s{i} - [r_f{i}.^0 r_f{i} r_f{i}.^2 r_f{i}.^3 r_f{i}.^4 r_f{i}.^5]*P').^2);
    end
end

%% absolute angle
for i=1:length(subject)
    for j=0:5
        if j==0
            P = [1 mean(theta_s{i} - theta_f{i})];
        else    
            P = polyfit(theta_f{i},theta_s{i},j);
        end
        P = fliplr([zeros(1,6-length(P)) P]);
        SSE_angular(i,j+1) = mean((theta_s{i} - [theta_f{i}.^0 theta_f{i} theta_f{i}.^2 theta_f{i}.^3 theta_f{i}.^4 theta_f{i}.^5]*P').^2);
    end
end

%% F-statistic
for i=1:length(subject)
    ntrials = length(r_f{i});
    for j=2:6
        F_radial(i,j) = (SSE_radial(i,j-1) - SSE_radial(i,j))/(SSE_radial(i,j)/(ntrials-j)); 
        p_radial(i,j) = 1-fcdf(F_radial(i,j),1,ntrials);
        
        F_angular(i,j) = (SSE_angular(i,j-1) - SSE_angular(i,j))/(SSE_angular(i,j)/(ntrials-j)); 
        p_angular(i,j) = 1-fcdf(F_angular(i,j),1,ntrials);
    end
end