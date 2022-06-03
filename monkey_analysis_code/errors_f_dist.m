function errors_f_dist(trials,condition,r_tar,r_sub,theta_tar,theta_sub)

%% plot radial and angular error as a function of target distance

if isempty(condition) % for all trials
    
%     for i = 1:length(trials)
%         taus(i) = trials(i).prs.tau;
%     end
    
    % radial error
    x = r_tar';
    X = [x ones(length(x),1)];
    y = (r_sub - r_tar)';
    [b]=regress(y,X);
    c = b(2);
    b = b(1);
    figure;
    subplot(2,1,1);plot(r_tar,r_sub-r_tar,'.');hline(0);ylim([-300 300]);
    hold on;plot(x,x*b + c,'r');
    xlabel('target distance (cm)');ylabel('radial error (cm)');title('radial error as a function of distance')
    
    % angular error
    x = r_tar';
    X = [x ones(length(x),1)];
    y = (abs(theta_sub) - abs(theta_tar))';
    [b]=regress(y,X);
    c = b(2);
    b = b(1);
    subplot(2,1,2);plot(r_tar,theta_sub-theta_tar,'.');hline(0);ylim([-30 30]);
    hold on;plot(x,x*b + c,'r');
    xlabel('target distance (cm)');ylabel('rotational error (deg)');title('angular error as a function of distance');
    suptitle(trials(1).prs.subject);
else
    condnames = fieldnames(r_tar);
    % name the conditions for the legend
    if any(strcmp(condnames,'s1bin1')) % check for intersection first
        for i = 1:length(condnames)
            conditionS{i} = condnames{i}(1:2);
            conditionJS{i} = condnames{i}(3:end);
        end
        if any(strcmp(conditionS,'s1')) && any(strcmp(conditionS,'s2'))
            conditionStim = [{'vestibular'} {'visual'} {'combined'} ];
        end
        if any(strcmp(conditionJS,'bin1')) && any(strcmp(conditionJS,'bin2'))
            for i = 1:length(condnames)
                conditionJS{i} = ['\tau ' conditionJS{i}];
            end
        end
        conditionJS = unique(conditionJS);
        count = 0;
        for j = 1:length(conditionStim)
            for n = 1:length(conditionJS)
                count = count + 1;
                cond{count} = [conditionStim{j} ' and ' conditionJS{n}];
            end
        end
    else
        if any(strcmp(condnames,'s1')) && any(strcmp(condnames,'s2'))
            cond = [{'vestibular'} {'visual'} {'combined'} ];
        elseif any(strcmp(condnames,'bin1')) && any(strcmp(condnames,'bin2'))
            for i = 1:length(condnames)
                cond{i} = ['\tau ' condnames{i}];
            end
        end
    end
    %%
    for j = 1:length(condnames)
%         taus = [];
        x = [];
        X = [];
        y = [];
%         for i = 1:length(condition.(condnames{j}))
%             taus(i) = trials(condition.(condnames{j})(i)).prs.tau;
%         end
        figure;
        % radial error
        x = r_tar.(condnames{j})';
        X = [x ones(length(x),1)];
        y = (r_sub.(condnames{j}) - r_tar.(condnames{j}))';
        [b]=regress(y,X);
        c = b(2);
        b = b(1);
        subplot(2,1,1);plot(r_tar.(condnames{j}), r_sub.(condnames{j}) - r_tar.(condnames{j}),'.');hline(0);ylim([-300 300]);
        hold on;plot(x,x*b + c,'r');
        xlabel('target distance (cm)');ylabel('radial error (cm)');title('radial error as a function of distance')
        
        % angular error
        x = r_tar.(condnames{j})';
        X = [x ones(length(x),1)];
        y = (abs(theta_sub.(condnames{j})) - abs(theta_tar.(condnames{j})))';
        [b]=regress(y,X);
        c = b(2);
        b = b(1);
        subplot(2,1,2);plot(r_tar.(condnames{j}), theta_sub.(condnames{j}) - theta_tar.(condnames{j}),'.');hline(0);ylim([-30 30]);
        hold on;plot(x,x*b + c,'r');
        xlabel('target distance (cm)');ylabel('rotational error (deg)');title('angular error as a function of distance');
        suptitle([trials(1).prs.subject ' - ' cond{j}]);
        % obtain figures to rearrange them on the screen later
        fig(j) = gcf;
        fig(j).Units = 'centimeters';
        
        
    end
    % place figures on Right screen in vertical order
    if any(strcmp(condnames,'s1'))
        count = 0;
        for i = 1:length(fig)
            figdim = fig(i).Position(3:4);
            fig(i).Position = [-1.5+(count*(figdim(1)-2)) 15 figdim];  % 50: RL position , 18: max top position
            count = count + 1;
        end
    elseif any(strcmp(condnames,'bin1'))
        count = 0;
        for i = 1:length(fig)
            figdim = fig(i).Position(3:4);
            fig(i).Position = [-.5+(count*(figdim(1)-.5)) 2 figdim]; % 0: max left position , 2: max down position
            count = count + 1;
        end
    else
    figposx =  [-1.5 -1.5 -1.5 ...
                -1.5 -1.5 -1.5...
                50 50 50 ...
                50 50 50];
    figposy =  [15 15 15 ...
                1 1 1 ...
                15 15 15 ...
                3 3 3];
        count1 = 0;
        for n = 1:length(fig)
            figdim = fig(n).Position(3:4);
            fig(n).Position = [figposx(n)+(count1*(figdim(1)-2)) figposy(n) figdim]; % figposx: initial LR position , figposy: initial top-down position
            %         fig(i).Position = [-1.5+(count*(figdim(1)-2)) 15 figdim];  % 50: RL position , 18: max top position
            count1 = count1 + 1;
            if ~rem(n,3); count1 = 0; end
        end
    end
end