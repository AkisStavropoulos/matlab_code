function [tau,JS,sub,start_pos,exp_sub,realindx,tar] = extract_data_for_model(subject,params)

%% extract data to be used in tau estimation models
% tar is already in different format

[poolindx,legend_input] = get_poolindx(subject,params);

[bias,skipsub] = get_bias(subject,params);

for i = 1:length(subject)
    for s = 1:size(poolindx,2)
        indx = poolindx{i,s};
        
        r_tar = []; th_tar = []; 
        % Remove trials without joystick data
        rmindx = [];
        for j = 1:length(indx)
            JSData = subject(i).trials(indx(j)).mc.JS_Yaw_Raw;
            if isempty(JSData) || (sum(isnan(JSData)) == length(JSData))
                rmindx = [rmindx j];
            end
        end
        indx(rmindx) = [];
        ntrls = length(indx);
        
        realindx{i,s} = indx;
        
        % extract data
        tau{i,s} = arrayfun(@(x) x.prs.tau, subject(i).trials(indx));
        JS{i,s}.x = arrayfun(@(x) x.mc.JS_Yaw_Raw/x.prs.wmax,subject(i).trials(indx),'uniformoutput',false);
        JS{i,s}.y = arrayfun(@(x) x.mc.JS_X_Raw/x.prs.vmax,subject(i).trials(indx),'uniformoutput',false);
        r_tar = arrayfun(@(x) x.prs.r_tar, subject(i).trials(indx));
        th_tar = arrayfun(@(x) x.prs.th_tar, subject(i).trials(indx));
        sub{i,s}.r = arrayfun(@(x) x.prs.r_sub, subject(i).trials(indx));
        sub{i,s}.th = arrayfun(@(x) x.prs.th_sub, subject(i).trials(indx));
        sub{i,s}.x = arrayfun(@(x) x.continuous.xmp(end), subject(i).trials(indx));
        sub{i,s}.y = arrayfun(@(x) x.continuous.ymp(end), subject(i).trials(indx));
        start_pos{i,s}.x = arrayfun(@(x) x.continuous.xmp(1), subject(i).trials(indx));
        start_pos{i,s}.y = arrayfun(@(x) x.continuous.ymp(1), subject(i).trials(indx));
        
        exp_sub{i,s}.r = bias.r(i,s)*r_tar;
        exp_sub{i,s}.th = bias.th(i,s)*th_tar;
        
        % notice the differet format
        tar.r{i,s} = r_tar;
        tar.th{i,s} = th_tar;
        tar.x{i,s} = arrayfun(@(x) x.prs.fireflyposx, subject(i).trials(indx));
        tar.y{i,s} = arrayfun(@(x) x.prs.fireflyposy, subject(i).trials(indx));
        
        % remove outlier trials (awful responses)
        
        
        
    end
end