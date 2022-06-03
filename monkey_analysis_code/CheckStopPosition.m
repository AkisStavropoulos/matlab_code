for i = 176
    plot(trials(i).mc.timestamp/trials(i).mc.timestamp(end),trials(i).mc.JS_X_Raw);hold on;
end
% find when subject stopped steering
for i = 1:length(trials)
    flipped_JS = flip(trials(i).mc.JS_X_Raw);
    indxjs{i} = find(flipped_JS ~= 0,1);
    indxjs{i} = length(trials(i).mc.JS_X_Raw) - indxjs{i};
    
    ts(i) = trials(i).mc.timestamp(indxjs{i});
    indxsmr(i) = find(trials(i).continuous.ts >= ts(i),1);
%     indxv{i} = find(trials(i).continuous.v(round(end/2):end) < 20);

end
% find when subject started braking
for i = 1:length(trials)
    indxjs{i} = find(trials(i).mc.JS_X_Raw(round(end/4):end) <= 0,1);
    indxjs{i} = round(length(trials(i).mc.JS_X_Raw)/4) + indxjs{i};
    
    ts(i) = trials(i).mc.timestamp(indxjs{i});
    indxsmr(i) = find(trials(i).continuous.ts >= ts(i),1);
%     indxv{i} = find(trials(i).continuous.v(round(end/2):end) < 20);

end
% find when subjects stopped braking for first time
for i = 1:length(trials)
    indx_temp = find(trials(i).mc.JS_X_Raw(round(end/4):end) <= 0,1);
    indx_temp = round(length(trials(i).mc.JS_X_Raw)/4) + indx_temp;
    indxjs{i} = find(trials(i).mc.JS_X_Raw(indx_temp:end) >= 0,1);
    indxjs{i} = indx_temp + indxjs{i};
    
    ts(i) = trials(i).mc.timestamp(indxjs{i});
    indxsmr(i) = find(trials(i).continuous.ts >= ts(i),1);
%     indxv{i} = find(trials(i).continuous.v(round(end/2):end) < 20);

end
%% calculate new distance and angle traveled
for i = 1:length(trials)
    x_s = trials(i).continuous.xmp(indxsmr(i));
    x_0 = trials(i).continuous.xmp(1);
    y_s = trials(i).continuous.ymp(indxsmr(i));
    y_0 = trials(i).continuous.ymp(1);
r_sub2(i) = sqrt((x_s - x_0).^2 + ...
    (y_s - y_0).^2);

th_sub2(i) = atan2d((x_s - x_0),(y_s - y_0));
end
trials = save_errors2(trials,r_sub2,th_sub2);

scattermod(trials,stimtype,tau,1);
scattermod_js(trials,stimtype,tau,1);

errors_vs_tau(trials,stimtype);
errors_vs_tau_js(trials,stimtype);
