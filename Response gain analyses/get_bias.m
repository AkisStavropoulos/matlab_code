function [b,skipsub] = get_bias(subject,params)

%% Get bias slope for either 1st button push or end of trial
    polyorder = 1;  intercept = 0;    plt = 0;
    [b.r,b.th] = ErrScatterFit(subject,params,polyorder,intercept,plt);
    skipsub = [];

