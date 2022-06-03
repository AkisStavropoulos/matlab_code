function y = logit(x)
%% single variable logistic regression

if any(x<0) || any(x>1); error('x range must be bounded within 0<x<1.'); end

y = log(x./(1-x));