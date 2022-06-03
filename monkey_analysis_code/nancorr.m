function [r,p] = nancorr(x,y)
% [rho,pval] = nancorr(x,y)
% correlation between two vectors with nans
indx = ~isnan(x) & ~isnan(y);

if sum(indx)==0
    r = nan;    p = nan;
else
    [r,p] = corr(x(indx),y(indx));
end
