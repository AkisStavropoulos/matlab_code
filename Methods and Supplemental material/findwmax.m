function wmax = findwmax(theta,T,tau)
% calculate vmax based on x,T,tau
if tau > 0
wmax = (theta/T)*( 1 / (-1 + 2*(tau/T) * log((1 + exp(T/tau))/2))) ;
elseif tau == 0
    wmax = theta/T;
end