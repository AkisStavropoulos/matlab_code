function vmax = findvmax(x,T,tau)

% calculate vmax based on x,T,tau
% x: in cm

if tau > 0
% vmax = (x/T)*( 1 / (-1 + 2*(tau/T) * log((1 + exp(T/tau))/2))) ;
vmax = (x/T)*( 1 / (2*(tau/T) * log(cosh(T/(2*tau))))) ;

elseif tau == 0
    vmax = x/T;
end