function sw = switchtime(tau, T)
% compute the time where you switch the joystick control from push to brake, to stop at the target location

for i = 1:length(tau)
    if tau(i)~=0
        sw(i) = tau(i).*log((1 + exp(T(i)./tau(i)))./2);
    elseif tau(i) == 0
        sw(i) = T(i);
    end
end

sw = sw(:);