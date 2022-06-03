function [yle,zle,yre,zre] = world2eye(x,y,z,delta)
% Transform world coordinates to eye coordinates

yle = atan2d(x + delta, sqrt(y.^2 + z^2));
yre = atan2d(x - delta, sqrt(y.^2 + z^2));
zle = atan2d(z , sqrt(y.^2 + (x + delta).^2));
zre = atan2d(z , sqrt(y.^2 + (x - delta).^2));

ver_mean = nanmean([zle , zre],2); % mean vertical eye position (of the two eyes)
hor_mean = nanmean([yle , yre],2); % mean horizontal eye position
