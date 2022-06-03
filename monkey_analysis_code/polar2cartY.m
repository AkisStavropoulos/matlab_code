function [x,y] = polar2cartY(r,theta)
% transform polar to cartesian coordinates with respect to Y axis

x = r.*sind(theta);
y = r.*cosd(theta);

