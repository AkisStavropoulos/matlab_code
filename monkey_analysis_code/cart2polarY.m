function [r,theta] = cart2polarY(x,y)

% Transform cartesian to polar coordinates with respect to Y axis

r = sqrt(x.^2 + y.^2);
theta = atan2d(x,y);