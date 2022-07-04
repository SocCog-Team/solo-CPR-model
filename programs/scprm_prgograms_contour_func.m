% Function that takes mean and variance matrix
% and compute the contours
% Output [X,Y] pair is centers+radii
% Input S: covariance matrix
% mu: mean
% c: std 1, 2, and so on. Or it can be chi square value
% For 95% confidence interval error ellipse, chisquare_val = 2.4477;

function [X,Y] = contour_func(mu,S,c)
%S = [1 0;0 1];
%c = 2.4;
%mu = [0;0];

[V,D] = eig(S*c);
th = linspace(0,2*pi,1000); % 1000 points around ellipse
a = sqrt(D(2,2)); % major axis length
b = sqrt(D(1,1)); % minor axis length
ellipse_x = a*cos(th);
ellipse_y = b*sin(th);
phi = atan2(V(2,2),V(2,1));
R = [cos(phi) sin(phi);-sin(phi) cos(phi)];
ellipse = [ellipse_x; ellipse_y]'*R; % rotated ellipse

X = ellipse(:,1)+mu(1);
Y = ellipse(:,2)+mu(2);