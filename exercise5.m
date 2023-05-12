% Panel method for a cylinder at incidence

clear;
close all;

% Domain
xmin = -2.5; xmax = 2.5; ymin = -2; ymax = 2;
nx = 51; ny = 41; % Discretisation

% Cylinder
np = 10; % Panel discretisation
theta = (0:np)*2*pi/np; % Angles
xs = cos(theta); ys = sin(theta); % Points on unit circle 

A = build_lhs(xs, ys);
b = build_rhs(xs, ys, 0);

det(A)

gam = A\b;