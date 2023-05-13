% Vortex sheets with aribtrary position and angle.
% Analytical solution.

clear;
close all;

% Domain
xmin = 0; xmax = 5; ymin = 0; ymax = 4;
nx = 51; ny = 41; % Domain discretisation

% Vortex sheet
xa = 4.1; ya = 1.3; xb = 2.2; yb = 2.9;

for i = 1:nx
    for j = 1:ny
        xm(i,j) = xmin + (i-1)*(xmax-xmin)/(nx-1);
        ym(i,j) = ymin + (j-1)*(ymax-ymin)/(ny-1);
        [a(i,j), b(i,j)] = panelinf(xa, ya, xb, yb, xm(i,j), ym(i,j)); % Streamfunction influence coefficients
    end
end

% Plot contours of influence coefficients
c = -0.15:0.05:0.15;
contour(xm, ym, a, c);
axis("equal")

figure(2);
contour(xm, ym, b, c);
axis("equal")