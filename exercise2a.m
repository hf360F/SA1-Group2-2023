% Analytical solution of vortex sheet
% with linearly varying strength.

clear;
close all;

% Domain
xmin = -2.5; xmax = 2.5;
ymin = -2.5; ymax = 2.5;
delta = 1.5; % Sheet length
nx = 51; ny = 41; % Domain discretisation

for i = 1:nx
    for j = 1:ny
        xm(i,j) = xmin + (i-1)*(xmax-xmin)/(nx-1);
        ym(i,j) = ymin + (j-1)*(ymax-ymin)/(ny-1);
        [a(i,j), b(i,j)] = refpaninf(delta, xm(i,j), ym(i,j)); % Streamfunction influence coefficients
    end
end

% Plot contours of influence coefficients
c = -0.15:0.05:0.15;
contour(xm, ym, a, c);
axis("equal")

figure(2);
contour(xm, ym, b, c);
axis("equal")
xlabel("x")
ylabel("y")

