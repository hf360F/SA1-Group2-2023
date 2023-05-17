% Vortex sheets with aribtrary position and angle.
% Numerical solution.

clear;
close all;

% Domain
xmin = 0; xmax = 5; ymin = 0; ymax = 4;
nx = 51; ny = 41; % Discretisation

% Vortex sheet
nv = 100; % Discretisation
gamma_a = 1; gamma_b = 1; % Start and end vortex sheet strength
xa = 4.1; ya = 1.3; xb = 2.2; yb = 2.9;
delta = ((xb - xa)^2 + (yb - ya)^2)^(1/2);

for i = 1:nx
    for j = 1:ny
        xm(i,j) = xmin + (i-1)*(xmax-xmin)/(nx-1);
        ym(i,j) = ymin + (j-1)*(ymax-ymin)/(ny-1);
        % Need to define psia here to add to it in for loop
        psia(i, j) = 0;
        psib(i, j) = 0;
        for n = 0:nv
            % Discrete vortex co-ordinates
            xc = xa + (xb-xa)*n/100 + (xb-xa)/(nv*2);
            yc = ya + (yb-ya)*n/100 + (yb-ya)/(nv*2);
            % Separate influence coefficients in terms of start and end sheet strength
            psia(i,j) = psia(i, j) + psipv(xc,yc,gamma_a*delta*(1-(n/nv))/nv,xm(i,j),ym(i,j));
            psib(i,j) = psib(i, j) + psipv(xc,yc,gamma_b*delta*n/(nv^2),xm(i,j),ym(i,j));
        end
    end
end

% Plot contours of influence coefficients
c = -0.15:0.05:0.15;
contour(xm, ym, psia, c);
axis("equal")
xlabel("x")
ylabel("y")

figure(2);
contour(xm, ym, psib, c);
axis("equal")
xlabel("x")
ylabel("y")