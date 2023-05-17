% Approximation of vortex sheet with linearly varying
% sheet strength as a sum of discrete line vortices.

clear;
close all;

xmin = -2.5; xmax = 2.5; ymin = -2.5; ymax = 2.5; % Domain
delta = 1.5; % Vortex sheet length
nx = 51; ny = 41; nv = 100; % Discretisation of domain and vortex sheet
gamma_a = 1; gamma_b = 1; % Start and end vortex sheet strength
xa = 0; ya = 0; % Co-ordinates of left hand edge of vortex sheet

for i = 1:nx
    for j = 1:ny
        xm(i,j) = xmin + (i-1)*(xmax-xmin)/(nx-1);
        ym(i,j) = ymin + (j-1)*(ymax-ymin)/(ny-1);
        % Need to define psia here to add to it in for loop
        psia(i, j) = 0;
        psib(i, j) = 0;
        for n = 0:nv
            % Discrete vortex co-ordinates
            xc = xa + delta*n/nv + delta/(nv*2);
            yc = ya;
            % Separate influence coefficients in terms of start and end sheet strength
            psia(i, j) = psia(i, j) + psipv(xc,yc,gamma_a*delta*(1-(n/nv))/nv,xm(i,j),ym(i,j));
            psib(i, j) = psib(i, j) + psipv(xc,yc,gamma_b*delta*n/(nv^2),xm(i,j),ym(i,j));
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