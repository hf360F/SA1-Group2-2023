% Panel method for a cylinder

clear;
close all;

% Domain
xmin = -2.5; xmax = 2.5; ymin = -2; ymax = 2;
nx = 51; ny = 41; % Discretisation

% Cylinder
np = 100; % Panel discretisation
theta = (0:np)*2*pi/np; % Angles
xs = cos(theta); ys = sin(theta); % Points on unit circle 
gamma = -2*sin(theta); % Variation of sheet strength around circumference

for i = 1:nx
    for j = 1:ny
        xm(i,j) = xmin + (i-1)*(xmax-xmin)/(nx-1);
        ym(i,j) = ymin + (j-1)*(ymax-ymin)/(ny-1);
        % Contribution of free-stream to streamfunction
        psi(i, j) = ym(i, j);
        for n = 0:np
            % Inf. coeffs. due to panel between i and i+1'th coordinate pairs
            [infa, infb] = panelinf(xs(i), ys(i), xs(i+1), ys(i+1), xm(i,j), ym(i,j));
            % Contribution of this panel to streamfunction
            psi(i,j) = psi(i, j) + infa*gamma(i) + infb*gamma(i+1);
        end
    end
end

% Plot contours of influence coefficients
c = -1.75:0.25:1.75
contour(xm, ym, psi, c);
hold on
plot(xs,ys)
hold off