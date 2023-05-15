% Panel method for a cylinder at incidence

clear;
close all;

% Domain
xmin = -2.5; xmax = 2.5; ymin = -2; ymax = 2;
nx = 101; ny = 81; % Discretisation

% Cylinder
np = 100; % Panel discretisation
theta = (0:np)*2*pi/np; % Angles
xs = cos(theta); ys = sin(theta); % Points on unit circle 
alpha = pi/20;

A = build_lhs(xs, ys);
b = build_rhs(xs, ys, alpha);

gam = A\b;

delta = ((xs(2) - xs(1))^2 + (ys(2) - ys(1))^2)^(1/2);

for i = 1:nx
    for j = 1:ny
        xm(i,j) = xmin + (i-1)*(xmax-xmin)/(nx-1);
        ym(i,j) = ymin + (j-1)*(ymax-ymin)/(ny-1);
        % Contribution of free-stream to streamfunction
        psi(i, j) = ym(i, j)*cos(alpha) - xm(i, j)*sin(alpha);
        for n = 1:np
            % Inf. coeffs. due to panel between i and i+1'th coordinate pairs
            [infa, infb] = panelinf(xs(n), ys(n), xs(n+1), ys(n+1), xm(i,j), ym(i,j));
            % Contribution of this panel to streamfunction
            psi(i,j) = psi(i, j) + infa*gam(n) + infb*gam(n+1);
        end
    end
end

% Plot contours of streamlines
c = -1.75:0.25:1.75;
contour(xm, ym, psi, c);
hold on
axis("equal")
xlabel("x")
ylabel("y")
plot(xs, ys)
hold off

figure(2);
plot(theta/pi, gam)
xlabel("Angle around cylinder /pi, radians")
ylabel("Vortex sheet strength")

totalCirc = 0;
for n = 1:(np-1)
    gamAvg = gam(n)+gam(n+1)/2;
    totalCirc = totalCirc + gamAvg * delta;
end

totalCirc