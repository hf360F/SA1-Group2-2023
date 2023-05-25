clear;
close all;

xmin = -2.5;
xmax = 2.5;
ymin = -2.5;
ymax = 2.5;
xc = 0.75;
yc = 0.50;
Gamma = 3.0;
nx = 51;
ny = 41;

for i = 1:nx
    for j = 1:ny
        xm(i,j) = xmin + (i-1)*(xmax-xmin)/(nx-1);
        ym(i,j) = ymin + (j-1)*(ymax-ymin)/(ny-1);
        [a(i,j), b(i,j)] = refpaninf(1.5, xm(i,j), ym(i,j));
    end
end

c = -0.15:0.05:0.15;
contour(xm, ym, a, c);

figure(2);
contour(xm, ym, b, c);

