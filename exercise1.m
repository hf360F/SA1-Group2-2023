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
        psi(i,j) = psipv(xc,yc,Gamma,xm(i,j),ym(i,j));
    end
end

c = -0.4:0.2:1.2;
contour(xm, ym, psi, c);