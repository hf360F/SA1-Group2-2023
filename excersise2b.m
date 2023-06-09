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
ga = 0.015;
gb = 0.015;


for i = 1:nx
    for j = 1:ny
        xm(i,j) = xmin + (i-1)*(xmax-xmin)/(nx-1);
        ym(i,j) = ymin + (j-1)*(ymax-ymin)/(ny-1);
        psia(i, j) = 0;
        psib(i, j) = 0;

        for n = 0:100
            Gamma = ga*n/100;
            xc = 1.5*n/100;
            psia(i,j) = psia(i, j) + psipv(xc,yc,ga*(1-n)/100,xm(i,j),ym(i,j));
            psib(i,j) = psib(i, j) + psipv(xc,yc,gb*n/100,xm(i,j),ym(i,j));
        end
    end
end

c = -0.15:0.05:0.15;
contour(xm, ym, psia, c);

figure(2);
contour(xm, ym, psib, c);

