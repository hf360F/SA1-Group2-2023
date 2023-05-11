function psixy = psipv(xc, yc, Gamma, x, y)
r = sqrt((x-xc)^2 + (y-yc)^2);
psixy = -Gamma*log(r^2)/(4*pi);
end