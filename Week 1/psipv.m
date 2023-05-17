function psixy = psipv(xc, yc, Gamma, x, y)
% Streamfunction contribution of a line vortex of strength Gamma locate
% at (xc, yc) on (x, y).

r = sqrt((x-xc)^2 + (y-yc)^2);
psixy = -Gamma*log(r^2)/(4*pi);

end