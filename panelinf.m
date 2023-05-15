function [infa infb] = panelinf(xa, ya, xb, yb, x, y)
% Influence coefficients of a general vortex panel with edges at (xa, ya)
% and (xb, yb), on a point (x, y)

% Co-ordinates relative to first edge of vortex sheet
r = [x - xa, y - ya];

% Tangent and normal vectors to vortex sheet
t = [xb - xa, yb - ya];
n = [ya - yb, xb - xa];

% Normalise to unit length
t = t/norm(t);
n = n/norm(n);

% Relative co-ordinates in frame aligned with vortex sheet
X = dot(r,t);
Y = dot(r,n);

% Length of vortex sheet from end co-ordinates
delta = power(power((xb - xa),2) + power((yb - ya),2),(1/2));

[infa, infb] = refpaninf(delta, X, Y);

end