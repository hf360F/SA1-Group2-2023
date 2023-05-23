function f = ueintbit(xa, ua, xb, ub)
% Approximate contribution to Thwaites boundary layer momentum thickness
% integral between xa and xb, with velocity ratio ua = ue/U at x = xa
% and ub = ue/U at x = xb.

%average velocity
ubar = (ua + ub)/2;
%change in velocity and x
deltau = ub - ua;
deltax = xb - xa;

f = (ubar^5 + 5*(ubar^3)*(deltau^2)/6 + ubar*(deltau^4)/16)*deltax;

end