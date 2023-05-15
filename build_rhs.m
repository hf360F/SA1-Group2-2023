function rhsmat = build_rhs(xs, ys, alpha)
% Construct the vector b for the equation A*gamma = b, where gamma is the
% vector of unknown vortex panel strengths and b is the vector of 
% streamfunction contributions from each panel needed to enforce no
% through flow at each panel.

np = length(xs) - 1;

rhsmat = zeros(np+1, 1);

rhsmat(2:np,:) = (ys(1:np-1) - ys(2:np))*cos(alpha) - (xs(1:np-1) - xs(2:np))*sin(alpha);

end