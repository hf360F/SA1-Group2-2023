function rhsmat = build_rhs(xs, ys, alpha)
% Construct the vector b for the equation A*gamma = b, where gamma is the
% vector of unknown vortex panel strengths and b is the vector of 
% streamfunction contributions from each panel needed to enforce no
% through flow at each panel.

np = length(xs) - 1;

rhsmat = zeros(np+1, 1);

for i = 1:(np-1)
    rhsmat(i+1) = (ys(i) - ys(i+1))*cos(alpha) - (xs(i) - xs(i+1))*sin(alpha);
end

end