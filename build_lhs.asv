function lhsmat = build_lhs(xs, ys)
% Construct the matrix A for the equation A*gamma = b, where gamma is the
% vector of unknown vortex panel strengths and b is the vector of 
% streamfunction contributions from each panel needed to enforce no
% through flow at each panel.

np = length(xs) - 1;
psip = zeros(np, np+1);

for  i = 1:(np)
    for j = 1:(np)
        [infa, infb] = panelinf(xs(j), ys(j), xs(j+1), ys(j+1), xs(i), ys(i));
        psip(i, j) = psip(i, j) + infa;
        psip(i, j+1) = infb;
    end
end

lhsmat = zeros(np+1, np+1);

lhsmat(1,2) = 0.5;
%lhsmat(1,3) = -0.5;
%lhsmat(1,np-1) = -0.5;
lhsmat(1,np-1) = -0.5;

lhsmat(np+1,1) = 0.5;
%lhsmat(np+1,3) = -0.5;
%lhsmat(np+1,np-1) = -0.5;
lhsmat(np+1,np) = 0.5;

lhsmat(2:np,:) = psip(2:np,:) - psip(1:(np-1),:);
lhsmat(np+1,np+1) = 1;

end