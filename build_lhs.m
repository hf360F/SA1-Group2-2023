function lhsmat = build_lhs(xs, ys)
% Construct the matrix A for the equation A*gamma = b, where gamma is the
% vector of unknown vortex panel strengths and b is the vector of 
% streamfunction contributions from each panel needed to enforce no
% through flow at each panel.

np = length(xs) - 1;
psip = zeros(np, np+1);

for  i = 1:(np)
    %for j = 1:(np)
     %   [infa, infb] = panelinf(xs(j), ys(j), xs(j+1), ys(j+1), xs(i), ys(i));
      %  psip(i, j) = psip(i, j) + infa;
       % psip(i, j+1) = infb;
    %end
    [infa, infb] = panelinf(xs(1:np), ys(1:np), xs(2:np+1), ys(2:np+1), xs(i), ys(i));
    psip(:,1:np) = psip(:,1:np) + infa;
    psip(:,2:np+1) = infb;
end

lhsmat = zeros(np+1, np+1);

lhsmat(1,1) = 1;
lhsmat(2:np,:) = psip(2:np,:) - psip(1:(np-1),:);
lhsmat(np+1,np+1) = 1;

end