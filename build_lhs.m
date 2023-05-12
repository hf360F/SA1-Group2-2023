function lhsmat = build_lhs(xs, ys)
% Construct the matrix A for the equation A*gamma = b, where gamma is the
% vector of unknown vortex panel strengths and b is the vector of 
% streamfunction contributions from each panel needed to enforce no
% through flow at each panel.

np = length(xs) - 1;
psip = zeros(np, np+1);

for  i = 1:np
    for j = 1:(np+1)
        if j == 1
            [infa, infb] = panelinf(xs(j), ys(j), xs(j+1), ys(j+1), xs(i), ys(i));
            psip(i, j) = infa;
         elseif j == np+1
             [infa, infb] = panelinf(xs(j-1), ys(j-1), xs(j), ys(j), xs(i), ys(i));
            psip(i, j) = infb;
        else
            [infa, infb2] = panelinf(xs(j), ys(j), xs(j+1), ys(j+1), xs(i), ys(i));
            [infa2, infb] = panelinf(xs(j-1), ys(j-1), xs(j), ys(j), xs(i), ys(i));
            psip(i, j) = infa + infb;
        end
    end
end

lhsmat = zeros(np+1, np+1);

for i = 1:(np-1)
    for j = 1:(np+1)
        lhsmat(i, j) = psip(i+1, j) - psip(i, j);
    end
end

lhsmat(np, :) = 0; % Stuck on last two rows of A
lhsmat(np, :) = 0;

end