function [infa, infb] = refpaninf(delta, X, Yin)
% Influence coefficients of a panel of length delta on a point (X, Yin)

% Numerical check on Y, to avoid evaluating atan(infinity)
if abs(Yin) < 1e-8
    Y = 1e-8;
else 
    Y = Yin;
end

I0 = -1/(4*pi)*(X*log(X^2 + Y^2) - (X-delta)*log((X-delta)^2 + Y^2) - 2*delta + 2*Y*(atan(X/Y) - atan((X-delta)/Y)));
I1 = 1/(8*pi)*((X^2 + Y^2)*log(X^2 + Y^2) - ((X - delta)^2 + Y^2)*log((X - delta)^2 + Y^2) - 2*delta*X + delta^2);

infa = (1- X/delta) * I0 - I1/delta;
infb = X/delta*I0 + I1/delta;

end