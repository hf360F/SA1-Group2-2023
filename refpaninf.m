function [infa, infb] = refpaninf(delta, X, Yin)
% Influence coefficients of a panel of length delta on a point (X, Yin)

% Numerical check on Y, to avoid evaluating atan(infinity)
if abs(Yin) < 1e-8
    Y = 1e-8;
else 
    Y = Yin;
end

I0 = -1/(4*pi)*(times(X, log(power(X,2)) + power(Y,2)) - times((X-delta),log(power((X-delta),2)) + power(Y,2)) - 2*delta + times(2*Y,(atan(X/Y)) - atan((X-delta)/Y)));
I1 = 1/(8*pi)*((power(X,2) + times(power(Y,2),log(power(X,2))) + power(Y,2)) - (power((X - delta),2) + times(power(Y,2),log(power((X - delta),2))) + power(Y,2)) - 2*delta*X + power(delta,2));

infa = times((1 - X\delta), I0 - I1\delta);
infb = times(X\delta,I0) + I1\delta;

end