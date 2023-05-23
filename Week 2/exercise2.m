% Test of ueintbit.m with varying ue and test for transition

clear;
close all;

ReL = 1E6;
L = 1;

nx = 101;

ue = linspace(1, 1-0.25, nx);
x = linspace(0, 1/L, nx);
theta = zeros(nx, 1);
thetaBlas = zeros(nx, 1);
ReTheta = zeros(nx, 1);

intTot = 0;
i = 1;
laminar = true;

while laminar && i < (nx-1)
    i = i + 1;

    ua = ue(i);
    ub = ue(i+1);
    xa = x(i);
    xb = x(i+1);

    % Momentum thickness by Blasius, Thwaites
    intTot = intTot + ueintbit(xa, ua, xb, ub); % Euler for Thwaites integral
    theta(i) = (0.45*(ua^(-6))*intTot/ReL)^(0.5);
    thetaBlas(i) = 0.664*(xa^0.5)/(ReL^(0.5));

    ReTheta(i) = ReL*ua*theta(i);

    % Find energy shape factor from Thwaites gradient factor
    m = -ReL*(theta(i))^2*(ub-ua)/(xb-xa);
    H = thwaites_lookup(m);
    He = laminar_He(H);
    
    % Test for transition
    if log(ReTheta(i)) >= 18.4*He - 21.74
        laminar = false;
        disp([xa ReTheta(i)/1000])
    end
end

plot(x, theta)
hold on
plot(x, thetaBlas)
xlabel("x/L")
ylabel("theta/L")
legend("Thwaites", "Blasius")
