% Test of ueintbit.m with constant ue

ReL = 2500;
L = 1;

nx = 101;

ue = linspace(1, 1, nx);
x = linspace(0, 1, nx);
theta = zeros(nx, 1);
thetaBlas = zeros(nx, 1);

intTot = 0;

for i = 1:(length(x)-1)
    ua = ue(i);
    ub = ue(i+1);
    xa = x(i);
    xb = x(i+1);

    intTot = intTot + ueintbit(xa, ua, xb, ub);
    theta(i) = L*(0.45*(ue(i)^(-6))*intTot/ReL)^(0.5);
    thetaBlas(i) = L*0.664*(xa^0.5)/(ReL^(0.5));
end


plot(x/L, theta/L)
hold on
plot(x/L, thetaBlas/L)
xlabel("x/L")
ylabel("theta/L")
legend("Thwaites", "Blasius")
