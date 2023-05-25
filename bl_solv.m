function [int ils itr its delstar theta] = bl_solv(x, cp)
global ReL ue0;

nx = length(x);
np = nx;

ReL = 1E5;
ue0 = 1;

x = [0, x];
ue = zeros(nx, 1);
duedx = zeros(nx, 1);
delstar = zeros(nx, 1);
theta = zeros(nx, 1);
thetaBlas = zeros(nx, 1);
ReTheta = zeros(nx, 1);
He = zeros(nx,1);
He(1) = 1.57258;

intTot = 0;
i = 1;
laminar = true;

int = 0;
ils = 0;
itr = 0;
its = 0;

% Initialise, from stagnation point to x(1)
ua = 0;
ub = sqrt(1 - cp(1));
xa = 0;
xb = x(1);

intTot = ueintbit(xa, ua, xb, ub); % Euler for Thwaites integral
theta(1) = (0.45*(ua^(-6))*intTot/ReL)^(0.5);

m1 = -ReL*(theta(1))^2*(ub-ua)/(xb-xa);
H1 = thwaites_lookup(m1);

delstar(1) = theta(1)*H1;

while laminar && i < (nx-1)
    i = i + 1;

    ua = sqrt(1 - cp(i));
    ub = sqrt(1 - cp(i+1));
    xa = x(i);
    xb = x(i+1);

    duedx(i) = (ub - ua)/(xb - xa);

    % Momentum thickness by Blasius, Thwaites
    intTot = intTot + ueintbit(xa, ua, xb, ub); % Euler for Thwaites integral
    theta(i) = (0.45*(ua^(-6))*intTot/ReL)^(0.5);
    thetaBlas(i) = 0.664*(xa^0.5)/(ReL^(0.5));

    ReTheta(i) = ReL*ua*theta(i);

    % Find energy shape factor from Thwaites gradient factor
    m = -ReL*(theta(i))^2*(ub-ua)/(xb-xa);
    H = thwaites_lookup(m);

    delstar(i) = H*theta(i);

    He(i) = laminar_He(H);
    
    if log(ReTheta(i)) >= 18.4*He(i) - 21.74
        laminar = false;
        int = i;

    elseif m >= 0.09
        laminar = false;
        ils = i;
        He(i) = 1.51509;
    end
end

delta = He.*theta;

while its == 0 && i < (nx)
    i = i+1;
    thick0(1) = theta(i-1);
    thick0(2) = delta(i-1);
    ue0 = ue(i-1);
    [delx, thickhist] = ode45(@thickdash,[0,x(i)-x(i-1)],thick0);
    theta(i) = thickhist(end,1);
    delta(i) = thickhist(end,2);
    He(i) = delta(i)/theta(i);

    if He(i) < 1.46
        its = i;
    end
    if He(i) > 1.58 && itr == 0 && int == 0
        itr = i;
    end
end

while i < (nx-1)
    i = i +1;
    He(i) = He(its);
    uea = ue0 + duedx*i/nx;
    ueb = uea + duedx/nx;

    theta(i) = theta(i-1)*(uea/ueb)^(He(i)+2);
end