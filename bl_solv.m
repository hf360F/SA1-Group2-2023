function [int ils itr its delstar theta] = bl_solv(x, cp)
global Re ue0 duedx;

ue0 = 1;
duedx = 0;

nx = length(x);
np = nx;

ue = zeros(nx, 1);
%duedx = zeros(nx, 1);
delstar = zeros(nx, 1);
theta = zeros(nx, 1);
delta = zeros(nx, 1);
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
theta(1) = (0.45*(ub^(-6))*intTot/Re)^(0.5);

m = -Re*(theta(1))^2*(ub-ua)/(xb-xa);
H = thwaites_lookup(m);

delstar(1) = theta(1)*H;

while laminar && i < (nx-1)
    i = i + 1;

    ua = sqrt(1 - cp(i-1));
    ub = sqrt(1 - cp(i));
    xa = x(i-1);
    xb = x(i);

    duedx = (ub - ua)/(xb - xa);

    % Momentum thickness by Blasius, Thwaites
    intTot = intTot + ueintbit(xa, ua, xb, ub); % Euler for Thwaites integral
    theta(i) = (0.45*(ub^(-6))*intTot/Re)^(0.5);

    % Find energy shape factor from Thwaites gradient factor
    m = -Re*(theta(i))^2*(ub-ua)/(xb-xa);
    H = thwaites_lookup(m);
    ReTheta(i) = Re*ub*theta(i);

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

delta(i) = He(i)*theta(i);

while its == 0 && i < (nx)
    i = i+1;
    
    ua = sqrt(1 - cp(i-1));
    ub = sqrt(1 - cp(i));
    xa = x(i-1);
    xb = x(i);

    thick0(1) = theta(i-1);
    thick0(2) = delta(i-1);
    
    ue0 = ua;
    duedx = (ub - ua)/(xb - xa);

    [delx, thickhist] = ode45(@thickdash,[0,xb-xa],thick0);
    theta(i) = thickhist(end,1);
    delta(i) = thickhist(end,2);
    He(i) = delta(i)/theta(i);

    if He(i) >= 1.46
        H = (11*He(i) + 15)/(48*He(i) - 59);
    else
        H = 2.803;
    end
    delstar(i) = H*theta(i);

    if He(i) < 1.46
        its = i;
    end
    if He(i) > 1.58 && itr == 0 && int == 0
        itr = i;
    end
end

while i < (nx)
    i = i +1;
    He(i) = He(its);
    
    uea = sqrt(1 - cp(i-1));
    ueb = sqrt(1 - cp(i));

    H = 2.803;
    theta(i) = theta(i-1)*(uea/ueb)^(H+2);
    delstar(i) = theta(i)*H;
end