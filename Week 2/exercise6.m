% Test of laminar boundary layer with separation or transition

clear;
%close all;

global ReL ue0 duedx;

ReL = 1E5;
ue0 = 1;
duedx = -0.5;

L = 1;

nx = 101;

x = linspace(0, 1/L, nx);
ue = 1 + x.*duedx;
theta = zeros(nx, 1);
thetaBlas = zeros(nx, 1);
ReTheta = zeros(nx, 1);

intTot = 0;
i = 1;
laminar = true;

int = 0;
ils = 0;
itr = 0;
its = 0;

He = zeros(nx,1);
He(1) = 1.57258;

while laminar && i < (nx)
    i = i + 1;

    ua = ue(i-1);
    ub = ue(i);
    xa = x(i-1);
    xb = x(i);

    % Momentum thickness by Blasius, Thwaites
    intTot = intTot + ueintbit(xa, ua, xb, ub); % Euler for Thwaites integral
    theta(i) = sqrt((0.45*(ub^(-6))*intTot/ReL));
    %thetaBlas(i) = 0.664*(xa^0.5)/(ReL^(0.5));

    ReTheta(i) = ReL*ub*theta(i);

    % Find energy shape factor from Thwaites gradient factor
    m = -ReL*(theta(i))^2*(ub-ua)/(xb-xa);
    H = thwaites_lookup(m);
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

while i < (nx)
    i = i + 1;
    He(i) = He(its);
    uea = ue0 + duedx*(i-2)/nx;
    ueb = uea + duedx/nx;

    theta(i) = theta(i-1)*(uea/ueb)^(He(i)+2);
end

if int ~= 0
    disp(['Natural transition at ' num2str(x(int)) ' with Rethet ' num2str(ReTheta(int))])
    % figure(1)
    % plot(x(int), theta(int), "x")
    % legend("Natural transition")
    % hold on;
    % figure(2)
    % plot(x(int), He(int), "x")
    % legend("Natural transition")
    % hold on;
end

if ils ~= 0
    disp(['Laminar separation at ' num2str(x(ils)) ' with Rethet ' num2str(ReTheta(ils))])
    % figure(1)
    % plot(x(ils), theta(ils), "o")
    % legend("Laminar separation")
    % hold on;
    % figure(2)
    % plot(x(ils), He(ils), "o")
    % legend("Laminar separation")
    % hold on;
end

if itr ~= 0
    disp(['Turbulent reattachment at ' num2str(x(itr))])
    % figure(1)
    % plot(x(itr), theta(itr), "+")
    % legend("Turbulent reattachment")
    % hold on;
    % figure(2)
    % plot(x(itr), He(itr), "+")
    % hold on;
end

if its ~= 0
    disp(['Turbulent separation at ' num2str(x(its))])
    % figure(1)
    % plot(x(its), theta(its), "*")
    % hold on;
    % figure(2)
    % plot(x(its), He(its), "*")
    % hold on;
end

% figure(1)
% plot(x(1:nx), theta(1:nx))
% xlabel("x/L")
% ylabel("\theta/L")
% legend("Laminar separation", "Turbulent separation", "Re = 1E4", "Laminar separation", "Turbulent reattachment", "Re = 1E5", "Natural transition", "Re = 1E6")
% hold on;
% 
% figure(2)
% plot(x(1:nx), He(1:nx))
% xlabel("x/L")
% ylabel("H_E")
% legend("Laminar separation", "Turbulent separation", "Re = 1E4", "Laminar separation", "Turbulent reattachment", "Re = 1E5", "Natural transition", "Re = 1E6")
% hold on;