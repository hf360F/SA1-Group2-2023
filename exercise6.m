% Test of laminar boundary layer with separation or transition

clear;
close all;

global ReL ue0 duedx;

ReL = 1E8;
ue0 = 1;
duedx = -0.5;

ReL = 3.2E5;
L = 1;

nx = 101;

ue = linspace(1, 0.5, nx);
x = linspace(0, 1/L, nx);
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
    m = -ReL*(theta(i)/L)^2*(ub-ua)/(xb-xa);
    H = thwaites_lookup(m);
    He = laminar_He(H);
    
    if log(ReTheta(i)) >= 18.4*He - 21.74
        laminar = false;
        disp([xa ReTheta(i)/1000])
        int = i;

    elseif m >= 0.9
        laminar = false;
        ils = i;
        He(i) = 1.51509;
    end
end

delta = times(He, theta);


while its == 0 && i < nx
    i = i+1;
    thick0(1) = 0.037*x(i)*(ReL*x(i))^(-1/5);
    thick0(2) = 1.80*thick0(1);
    [delx, thickhist] = ode45(@thickdash,[0,x(i)-x(i-1)],thick0);
    theta = thickhist( ...
        i,1);
    delta = thickhist(i,2);
    He(i) = delta./theta;
    if He(i) < 1.46 && its == 0
        its = i;
    end
    if He(i) > 1.58 && itr == 0
        itr = i;
    end
end

if int ~= 0
    disp(['Natural transition at ' num2str(x(int)) ' with Rethet ' num2str(ReTheta(int))])
end

if ils ~= 0
    disp(['Separation at ' num2str(x(ils)) ' with Rethet ' num2str(ReTheta(ils))])
end

if itr ~= 0
    disp(['Turbulent reattachment ' num2str(x(itr)) ' with Rethet ' num2str(ReTheta(ils))])
end

if its ~= 0
    disp(['Turbulent separation at ' num2str(x(its)) ' with Rethet ' num2str(ReTheta(int))])
end

plot(x, theta)
hold on
plot(x, thetaBlas)
xlabel("x/L")
ylabel("theta/L")
legend("Thwaites", "Blasius")
