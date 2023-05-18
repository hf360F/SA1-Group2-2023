% Zero pressure graident, turbulent boundary layer

clear;
close all;

global ReL ue0 duedx;

ReL = 1E7;
ue0 = 1;
duedx = 0;

x0 = 0.01;
thick0(1) = 0.037*x0*(ReL*x0)^(-1/5);
thick0(2) = 1.80*thick0(1);

[delx, thickhist] = ode45(@thickdash,[0 0.99],thick0);

theta = thickhist(:,1);
theta7 = 0.037*times(delx,power(ReL*delx, -1/5));
theta9 = 0.023*times(delx,power(ReL*delx, -1/6));

He = thickhist(:,2)/theta;

for i =1:length(delx)
    if He(i) > 1.46
        delx(i)
        break
    end
end

%plot(delx, theta)
%hold on
plot(delx, He)
yline(1.46)
xlabel("x/L")
ylabel("theta")
legend("\theta", "\He")
