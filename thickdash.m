function dthickdx = thickdash(xmx0, thick)
   global ReL ue0 duedx;
   
   %dimensionless boundary layer edge velocity at given x
   ue = ue0 + duedx*xmx0;
   ReTheta = ReL * ue * thick(1);
   %He = deltae/theta is energy shape factor
   He = thick(2)/thick(1);
   
   %Eppler-Somers
   if He >= 1.46
       H = (11*He + 15)/(48*He - 59);
   else
       H = 2.803;
   end
   
   %friction coefficient
   cf = 0.091448*(((H - 1)*ReTheta)^(-0.232))*exp(-1.26*H);
   %dissipation coefficient
   cdiss = 0.010023*((H - 1)*ReTheta)^(-1/6);

   %from momentum and energy integral equations
   dthickdx = zeros(2, 1);
   H
   ue
   dthickdx(1) = cf/2 - (H+2)/ue * duedx * thick(1);
   dthickdx(2) = cdiss - 3/ue * duedx * thick(2);

end