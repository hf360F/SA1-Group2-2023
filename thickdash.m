function dthickdx = thickdash(xmx0, thick)
   global ReL ue0 duedx;
   
   ReTheta = ReL * ue0 * thick(1);
   He = thick(2)/thick(1);
   
   if He >= 1.46
       H = (11*He + 15)/(48*He - 59);
   else
       H = 2.803;
   end

   cf = 0.09148*(((H - 1)*ReTheta)^(-0.232))*exp(-1.26*H);
   cdiss = 0.010023*((H - 1)*ReTheta)^(-1/6);

    dthickdx = zeros(2, 1);
    dthickdx(1) = cf/2 - (H+2)/ue0 * duedx * thick(1);
    dthickdx(2) = cdiss - 3/ue0 * duedx * thick(2);

end