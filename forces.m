% calculate list and drag

function [cl, cd] = forces(circ,cp,delstarl,thetal,delstaru,thetau)

    cl = -2*circ;

    thetate = thetal(end) + thetau(end);
    delstarte = delstarl(end) + delstaru(end);
    Hte = delstarte/thetate;
    uete = sqrt(1 - cp(end));

    thetainf = thetate*uete^((Hte + 5)/2);
    cd = 2*thetainf;
end
