function R = rotx(t)
    ct = cosd(t);
    st = sind(t);
    R = [1 0 0; 0 ct st; 0 -st ct];