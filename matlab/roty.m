function R = roty(t)
    ct = cosd(t);
    st = sind(t);
    R = [ct 0 -st; 0 1 0; st 0 ct];