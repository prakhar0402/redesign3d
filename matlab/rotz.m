function R = rotz(t)
    ct = cosd(t);
    st = sind(t);
    R = [ct st 0; -st ct 0; 0 0 1];