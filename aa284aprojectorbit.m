u_mars = 4.282837 * 10^13;
r0 = 3389500; %m
i_init = 30; %deg
i_final = 0;
h_init = 0;
h_final = 200000 %m
mars_rot = 241.17; %m/s
deg2rad = 2 * pi / 360;

V_co = sqrt(u_mars / (r0 + h_final)) %km/s



V_rot = mars_rot * cos(deg2rad * i_init) %m/s



deltavhoh = (2 * u_mars * (1 / r0 - 1 / ( r0 + r0 + h_final) ) ) ^ .5 - V_rot - ...
    (2 * u_mars * (1 / (r0 + h_final) - 1 / ( r0 + r0 + h_final) ) ) ^ .5 + (u_mars / (r0 + h_final)) ^ .5
V_incchange = 2 * V_co * sin(.5 * deg2rad * abs(i_final - i_init)) %m/s