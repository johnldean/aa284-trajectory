r_mars = 6779000/2; %radius mars [m]
m_mars = 6.39E23; %mass mars [kg]
u_mars = 4.282837E13; %standard gravitational parameter [m^3/s^2]
h = 210000;
v = (u_mars/(h+r_mars))^0.5;
A = pi*0.75^2;
get_drag(h,v,2,1)
