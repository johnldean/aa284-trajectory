% John Dean, V1, 11-1
% Robbie, V2, 11-5

function drag = get_drag(h,v,Cd,A)

%Gets the drag force on the rocket
scaleHeight = 0.00009; % mars []
p = 0.699 * exp(-scaleHeight * h);
% if h < 7000
%     T = -31 - 0.000998*h; % [C]
% else
%     T = -23.4 - 0.00222*h;
% end
T = -31; %surface temp [C]
rho = p/(0.1921 * (T + 273.1)); % change to Kelvin

drag = Cd * rho * A * v^2/2;

end
