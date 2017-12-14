% John Dean, V1, 11-1
% Robbie, V2, 11-5

function drag = get_drag(h,v,Cd,A)
%Gets the drag force on the rocket
y = 1.289;
R = 189;
if h > 10000
    rho0 = .015; %nasa data

    coeff = -8.3781 * 10^-5; % linearly extrapolated nasa data for logrho vs h from h = 0 to 10000m 

    rho = rho0 * exp(h * coeff);
else
    scaleHeight = 0.00009; % mars []
    p = 0.699 * exp(-scaleHeight * h);
    if h < 7000
        T = -31 - 0.000998 * h; % [C]
    else
        T = -23.4 - 0.00222 * h;
    end
    rho = p/(0.1921 * (T + 273.15));
end
if h < 7000
    T = -31 - 0.000998 * h; % [C]
else
    T = -23.4 - 0.00222 * h;
end
M = sqrt(y * R * T);
if T <= 0
    Cd = .11;
elseif M < .85
    Cd = .2;
else
    Cd = .11 + .82 / M^2 - .55 /  M^4;
end
drag = Cd * rho * A * v^2/2 ;
end
