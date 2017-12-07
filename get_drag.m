% John Dean, V1, 11-1
% Robbie, V2, 11-5

function drag = get_drag(h,v,Cd,A)
%Gets the drag force on the rocket

rho0 = .015; %nasa data

coeff = -8.3781 * 10^-5; % linearly extrapolated nasa data for logrho vs h from h = 0 to 10000m 

rho = rho0 * exp(h * coeff);

drag = Cd * rho * A * v^2/2;

end
