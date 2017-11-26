order = 6;
delta = 0.01;

coeffs = [-0.0036   -0.0047   -0.0058   -0.0070    0.0934   -0.0034];
coeffs = [0  0  0  0  1  0];
for i = 1:100
    ddv = zeros(1,order);
    for j = 1:length(coeffs)
        dc = zeros(1,order);
        dc(j) = delta;
        ddv(j) = (trajectory_calcs(coeffs + dc) - trajectory_calcs(coeffs))/delta;
    end
    trajectory_calcs(coeffs + dc)
    coeffs = coeffs - ddv*delta*0.000004
end