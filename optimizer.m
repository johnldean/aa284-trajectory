order = 6;
delta = 0.001;
Dv_prev = 10000
lr = 1e-10
%coeffs = [-0.0036   -0.0047   -0.0058   -0.0070    0.0934   -0.0034]
coeffs = [0   0   0   0    1   0]
for i = 1:2
    ddv = zeros(1,order);
    for j = 1:length(coeffs)
        dc = zeros(1,order);
        dc(j) = delta;
        ddv(j) = (trajectory_calcs(coeffs + dc) - trajectory_calcs(coeffs))/delta;
    end
    coeffs = coeffs - ddv/delta*lr;
    Dv = trajectory_calcs(coeffs);
    if Dv > Dv_prev
    	lr = lr/2;
    else
    	lr = lr*1.1;
    end
    Dv_prev = Dv;
    [Dv, lr]
end

