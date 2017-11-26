format long
order = 6;
delta = 0.000001;
lr = 1e-7;
coeffs = [-0.0036   -0.0047   -0.0058   -0.0070    0.0934   -0.0034]
coeffs = [0   0   0   0    0.2   0]
Dv_prev = trajectory_calcs(coeffs);

%get best first order tragectory

for i = 1:50
    ddv = zeros(1,order);
    for j = 1:length(coeffs)
        dc = zeros(1,order);
        dc(j) = delta;
        ddv(j) = (trajectory_calcs(coeffs + dc) - trajectory_calcs(coeffs))/delta;
    end
    coeffs_new = coeffs - ddv*lr;
    Dv_new = trajectory_calcs(coeffs_new);
    if Dv_new >= Dv_prev
    	lr = lr/2;
    else
    	lr = lr*1.1;
        coeffs = coeffs_new;
        Dv_prev = Dv_new;
    end
    coeffs_new
    coeffs
    {Dv_prev, lr}
    
end

