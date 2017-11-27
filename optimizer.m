%%
minDv_ = [];
steps = 1e4;
coeffs = [0.1, 0];
for order = 2:10
    coeffs = [0 , coeffs];
    %coeffs = [-0.0036   -0.0047   -0.0058   -0.0070    0.0934   -0.0034]
    delta = 0.0005;
    lr = 1e-3;
    Dv_prev = trajectory_calcs(coeffs,steps);
    grad = traj_gradient(coeffs,delta,steps);
    Dv_ = [];
    coeffs_ = [];
    lr_ = [];
    while lr > 1e-8
        coeffs_new = coeffs - grad*lr;
        Dv_new = trajectory_calcs(coeffs_new,steps);
        if Dv_new >= Dv_prev
            lr = lr/2;
        else
            lr = lr*1.1;
            coeffs = coeffs_new;
            Dv_prev = Dv_new;
            grad = traj_gradient(coeffs,delta,steps);
        end
        coeffs;

        {lr, order}
        coeffs_ = [coeffs_,coeffs'];
        Dv_ = [Dv_, Dv_prev];
        lr_ = [lr_, lr];
    end
    order
    Dv_prev
    minDv_ = [minDv_, Dv_prev];
end
plot(2:10,minDv_)
%%
    subplot(121)
    plot(Dv_/1000)
    ylabel('Delta V (km/s)')
    xlabel('iteration number')
    subplot(122)
    semilogy(lr_)
    ylabel('learning rate')
    xlabel('iteration number')
    