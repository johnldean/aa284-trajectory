%%
minDv_ = [];
delta = 5e-4;%/(2^i);
steps = 5e4;
coeffs = [.5];
data = cell(1,9);

for order = 2:5%10
    coeffs = [0 , coeffs];
    %coeffs = [-0.0957   -0.1129   -0.0426    0.5026    0.6756]
    Dv_prev = trajectory_calcs(coeffs,steps);
    %aaa = bbb
    grad = traj_gradient(coeffs,delta,steps);
    Dv_ = [];
    coeffs_ = [];
    lr_ = [];
    for i = 0:0
        lr = 1e-0;
        delta = 5e-4/(2^i);
        steps = 1e4*2^i;
        grad = traj_gradient(coeffs,delta,steps);
        %aaa = bbb
        while lr > 1e-8
            coeffs_new = coeffs - grad*lr;
            Dv_new = trajectory_calcs(coeffs_new,steps);
            %if Dv_prev - Dv_new < .0001 && Dv_new < Dv_prev
            %    break;
            %end
            if Dv_new >= Dv_prev
                %lr = lr*1.1;
                lr = lr / 1.1;
            else
                %lr = lr/2;
                lr = lr * 1.5;
                coeffs = coeffs_new;
                Dv_prev = Dv_new;
                grad = traj_gradient(coeffs,delta,steps);
            end
            lr
            {lr, order,i}
            coeffs_ = [coeffs_,coeffs'];
            Dv_ = [Dv_, Dv_prev];
            lr_ = [lr_, lr];
        end
    end
    order
    coeffs
    Dv_prev
    minDv_ = [minDv_, Dv_prev];
    data(order-1) = {[lr_ ;Dv_]};
end
%%
clf()
plot(2:10,minDv_)
xlabel('order')
ylabel('delta v (m/s)')
%%

for i = 1:9
    subplot(9,2,(i-1)*2+1)
    lr_ = data{i}(1,:);
    Dv_ = data{i}(2,:);
    plot(Dv_/1000)
    ylabel('Delta V (km/s)')
    xlabel('iteration number')
    subplot(9,2,(i-1)*2+2)
    semilogy(lr_)
    ylabel('learning rate')
    xlabel('iteration number')
end    

%%
step_ = [];
dv_ = [];
for i = 3:7
    step_ = [step_, 10^i];
    dv_ = [dv_, trajectory_calcs(coeffs,10^i)];
    i
end
clf()
semilogx(step_,dv_)
xlabel('number of steps')
ylabel('delta v (m/s)')