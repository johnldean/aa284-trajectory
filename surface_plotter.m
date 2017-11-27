res1 = 2;
res2 = 2;
dv = zeros(res1,res2);
a1_ = linspace(0.08,0.11,res1);
a2_ = linspace(-0.020,0.0,res2);
i = 1;
dv_ = zeros(res1,res2,4);
for steps = [1e3,1e4,1e5,1e6]
    ii = 1;
    for a1 = a1_
        jj=1;
        for a2 = a2_
            v = trajectory_calcs([a2,a1,0],steps);
            if v > 8000
                v = 8000;
            end
            dv(ii,jj) = v;
            jj = jj+1;
        end
        ii = ii +1;
        ii
    end
    dv_(:,:,i) = dv;
    i = i+1;
end
%%

[X, Y] = meshgrid(a2_,a1_);
clf()

for i = 1:4
surf(Y,X,dv_(:,:,i))
hold on
alpha 0.5

end
%plot3(coeffs_(2,:),coeffs_(1,:),Dv_,'r','linewidth',2)
%xlabel('a1')
%ylabel('a2')
%zlabel('DeltaV (m/s)')