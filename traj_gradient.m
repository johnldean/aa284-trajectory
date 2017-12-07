function grad = traj_gradient(coeffs, delta,steps)
    order = length(coeffs)-1;
    grad = zeros(1,order+1);
    
    parfor j = 1:length(coeffs)
        dc = zeros(1,order+1);
        dc(j) = delta;
        grad(j) = (trajectory_calcs(coeffs + dc,steps) - trajectory_calcs(coeffs-dc,steps));
    end
    grad(end) = 0;
    grad = grad/norm(grad);
end