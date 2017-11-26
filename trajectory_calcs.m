%John Dean, V1, 11-3
%Robbie, V2, 11-5
%James, V3, 11-6
%clear all
%Inputs are thrust(given from evaluating thrust eqn. at time vector) and
function deltaV = trajectory_calcs(p_coeffs)%inputpolycoeffs)

% Mars characteristics
r_mars = 6779000/2; %radius mars [m]
m_mars = 6.39E23; %mass mars [kg]
u_mars = 4.282837E13; %standard gravitational parameter [m^3/s^2]

alt_final = 600000; %final height [km]
orbit_inject = sqrt(u_mars /(alt_final + r_mars));

v_final = orbit_inject;

bTimeMax = 100*2; % [s]
steps = 1000000;
t = linspace(0,bTimeMax,steps);
dT = bTimeMax/steps; 
thrust_ = 10000 * ones(1, length(t));%thrust eq, eval at time %(see graph)/cantwell 283 eqns-revise

%inital conditions
r_vel_init = 0;
theta_vel_init = 0;

%pre-allocate

delta_v = zeros(1,length(t));
%velocity
r_vel = zeros(1,length(t)); %y-vel
r_vel(1) = r_vel_init;

theta_vel = zeros(1,length(t));%x-vel
theta_vel(1) = theta_vel_init;

%position
r_ = zeros(1,length(t));
r_(1) = r_mars;
theta = zeros(1,length(t));



Cd = 2; %limit is 2..see eqn from slides
A = .1;
%M = 350 * ones(1, length(t)); %mass rocket [kg] %losing mass. ADD IN
M = 350 * (1 - .95 * t / bTimeMax);


%note: make this shit another function. implement RK4/trapezoidal in long
%run

TV_th = zeros(1,length(t)); %trust vector angle
TV_th_init = pi/2; % takeoff at vertical [rad]
TV_th(1) = TV_th_init;
for i = 2: length(t)
    thrust = thrust_(i);
    TV_th_act = TV_th(i-1);
    r_act = r_(i-1);
    alt_act = r_act - r_mars;  %from polynomial..what we want
    theta_vel_act = theta_vel(i-1);
    theta_vel_norm = theta_vel_act / orbit_inject;
    r_vel_act = r_vel(i-1);
    alt_targ_norm = polyval(p_coeffs,theta_vel_norm); % in time this will be polyval of function at the v_theta value
    %alt_targ_norm = (1 - exp(-5 * theta_vel_norm));
    alt_targ = alt_targ_norm * alt_final;
    
    %set lmits of 90 and 0..quadrant so don;t correct over 90
    %factor needs to be more than a first order controller..take into
    %account the slope of last correction
    alt_diff_lim = 1000;
    alt_diff = alt_targ - alt_act;
    if alt_diff > alt_diff_lim
        alt_diff_norm = alt_diff_lim; % if alt. diff is above the limit, set to the limit
    elseif alt_diff < -alt_diff_lim
        alt_diff_norm = -alt_diff_lim;
    else
        alt_diff_norm = alt_diff;
    end
    alt_diff_norm = alt_diff_norm / alt_diff_lim;
    FACtor = 5/40; %theta changes by maximum ~1deg
    TV_th_new = TV_th_act + FACtor * alt_diff_norm;
    if TV_th_new < 0%-1
        TV_th_new = 0;%-1;
    end
    if TV_th_new > 1.57%2.57
        TV_th_new = 1.57;%2.57;
    end
    
    %d = (r^2 + r_prev^2 - 2*r*r_prev*cos(th-th_prev))^.5;   %distance between points
    %velocity vector
    v_mag = sqrt(theta_vel_act^2 + r_vel_act^2);                                     %altitude
    Fd_norm = get_drag(alt_act,norm(v_mag),Cd,A); %drag
    if v_mag > 0
        Drag_r = Fd_norm * r_vel_act / v_mag;
        Drag_theta = Fd_norm * theta_vel_act / v_mag;
    else
        Drag_r = 0;
        Drag_theta = 0;
    end
    
    F_g = M(i) * u_mars/r_act^2;     %gravity..a;lways r-dir; F = ma = GMm/r^2 => F = u/r^2
    
    r_accel = (thrust * sin(TV_th_new) - Drag_r - F_g)/M(i);
    theta_accel = (thrust * cos(TV_th_new) - Drag_theta)/M(i);
    
    theta_vel_new = theta_vel_act +  theta_accel * dT;
    r_vel_new = r_vel_act +  r_accel * dT;
    
    r_new = r_act + dT * r_vel_act + r_accel * dT^2 / 2;
    theta_new = theta(i-1) + theta_vel_act / r_act * dT + theta_accel /2 * dT^2 / r_act;
    r_(i) = r_new;
    theta(i) = theta_new;
    TV_th(i) = TV_th_new;
    del_theta = theta_vel_new / r_act * dT;
    r_vel_rotate = r_vel_new * cos(del_theta) + theta_vel_new * sin(del_theta);
    theta_vel_rotate = r_vel_new * sin(del_theta) + theta_vel_new * cos(del_theta);
    theta_vel(i) = theta_vel_rotate;
    r_vel(i) = r_vel_rotate;
    total_v = sqrt(r_vel_rotate^2 + theta_vel_rotate^2);
    energy = total_v^2 / 2 - u_mars / r_new;
    semi_major = -u_mars / (2 * energy);
    eccentricity = total_v^2 * r_new / u_mars * [1,0] - ...
        r_new * r_vel_rotate * total_v / u_mars * [r_vel_rotate, theta_vel_rotate] / total_v - ...
        [1, 0];
    alt_ap = semi_major * ( 1 + norm(eccentricity)) - r_mars;
    delta_v(i) = thrust / M(i) * dT;
    if alt_ap - 1000 > alt_final
        cutoff_time = i;
        break;
    end
end

subplot(311)
plot(t(1:cutoff_time), (r_(1:cutoff_time) - r_mars)/1000) %alt in km
xlabel('t (s)');
ylabel('alt (km');
subplot(312)
plot(t(1:cutoff_time), theta_vel(1:cutoff_time)) %alt in km
xlabel('t (s)');
ylabel('theta vel (km');
subplot(313)
plot(theta_vel(1:cutoff_time), (r_(1:cutoff_time) - r_mars)/1000)
xlabel('theta vel (m/s)');
ylabel('alt (km');
%subplot(414)
%sx = [-1, -1, -1, 1, 1, 1]';% + 2;
%sx2 = [1, -1, -1, 1, -1, -1, 1]';% + 5;
%sx3 = [1, 0, -1, 0, 0]' + 8;
%ys = [1, -1, 0, 0, 1, -1]';% + 2;
%ys2 = [1, 1, 0, 0, 0, -1, -1]';% + 2;
%ys3 = [1, 0, 1, 0, -1]' + 2;
%reps = 10;
%for i = 1:reps
%    plot(sx + 2 + 8 * (i-1), ys + 5, 'DisplayName', 'H');
%    hold on;
%    plot(sx2 + 6 + 8 * (i-1), ys2 + 5, 'DisplayName', 'E');
%    hold on;
%end
%legend('show');
%xlim([.5 .5 + 8 * reps])
%ylim([.5, 9.5]);
%plot(theta_vel(1:cutoff_time) / orbit_inject, (1 - exp(-5 * theta_vel(1:cutoff_time)/orbit_inject) ));
delta_v_approx = sum(delta_v) + orbit_inject - max(theta_vel);
deltaV = delta_v_approx;
end