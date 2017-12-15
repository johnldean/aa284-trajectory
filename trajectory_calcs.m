%John Dean, V1, 11-3
%Robbie, V2, 11-5
%James, V3, 11-6
%Robbie, V4, 12-3, Drag after cutoff
%clear all
%Inputs are thrust(given from evaluating thrust eqn. at time vector) and
%maxq = 4.4418 kPa
function deltaV = trajectory_calcs(p_coeffs,steps)%inputpolycoeffs)
steps;
p_coeffs
% Mars characteristics
r_mars = 6779000/2; %radius mars [m]
m_mars = 6.39E23; %mass mars [kg]
u_mars = 4.282837E13; %standard gravitational parameter [m^3/s^2]

alt_final = 210000; %final height [m]
orbit_inject = sqrt(u_mars /(alt_final + r_mars));

v_final = orbit_inject;

TimeMax = 1000;% [s]
t = linspace(0,TimeMax,steps);
dT = TimeMax/steps;
M_payl = 61.4;
M2_prop = 64.7815;
M2_struc = 8.8338;
M2_t = M_payl + M2_prop + M2_struc;
M1_prop = 118.2539;
M1_struc = 17.6701;
M1_t = M2_t + M1_prop + M1_struc;

isp1 = 316.41 * 9.81;%* u_mars / r_mars^2; %m/s
isp2 = 316.41 * 9.81;%* u_mars / r_mars^2; %m/s
thrust1 = 9400;
thrust2 = 1400;
bTime1 = M1_prop / (thrust1 / isp1);
bSteps1 = floor(bTime1 / dT);
bTime2 = M2_prop / (thrust2 / isp2);
bSteps2 = floor(bTime2 / dT);
timeSteps = floor(TimeMax / dT);
stageTime = 2;
stageSteps = floor(stageTime / dT);

steps_stage = bSteps1 + stageSteps;
T1 = thrust1 * ones(bSteps1, 1);
Tstaging = zeros(stageSteps, 1);
T2 = thrust2 * ones(timeSteps - (bSteps1 + stageSteps), 1);
thrust_ = vertcat(T1, Tstaging, T2);
I1 = isp1 * ones(bSteps1, 1);
Istaging = isp1 * ones(stageSteps, 1);
I2 = isp2 * ones(timeSteps - (bSteps1 + stageSteps), 1);
isp_ = vertcat(I1, Istaging, I2);

%thrust_ = 9400 * ones(1, length(t));%thrust eq, eval at time %(see graph)/cantwell 283 eqns-revise

mars_rot = 241.17; %m/s
deg2rad = 2 * pi / 360;
i_init = 30;
v_rot = mars_rot * cos(deg2rad * i_init);
%inital conditions
r_vel_init = 0;
theta_vel_init = v_rot;

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
A = pi / 4 * .64^2;
%M = 350 * ones(1, length(t)); %mass rocket [kg] %losing mass. ADD IN
M = zeros(TimeMax, 1);
Minit = M1_t;
M(1) = Minit;
%ut = ti


%note: make this shit another function. implement RK4/trapezoidal in long
%run

TV_th = zeros(1,length(t)); %trust vector angle
TV_th_init = pi/2; % takeoff at vertical [rad]
TV_th(1) = TV_th_init;
cutoff_time = t(end);
drag_loss = zeros(TimeMax, 1);
grav_loss = zeros(TimeMax, 1);
TVC_loss = zeros(TimeMax, 1);
for i = 2: length(t)
    %i
    thrust = thrust_(i - 1);
%     M(i) = M(i - 1) - thrust / isp_(i - 1) * dT;
%     M_curr = M(i);
%     if M(i) < M2_t - M2_prop
%         M(i) = M2_t - M2_prop;
%         %output = 'burnout'
%     end
%     if i == steps_stage + 1
%         output = 'staging'
%         hehehe = sum(delta_v);
%         M(i) = M2_t;
%     end
    %TV_th_act = TV_th(i-1);
    r_act = r_(i-1);
    alt_act = r_act - r_mars;  %from polynomial..what we want
    theta_vel_act = theta_vel(i-1);
    r_vel_act = r_vel(i-1);
    %alt_targ_norm = polyval(p_coeffs,theta_vel_norm); % in time this will be polyval of function at the v_theta value

    if theta_vel_act ~= 0
        veh_angle = atan(r_vel_act / theta_vel_act);
    else
        veh_angle = pi/2;
    end
    %t_curr = i * dT;
    t_norm = i / (bSteps1 + bSteps2 + stageSteps);
    inv_TV_new_norm = (polyval(p_coeffs, t_norm));
    TV_th_new = pi / 2 * (1 - inv_TV_new_norm);
    %TV_th_new = veh_angle + FACtor * alt_diff_norm^5;% * dT
    if TV_th_new < 0
        TV_th_new = 0;
    end
    if TV_th_new > pi/2
        TV_th_new = pi/2;
    end
    %velocity vector
    %v_mag = sqrt(( theta_vel_act - v_rot * r_act / r_mars)^2 + r_vel_act^2);  %burn happens at about same latitude as launch
    v_mag = sqrt(( theta_vel_act - v_rot)^2 + r_vel_act^2);
    %delta theta ~.06radians
    Fd_norm = get_drag(alt_act,norm(v_mag),Cd,A); %drag
    if v_mag > 0
        Drag_r = Fd_norm * r_vel_act / v_mag;
        %Drag_theta = Fd_norm * (theta_vel_act - v_rot * r_act / r_mars) / v_mag;
        Drag_theta = Fd_norm * (theta_vel_act - v_rot) / v_mag;
    else
        Drag_r = 0;
        Drag_theta = 0;
    end
    TVC_angle = TV_th_new - veh_angle;
    TVC_loss(i) = thrust * sin(abs(TVC_angle)) * dT / M(i - 1);
    %if thrust < 1500
        %TV_th_new = 1*veh_angle;
    %end
    F_g = M(i - 1) * u_mars/r_act^2;     %gravity..a;lways r-dir; F = ma = GMm/r^2 => F = u * m/r^2
    r_accel = (thrust * sin(TV_th_new) - Drag_r - F_g)/M(i - 1); %%%
    theta_accel = (thrust * cos(TV_th_new) - Drag_theta)/M(i - 1); %%%
    grav_loss(i) = F_g / M(i - 1) * dT * sin(TV_th_new); %%%
    if thrust == 0
        grav_loss(i) = 0;
    end
    drag_loss(i) = Fd_norm / M(i - 1) * dT; %%%
    theta_vel_new = theta_vel_act +  theta_accel * dT;
    r_vel_new = r_vel_act +  r_accel * dT;
    
    r_new = r_act + dT * r_vel_act + r_accel * dT^2 / 2;
    theta_new = theta(i-1) + theta_vel_act / r_act * dT + theta_accel /2 * dT^2 / r_act;
    
    theta(i) = theta_new;
    TV_th(i) = TV_th_new;
    del_theta = -atan(theta_vel_new * dT / r_act);
    r_(i) = r_new;
    %del_theta = pi / 3;
    %r_vel_rotate = r_vel_new * cos(del_theta) + theta_vel_new * sin(del_theta);
    %theta_vel_rotate = -r_vel_new * sin(del_theta) + theta_vel_new * cos(del_theta);
    point1 = [r_new, 0];
    point2 = [r_new + r_vel_new, theta_vel_new];
    rotation_mat = [cos(del_theta), sin(del_theta); -sin(del_theta), cos(del_theta)];
    point1_rot = point1 * rotation_mat;
    point2_rot = point2 * rotation_mat;
    rot_vel = point2_rot - point1_rot;
    r_vel_rotate = rot_vel(1);
    theta_vel_rotate = rot_vel(2);
    theta_vel(i) = theta_vel_rotate;
    r_vel(i) = r_vel_rotate;
    total_v = sqrt(r_vel_rotate^2 + theta_vel_rotate^2);
    energy = total_v^2 / 2 - u_mars / r_new;
    semi_major = -u_mars / (2 * energy);
    eccentricity = total_v^2 * r_new / u_mars * [1,0] - ...
        r_new * r_vel_rotate * total_v / u_mars * [r_vel_rotate, theta_vel_rotate] / total_v - ...
        [1, 0];
    alt_ap = semi_major * ( 1 + norm(eccentricity)) - r_mars;
    delta_v(i) = thrust / M(i - 1) * dT;
    M(i) = M(i - 1) - thrust / isp_(i - 1) * dT;
    M_curr = M(i);
    if M(i) < M2_t - M2_prop
        M(i) = M2_t - M2_prop;
        %output = 'burnout'
    end
    if i == steps_stage + 1
        output = 'staging';
        hehehe = sum(delta_v);
        M(i) = M2_t;
    end
    res = 2;
    addBurn = zeros(1,res);%additional burn after cutoff to account for drag loss
    delta_v_addBurn = [];
    j = i;
    if alt_ap > alt_final % this is not accounting for drag after cutoff
        %sum(delta_v)
        note = 'arrived';
        sum(delta_v);
        for h = 1:res
            cutoff_time = i + sum(addBurn);
            %delta_v_cutoff = sum(delta_v) + sum(delta_v_addBurn)
            delta_v_toAlt = [];
            delta_V_cutToFinAlt = 0;
            for j = cutoff_time + 1:length(t)
                thrust = 0;
                r_act = r_(j-1);
                alt_act = r_act - r_mars;  %from polynomial..what we want
                theta_vel_act = theta_vel(j-1);
                r_vel_act = r_vel(j-1);
                v_mag = sqrt(theta_vel_act^2 + r_vel_act^2);                                     %altitude
                Fd_norm = get_drag(alt_act,norm(v_mag),Cd,A); %drag
                drag_loss(i + j) = Fd_norm / M(i) * dT;
                drag_dv = Fd_norm / M(i) * dT;
                if v_mag > 0
                    Drag_r = Fd_norm * r_vel_act / v_mag;
                    Drag_theta = Fd_norm * theta_vel_act / v_mag;
                else
                    Drag_r = 0;
                    Drag_theta = 0;
                end
                F_g = M(i) * u_mars/r_act^2;     %gravity..a;lways r-dir; F = ma = GMm/r^2 => F = u/r^2
                
                %accel
                r_accel = (thrust * sin(TV_th_new) - Drag_r - F_g)/M(i);
                theta_accel = (thrust * cos(TV_th_new) - Drag_theta)/M(i);
                %vel
                theta_vel_new = theta_vel_act +  theta_accel * dT;
                r_vel_new = r_vel_act +  r_accel * dT;
                %posn
                
                r_new = r_act + dT * r_vel_act + r_accel * dT^2 / 2;
                
                theta_new = theta(j-1) + theta_vel_act / r_act * dT + theta_accel /2 * dT^2 / r_act;
                
                %assign to vectors for the timestep
                r_(j) = r_new;
                theta(j) = theta_new;
                TV_th(j) = TV_th_new;
                %
                del_theta = atan(theta_vel_new * dT / r_act);
                r_vel_rotate = r_vel_new * cos(del_theta) + theta_vel_new * sin(del_theta);
                %eheh = heh
                theta_vel_rotate = -r_vel_new * sin(del_theta) + theta_vel_new * cos(del_theta);
                theta_vel(j) = theta_vel_rotate;
                r_vel(j) = r_vel_rotate;
                if  r_vel_new < 0
                    r_new;
                    delta_V_cutToFinAlt = sum(delta_v_toAlt);
                    break;
                end
                delta_v_toAlt(j-i+1) =  drag_dv; %accel mag * dt = delvmag
            end
            delta_V_cutToFinAlt;
            %delta_v_drag = delta_V_final - delta_v_cutoff;
            delV_addBurn = [];
            mass_start = M(i);
            mass_new = mass_start;
            while  sum(delV_addBurn) <  delta_V_cutToFinAlt
                addBurn(h) = addBurn(h) + 1;
                delV_addBurn(addBurn(h)) = thrust_(i - 1 + addBurn(h)) / mass_new * dT;
                mass_new = mass_new - thrust_(i - 1 + addBurn(h)) / isp_(i - 1 + addBurn(h)) * dT;
                if mass_new < M2_t - M2_prop
                    mass_new = M2_t - M2_prop;
                end
            end
            
            delta_v_addBurn(h) = sum(delV_addBurn);
        end

        break;
    end

end

if 1
    cut = (j);
    subplot(511)
    length(t);
    length(r_);
    plot(t(1:cut), (r_(1:cut) - r_mars)/1000) %alt in km
    xlabel('t (s)');
    ylabel('alt (km)');
    subplot(512)
    plot(t(1:cut), theta_vel(1:cut)) %alt in km
    xlabel('t (s)');
    ylabel('theta vel (m/s');
    subplot(513)
    plot(theta_vel(1:cut), (r_(1:cut) - r_mars)/1000)
    xlabel('theta vel (m/s)');
    ylabel('alt (km)');
    subplot(514)
    plot(t(1:cut), r_vel(1:cut))
    xlabel('t (s)');
    ylabel('r vel (m/s)');
    subplot(515)
    plot(t(1:cut), TV_th(1:cut))
    xlabel('t (s)');
    ylabel('TV (rad)');
    drawnow()
end

drag_loss_t = sum(drag_loss)
grav_loss_t = sum(grav_loss)
TVC_loss_t = sum(TVC_loss);
TVC_loss;
cutoff_time;
DEL_V_DRAG_COMP = sum(delta_v_addBurn);
PERC_Del_V = DEL_V_DRAG_COMP/sum(delta_v) * 100;
orbit_inject;
theta_vel(j);
deltaV = sum(delta_v)+ sum(delta_v_addBurn)+ orbit_inject - theta_vel(j)
orbit_inject - theta_vel(j)
del_v_t = sum(delta_v);
sum(delta_v_addBurn);
p_coeffs;
end