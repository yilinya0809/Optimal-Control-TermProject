%% Optimal Bounded Low-Thrust Rendezvous with Fixed Terminal-Approach Direction
% Bounded Thrust case for 2D example

clc; clear; close all;

% cosmic velocity
height = 480; % [km]
GM = 398600.4418; % [km^3/s^2]
R_earth = 6371; % [km]
cosmic_velocity = sqrt(GM /(R_earth + height)); % [km/s]
orbital_period = 2*pi*(R_earth+height) / cosmic_velocity; % [sec]
angular_velocity = 2*pi/orbital_period;

% time normalization 
% tau = angular_velocity * 1sec
% 1 tsec = 899.5542 sec

norm = angular_velocity; 
t0 = 0; t1 = norm*11*orbital_period; 
t2 = norm*1885; tf = t1+t2;

% boundary conditions 
r0 = [15000; 0; 0]; v0 = [-10; 0; 0]/norm; x0 = [r0; v0];
r1 = [-300; 0; 0]; v1 = [0.2; 0; 0]/norm;  x1 = [r1; v1];
rf = [-1e-8; 0; 0]; vf = [1e-8; 0; 0]; xf = [rf; vf];

% Thrust constraint
Gamma_max = 5*1e-4;
gm_norm = 5*1e-4/norm^2;

%% Step 1. Bounded Thrust Gradient method
alpha = 0;

% costate initial guess
lamb0 = UBT(alpha, t0, t1, x0, x1);
lamb0_original = lamb0;
% desired costate from fsolve
% lamb0_d = [5.9927; 2190.2; 0; -1177.9; 25.1969; 0];
lamb0_d = fsolve(@(lamb0) F(lamb0, t0, t1, x0, x1, gm_norm), lamb0);


itr = 1;
iter_num = 45;
lamb0_hist = zeros(iter_num, 6);
err_hist = zeros(iter_num,1);

while itr <= iter_num
    if itr == 1
        [time_1, x] = ode45(@(t, x) BT_dyn(t, x, gm_norm, lamb0), [t0 t1], x0);
        x1_ode = x(end, 1:6);
    end

    lamb0_hist(itr,:) = lamb0';
    
    % define error function
    L = vecnorm(x1_ode'-x1);
    err_hist(itr) = L;

   
    % direction vector
    d = lamb0_original - lamb0_d;

    % update costate
    k = d(1)/100;
    lamb0_up = lamb0 + k*d;

    [~, z_new] = ode45(@(t, x) BT_dyn(t, x, gm_norm, lamb0_up), [t0 t1], x0);

    x1_ode_new = z_new(end,1:6);
    err = vecnorm(x1_ode_new'-x1);

    lamb0 = lamb0_up;
    x1_ode = x1_ode_new;
    itr = itr+1;

end

figure()
plot(err_hist);
lamb0_fs_grad = lamb0;
x1_fs_grad = x1_ode;


%% find lamb0 by grid shooting
itr=1;
itr_num = 20;
lamb0_hist2 = zeros(itr_num,6);
err_hist2 = zeros(itr_num,1);
lamb0 = lamb0_fs_grad;

while itr <= itr_num
    lamb0_hist2(itr,:) = lamb0';
    if itr == 1
        [time_1, z] = ode45(@(t, z) BT_eqn(t, z, alpha, gm_norm), [t0 t1], [x0;lamb0]);
        x1_ode = z(end, 1:6);
        L = vecnorm(x1_ode'-x1);
    end
    
    err_hist2(itr) = L;

    grid = diag(lamb0.*0.2);

    L_new = repmat(L, [2 2 2 2]);
    for a = [1,2]
        rx_grid(:,a) = (-1)^a*grid(:,1);
        for b = [1,2]
            ry_grid(:,b) = (-1)^b*grid(:,2);
            for c= [1,2]
                vx_grid(:,c) = (-1)^c*grid(:,4);
                for d = [1,2]
                    vy_grid(:,d) = (-1)^d*grid(:,5);
                    lamb0_new = lamb0 + rx_grid(:,a) + ry_grid(:,b) + vx_grid(:,c) + vy_grid(:,d);
                    [~, z_new] = ode45(@(t, z) BT_eqn(t, z, alpha, gm_norm), [t0 t1], [x0;lamb0_new]);
                    x1_new = z_new(end,1:6);
                    L_new(a,b,c,d) = vecnorm(x1_new' - x1);
                end
            end
        end

    end

    L_new_rs = reshape(L_new, 16,1);
    [L_min, index] = min(L_new_rs);

    if L_min > L
        break;
    end

    if (fix((index-1)/4) == 0)
        c=1; d=1;
    elseif (fix((index-1)/4) == 1)
        c=2; d=1;
    elseif (fix((index-1)/4) == 2)
        c=1; d=2;
    elseif (fix((index-1)/4) == 3)
        c=2; d=2;
    end

    if(mod(index, 4) == 1)
        a=1; b=1;
    elseif(mod(index,4) == 2)
        a=2; b=1;
    elseif(mod(index,4) == 3)
        a=1; b=2;
    elseif(mod(index,4) == 0)
        a=2; b=2;
    end

    lamb0 = lamb0 + rx_grid(:,a) + ry_grid(:,b) + vx_grid(:,c) + vy_grid(:,d);
    L = L_min;  
    itr = itr+1;
    index_hist(itr) = index;

end
lamb0_grid = lamb0;
[~, x_grid] = ode45(@(t,z) BT_eqn(t, z, alpha, gm_norm), [t0 t1], [x0;lamb0_grid]);
x1_grid = x_grid(end,1:6);

figure()
plot(err_hist2);

%% after finding lamb0
lamb0_final = lamb0_fs_grad;
z0_final = [x0; lamb0_final];
[time_step1, z_step1] = ode45(@(t,z) BT_eqn(t,z,alpha,gm_norm), [t0 t1], z0_final);
x1_step1 = z_step1(end,1:6);

x_traj_step1 = z_step1(:,1);
y_traj_step1 = z_step1(:,2);
z_traj_step1 = z_step1(:,3);
lamb_v_step1 = z_step1(:,10:12);
gamma_x_step1 = -norm^2*lamb_v_step1(:,1);
gamma_y_step1 = -norm^2*lamb_v_step1(:,2);
gamma_z_step1 = -norm^2*lamb_v_step1(:,3);
gamma_step1 = sqrt(gamma_x_step1.^2 + gamma_y_step1.^2);
for i = 1:length(gamma_step1)
    if (gamma_step1(i) > Gamma_max)
        gamma_step1(i) = Gamma_max;
    end
end



time_step1 = time_step1/norm/3600; % [hour]

figure()
plot(time_step1, gamma_step1);
hold on
Gamma_max_traj = Gamma_max*ones(length(time_step1));
plot(time_step1, Gamma_max_traj);
ylim([0 6*1e-4]);
ylabel('|\Gamma| [m/sec^2]'); xlabel('Time [hours]');
title('First-stage optimal thrust acceleration')

%% Step 2. Unsaturated
alpha = 5000;

lamb1 = UBT(alpha, t1, tf, x1, xf);
z1 = double([x1; lamb1]);

[time_step2, z_step2] = ode45(@(t,z) UBT_eqn(t,z,alpha), [t1 tf], z1);

x_traj_step2 = z_step2(:,1);
y_traj_step2 = z_step2(:,2);
z_traj_step2 = z_step2(:,3);
gamma_x_step2 = -norm^2*z_step2(:,10);
gamma_y_step2 = -norm^2*z_step2(:,11);
gamma_z_step2 = -norm^2*z_step2(:,12);
gamma_step2 = sqrt(gamma_x_step2.^2 + gamma_y_step2.^2);

time_step2 = (time_step2 - time_step2(1))/norm; % [sec]

figure()
plot(time_step2, gamma_step2);
ylim([0 6*1e-4]);
ylabel('|\Gamma| [m/sec^2]'); xlabel('Time [hours]');
title('Second-stage optimal thrust acceleration')
grid on

%% Plot trajectories
figure()
plot(x_traj_step1, y_traj_step1);
grid on
hold on
plot(x_traj_step2, y_traj_step2);
