%% Optimal Bounded Low-Thrust Rendezvous with Fixed Terminal-Approach Direction
% State transition approach for unbounded thrust 3 dimensional case 

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
% r0 = [15000; 0; 2000]; v0 = [-10; 0; -2]/norm; x0 = [r0; v0];
r0 = [15000; 0; 0]; v0 = [-10; 0; 0]/norm; x0 = [r0; v0];
r1 = [-300; 0; 0]; v1 = [0.2; 0; 0]/norm;  x1 = [r1; v1];
rf = [-1e-8; 0; 0]; vf = [1e-8; 0; 0]; xf = [rf; vf];

% bounded thrust
Gamma_max = 5*1e-4/norm^2;


% step 1.
alpha = 0;

% costate initial guess
lamb0 = UBT(alpha, t0, t1, x0, x1);
lamb_r0 = lamb0(1:3); lamb_v0 = lamb0(4:6);

itr = 1;
iter_num = 100;
lamb0_hist = zeros(iter_num, 6);
err_hist = zeros(iter_num);

while itr <= iter_num
    if itr == 1

        [time_1, x] = ode45(@(t, x) BT_dyn(t, x, Gamma_max, lamb0), [t0 t1], x0);
        x1_ode = x(end, 1:6);
    end

    lamb0_hist(itr,:) = lamb0';
    
    % define error function
    L = vecnorm(x1_ode'-x1);
    err_hist(itr) = L;

    % calculate 1st order gradient 
    eps = 1e-5;
    del_lamb = diag(lamb0.*eps);
    lamb0_del = zeros(6,6);
    x1_del = zeros(6,6);
    L_del = zeros(6,1);
    dLdlamb = zeros(6,1);
    for i = [1,2,4,5]
        lamb0_del(:, i) = lamb0 + del_lamb(:,i);
        [~, x_del_rx] = ode45(@(t, x) BT_dyn(t, x, Gamma_max, lamb0_del(:,1)), [t0 t1], x0);
        [~, x_del_ry] = ode45(@(t, x) BT_dyn(t, x, Gamma_max, lamb0_del(:,2)), [t0 t1], x0);
        [~, x_del_vx] = ode45(@(t, x) BT_dyn(t, x, Gamma_max, lamb0_del(:,4)), [t0 t1], x0);
        [~, x_del_vy] = ode45(@(t, x) BT_dyn(t, x, Gamma_max, lamb0_del(:,5)), [t0 t1], x0);

        x1_del(1, :) = x_del_rx(end, 1:6);
        x1_del(2, :) = x_del_ry(end, 1:6);
        x1_del(4, :) = x_del_vx(end, 1:6);
        x1_del(5, :) = x_del_vy(end, 1:6);
        
        L_del(i) = vecnorm(x1_del(i,:)' - x1);
        dLdlamb(i) = (L_del(i) - L)/del_lamb(i,i);
    end

    if vecnorm(dLdlamb) < 1e-0
        disp('Solution converged')
        break;
    end

    % lamb0 update
    k = 1e-7;
    lamb0_new = lamb0 - k*dLdlamb;
    [~, x_new] = ode45(@(t, x) BT_dyn(t, x, Gamma_max, lamb0_new), [t0 t1], x0);

    x1_ode_new = x_new(end,1:6);
    err = vecnorm(x1_ode_new'-x1);


    lamb0 = lamb0_new;
    x1_ode = x1_ode_new;
    itr = itr+1;

end


figure()
plot(err_hist)

%%
z0_final = [x0; lamb0];
[time_step1, z_step1] = ode45(@(t,z) BT_eqn(t,z,alpha,Gamma_max), [t0 t1], z0_final);

lamb_v = z_step1(:,10:12);
if (vecnorm(lamb_v,2,2) <= Gamma_max )
    gamma = - lamb_v;
else
    gamma = - Gamma_max ./ vecnorm(lamb_v,2,2) .* lamb_v;
end

x_traj_step1 = z_step1(:,1);
y_traj_step1 = z_step1(:,2);
z_traj_step1 = z_step1(:,3);
gamma_x_step1 = norm^2*gamma(:,1);
gamma_y_step1 = norm^2*gamma(:,2);
gamma_z_step1 = norm^2*gamma(:,3);
gamma_step1 = sqrt(gamma_x_step1.^2 + gamma_y_step1.^2);

time_step1 = time_step1/norm; % [sec]
figure()
subplot(1,2,1)
plot(x_traj_step1, y_traj_step1);
grid on

subplot(1,2,2)
plot(time_step1, [gamma_x_step1, gamma_y_step1, gamma_z_step1]);
