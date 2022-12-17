%% Optimal Bounded Low-Thrust Rendezvous with Fixed Terminal-Approach Direction
% Unbounded Thrust case for 3D example

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
t0 = 0; t1 = norm*2*orbital_period; 
t2 = norm*1400; tf = t1+t2;

% boundary conditions 
r0 = [15000; 0; 0]; v0 = [-10; 0; -2]/norm; x0 = [r0; v0];
r1 = [-300; 0; 0]; v1 = [0.35; 0; 0]/norm;  x1 = [r1; v1];
rf = [-1e-8; 0; 0]; vf = [1e-8; 0; 0]; xf = [rf; vf];

% Thrust constraint
Gamma_max = 5*1e-4;
gm_norm = 5*1e-4/norm^2;

%% Step 1. 
alpha = 0;

lamb0 = UBT(alpha, t0, t1, x0, x1);

z0_final = [x0; lamb0];
[time_step1, z_step1] = ode45(@(t,z) UBT_eqn(t,z,alpha), [t0 t1], z0_final);

x_traj_step1 = z_step1(:,1);
y_traj_step1 = z_step1(:,2);
z_traj_step1 = z_step1(:,3);
lamb_v_step1 = z_step1(:,10:12);
gamma_x_step1 = -norm^2*lamb_v_step1(:,1);
gamma_y_step1 = -norm^2*lamb_v_step1(:,2);
gamma_z_step1 = -norm^2*lamb_v_step1(:,3);

time_step1 = time_step1/norm; % [sec]

figure()
plot(time_step1, [gamma_x_step1, gamma_y_step1, gamma_z_step1]);
legend('\Gamma_x', '\Gamma_y','\Gamma_z');
ylim([-3*1e-3 6*1e-3]);
ylabel('\Gamma [m/sec^2]'); xlabel('Time [sec]');
title('Three-dimensional first-stage optimal thrust acceleration')
grid on

%% Step 2.
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
plot(time_step2, [gamma_x_step2, gamma_y_step2, gamma_z_step2]);
ylim([-9*1e-4 9*1e-4]);
ylabel('\Gamma [m/sec^2]'); xlabel('Time [sec]');
legend('\Gamma_x', '\Gamma_y','\Gamma_z');
title('Three-dimensional second-stage optimal thrust acceleration')
grid on


%% plot trajectories
figure()
plot3(x_traj_step1/1000, y_traj_step1/1000, z_traj_step1/1000);
grid on
hold on
plot3(x_traj_step2/1000, y_traj_step2/1000, z_traj_step2/1000);

axis square
daspect([1 1 1]);
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
title('Unbounded Thrust 3D Optimal Trajectory')
