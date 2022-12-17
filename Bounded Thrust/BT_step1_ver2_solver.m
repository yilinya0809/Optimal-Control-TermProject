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

options = optimoptions('fsolve', 'Display', 'iter');
[lamb, fval, exitflag, output] = fsolve(@(lamb0)F(lamb0, t0, t1, x0, x1, Gamma_max), ...
    lamb0, options);

lamb0 = lamb;

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
