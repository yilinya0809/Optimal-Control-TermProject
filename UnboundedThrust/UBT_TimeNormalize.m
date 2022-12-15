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
t0 = 0; t1 = norm*2*orbital_period; 
t2 = norm*1885; tf = t1+t2;

% boundary conditions 
r0 = [15000; 0; 2000]; v0 = norm*[-10; 0; -2]; x0 = [r0; v0];
r1 = [-300; 0; 0]; v1 = norm*[0.35; 0; 0];  x1 = [r1; v1];
rf = [-1e-8; 0; 0]; vf = [1e-8; 0; 0]; xf = [rf; vf];


% Finding Initial guess of lamb0 at each stage 
% Step 1

syms s t

A1 = [0 0 0; 0 3 0; 0 0 -1];
A2 = [0 2 0; -2 0 0 ; 0 0 0];
A_step1 = [zeros(3,3), eye(3), zeros(3,6);
           A1, A2, zeros(3,3), -eye(3);
           zeros(3,9), -A1.';
           zeros(3,6), -eye(3), -A2.'];


C1 = eye(6); C2 = zeros(6,6);

% Phi_step1 = ilaplace(inv(s*eye(12)-A_step1));
Phi_step1 = expm(A_step1*(t1));
Phi11 = Phi_step1(1:6,1:6); Phi12 = Phi_step1(1:6, 7:12); 
Phi21 = Phi_step1(7:12, 1:6); Phi22 = Phi_step1(7:12, 7:12); 
% Phi11 = Phi, Phi12 = Psi
% Phi21 = 0,   Phi22 = Phi_lamb

% costate initial value
lamb0 = Phi12\(x1-Phi11*x0);
z0 = double([x0; lamb0]);


% find trajectories
opts = odeset('RelTol',1e-5,'AbsTol',1e-8);
[time1,z_step1] = ode45(@(t,z) UBT_eqn(t,z,0), [t0 t1], z0, opts);

x_traj_step1 = z_step1(:,1);
y_traj_step1 = z_step1(:,2);
z_traj_step1 = z_step1(:,3);
gamma_x_traj_step1 = -norm^2*z_step1(:,10);
gamma_y_traj_step1 = -norm^2*z_step1(:,11);
gamma_z_traj_step1 = -norm^2*z_step1(:,12);

% check boundary conditions
pos_t1 = z_step1(end,1:3);
vel_t1 = z_step1(end,4:6);
disp(pos_t1);
disp(vel_t1);
disp(v1);
%%
% plot trajectories
figure()
subplot(1,2,1)
plot(x_traj_step1, y_traj_step1);
grid on


subplot(1,2,2)
plot3(x_traj_step1, y_traj_step1, z_traj_step1);
% daspect([1 1 1]);
grid on
zlim([-20000 20000]);

time1_n = time1/norm;

figure()
plot(time1_n, [gamma_x_traj_step1, gamma_y_traj_step1, gamma_z_traj_step1]);
legend('\Gamma_x','\Gamma_y','\Gamma_z');




%% Step 2

alpha = 5000;
D = [0 0 0; 0 1 0; 0 0 1];
A_step2 = [zeros(3,3), eye(3), zeros(3,6);
           A1, A2, zeros(3,3), -eye(3);
           -alpha*D, zeros(3,6), -A1.';
           zeros(3,6), -eye(3), -A2.'];

% Phi_step2 = ilaplace(inv(s*eye(12)-A_step2));
Phi_step2 = expm(A_step2*(t2));
Phi11_step2 = Phi_step2(1:6, 1:6); Phi12_step2 = Phi_step2(1:6, 7:12);

% lamb1 = vpa(subs(-Phi12_step2\(Phi11_step2*x1 - xf), t, t2));
lamb1 = Phi12_step2\(xf - Phi11_step2*x1);
z1 = double([x1; lamb1]);

[time2, z_step2] = ode45(@(t,z) UBT2_eqn(t, z, alpha), [t1 tf], z1);
x_traj_step2 = z_step2(:,1);
y_traj_step2 = z_step2(:,2);
z_traj_step2 = z_step2(:,3);
gamma_x_traj_step2 = -norm^2*z_step2(:,10);
gamma_y_traj_step2 = -norm^2*z_step2(:,11);
gamma_z_traj_step2 = -norm^2*z_step2(:,12);

time2 = time2/norm;

figure()
subplot(1,2,1)
plot(x_traj_step2, y_traj_step2);
grid on
ylim([0 10]); 

subplot(1,2,2)
plot3(x_traj_step2, y_traj_step2, z_traj_step2);
grid on
ylim([0 10]); 

figure()
plot(time2, [gamma_x_traj_step2, gamma_y_traj_step2, gamma_z_traj_step2]);
legend('\Gamma_x','\Gamma_y','\Gamma_z');


%%
figure()
plot3(x_traj_step1, y_traj_step1, z_traj_step1);
grid on
hold on
plot3(x_traj_step2, y_traj_step2, z_traj_step2);

