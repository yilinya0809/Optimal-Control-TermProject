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
r0 = [15000; 0]; v0 = [-10; 0]/norm; x0 = [r0; v0];
r1 = [-300; 0]; v1 = [0.2; 0]/norm; x1 = [r1; v1];
rf = [-1e-8; 0]; vf = [1e-8; 0]; xf = [rf; vf];

% Thrust constraint
Gamma_max = 5*1e-4;
gm_norm = Gamma_max/norm^2;

%% Step 1. Bounded Thrust Gradient Method
alpha = 0;

% costate initial guess from unbounded case
lamb0 = UBT_2d(alpha, t0, t1, x0, x1);

itr = 1;
itr_num = 5;
lamb0_hist = zeros(itr_num, 4);
err_hist = zeros(itr_num);

while itr <= itr_num
    if itr == 1
        [time_1, x] = ode45(@(t,x) BT_dyn_2d(t,x,gm_norm, lamb0),[t0 t1],x0);
        x1_ode = x(end,1:4);
    end
    lamb0_hist(itr,:) = lamb0';

    % define err func
    L = vecnorm(x1_ode'-x1);
    err_hist(itr) = L;

    % calculate gradient
    eps1 = 1e-6;
    dlam1 = diag(lamb0.*eps1);
    lamb0_del = zeros(4,4);
    x1_del = zeros(4,4);
    L_del = zeros(4,1);
    dLdlamb = zeros(4,1);
    for i = 1:1:4
        lamb0_del(:,i) = lamb0 + dlam1(:,i);
        [~, x_del] = ode45(@(t,x) BT_dyn_2d(t,x,gm_norm,lamb0_del(:,i)), [t0 t1], x0);
        x1_del(i,:) = x_del(end,1:4);

        L_del(i) = vecnorm(x1_del(i,:)'-x1);
        dLdlamb(i) = (L_del(i) - L) / dlam1(i,i);

    
        % 2nd order gradient
        eps2 = 1e-6;
        dlam2 = diag(lamb0.*eps2);
%         lamb0_d2d1 = zeros(4,4); lamb0_d2 = zeros(4,4);
%         x1_d2d1 = zeros(4,4); x1_d2 = zeros(4,4);
%         L_d2d1 = zeros(4,1); L_d2 = zeros(4,1);
%         dLdlamb2 = zeros(4,1);
%         d2Ldlamb2 = zeros(4,4);
        
        for j = 1:1:4
            lamb0_d2d1(:,j) = lamb0 + dlam1(:,i) + dlam2(:,j);
            lamb0_d2(:,j) = lamb0 + dlam2(:,j);
            
            [~,x_d2d1] = ode45(@(t,x) BT_dyn_2d(t,x,gm_norm,lamb0_d2d1(:,j)),[t0 t1], x0);
            [~,x_d2] = ode45(@(t,x) BT_dyn_2d(t,x,gm_norm,lamb0_d2(:,j)),[t0 t1], x0);
            x1_d2d1(j,:) = x_d2d1(end,1:4);
            x1_d2(j,:) = x_d2(end,1:4);

            L_d2d1(j) = vecnorm(x1_d2d1(j,:)'-x1);
            L_d2(j) = vecnorm(x1_d2(j,:)'-x1);
            dLdlamb2(j) = (L_d2d1(j) - L_d2(j))/dlam1(i,i);
            d2Ldlamb2(i,j) = (dLdlamb2(j) - dLdlamb(i))/dlam2(j,j);



%             [~, x_d2d1_rx] = ode45(@(t, x) BT_dyn_2d(t, x, gm_norm, lamb0_d2d1(:,1)), [t0 t1], x0);
%             [~, x_d2_rx] = ode45(@(t, x) BT_dyn_2d(t, x, gm_norm, lamb0_d2(:,1)), [t0 t1], x0);
%             [~, x_d2d1_ry] = ode45(@(t, x) BT_dyn_2d(t, x, gm_norm, lamb0_d2d1(:,2)), [t0 t1], x0);
%             [~, x_d2_ry] = ode45(@(t, x) BT_dyn_2d(t, x, gm_norm, lamb0_d2(:,2)), [t0 t1], x0);
%             [~, x_d2d1_vx] = ode45(@(t, x) BT_dyn_2d(t, x, gm_norm, lamb0_d2d1(:,3)), [t0 t1], x0);
%             [~, x_d2_vx] = ode45(@(t, x) BT_dyn_2d(t, x, gm_norm, lamb0_d2(:,3)), [t0 t1], x0);
%             [~, x_d2d1_vy] = ode45(@(t, x) BT_dyn_2d(t, x, gm_norm, lamb0_d2d1(:,4)), [t0 t1], x0);
%             [~, x_d2_vy] = ode45(@(t, x) BT_dyn_2d(t, x, gm_norm, lamb0_d2(:,4)), [t0 t1], x0);
% 
%             x1_d2d1(1,:) = x_d2d1_rx(end,1:4);
%             x1_d2d1(2,:) = x_d2d1_ry(end,1:4);
%             x1_d2d1(3,:) = x_d2d1_vx(end,1:4);
%             x1_d2d1(4,:) = x_d2d1_vy(end,1:4);
%             
%             x1_d2(1,:) = x_d2_rx(end,1:4);
%             x1_d2(2,:) = x_d2_ry(end,1:4);
%             x1_d2(3,:) = x_d2_vx(end,1:4);
%             x1_d2(4,:) = x_d2_vy(end,1:4);

%             L_d2d1(j) = vecnorm(x1_d2d1(j,:)' - x1);      
%             L_d2(j) = vecnorm(x1_d2(j,:)' - x1);
%          
%             dLdlamb2(j) = (L_d2d1(j) - L_d2(j))/dlam1(i,i);
% 
%             d2Ldlamb2(i,j) = (dLdlamb2(j) - dLdlamb(i))/dlam2(j,j);

        end


    end


    if vecnorm(dLdlamb) < 1e-50
        disp('Solution converged')
        break;
    end

    % lamb0 update
    k = 1e-6;
    lamb0_new = lamb0 - k*d2Ldlamb2\dLdlamb;
    % lamb0_new = lamb0 - k*dLdlamb;
    [~, x_new] = ode45(@(t, x) BT_dyn_2d(t, x, gm_norm, lamb0_new), [t0 t1], x0);

    x1_ode_new = x_new(end,1:4);
    err = vecnorm(x1_ode_new'-x1);


    lamb0 = lamb0_new;
    x1_ode = x1_ode_new;
    itr = itr+1;


end

figure()
plot(err_hist)
