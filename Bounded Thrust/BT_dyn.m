function xdot = BT_dyn(t, x, Gamma_max, lamb0)

philamb11 = [1 0 0; 6*(sin(t)-t) 4-3*cos(t) 0; 0 0 cos(t)];
philamb12 = [0 0 0; 6*(1-cos(t)) -3*sin(t) 0; 0 0 sin(t)];
philamb21 = [3*t-4*sin(t) -2*(1-cos(t)) 0; 2*(1-cos(t)) -sin(t) 0; 0 0 -sin(t)];
philamb22 = [4*cos(t)-3 2*sin(t) 0; -2*sin(t) cos(t) 0; 0 0 cos(t)];

philamb = [philamb11, philamb12;
           philamb21, philamb22];

lambda = philamb*lamb0;

lamb_v = lambda(4:6);

if(vecnorm(lamb_v) <= Gamma_max)
    h = 1;
else
    h = Gamma_max/vecnorm(lamb_v);
end
A1 = [0 0 0; 0 3 0; 0 0 -1];
A2 = [0 2 0; -2 0 0 ; 0 0 0];

A = [zeros(3,3), eye(3);
    A1, A2];

xdot = A * x - [zeros(3, 1); h * lamb_v];
end