function zdot = BT_eqn(t,z, alpha, Gamma_max)
lamb_v = z(10:12);

if(vecnorm(lamb_v) <= Gamma_max)
    h = 1;
else
    h = Gamma_max/vecnorm(lamb_v);
end
A1 = [0 0 0; 0 3 0; 0 0 -1];
A2 = [0 2 0; -2 0 0 ; 0 0 0];
D = [0 0 0; 0 1 0; 0 0 1];

A = [zeros(3,3), eye(3), zeros(3,6);
    A1, A2, zeros(3,3), -h*eye(3);
    -alpha*D, zeros(3,6), -A1.';
    zeros(3,6), -eye(3), -A2.'];

zdot = A*z;
end