function lamb0 = UBT(alpha, t0, tf, x0, xf)
% find initial costate value for each step

A1 = [0 0 0; 0 3 0; 0 0 -1];
A2 = [0 2 0; -2 0 0 ; 0 0 0];
D = [0 0 0; 0 1 0; 0 0 1];
A = [zeros(3,3), eye(3), zeros(3,6);
    A1, A2, zeros(3,3), -eye(3);
    -alpha*D, zeros(3,6), -A1.';
    zeros(3,6), -eye(3), -A2.'];

Phi = expm(A*(tf-t0));
Phi11 = Phi(1:6,1:6); Phi12 = Phi(1:6,7:12);

lamb0 = Phi12\(xf-Phi11*x0);


end