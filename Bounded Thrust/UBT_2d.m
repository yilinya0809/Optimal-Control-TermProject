function lamb0 = UBT_2d(alpha, t0, tf, x0, xf)

A1 = [0 0; 0 3];
A2 = [0 2; -2 0];
D = [0 0; 0 1];

A = [zeros(2,2), eye(2), zeros(2,4);
    A1, A2, zeros(2,2), -eye(2);
    -alpha*D, zeros(2,4), -A1.';
    zeros(2,4), -eye(2), -A2.'];

Phi = expm(A*(tf-t0));
Phi11 = Phi(1:4, 1:4); Phi12 = Phi(1:4, 5:8);

lamb0 = Phi12\(xf-Phi11*x0);


end