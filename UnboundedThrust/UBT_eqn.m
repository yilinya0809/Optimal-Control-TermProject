function UBT_eqn = UBT_eqn(t, z, alpha)

A1 = [0 0 0; 0 3 0; 0 0 -1];
A2 = [0 2 0; -2 0 0 ; 0 0 0];

D = [0 0 0; 0 1 0; 0 0 1];
A = [zeros(3,3), eye(3), zeros(3,6);
    A1, A2, zeros(3,3), -eye(3);
    -alpha*D, zeros(3,6), -A1.';
    zeros(3,6), -eye(3), -A2.'];

UBT_eqn = A*z;

end