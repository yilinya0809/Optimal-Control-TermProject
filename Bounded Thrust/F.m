function F = F(lamb0, t0, t1, x0, x1, Gamma_max)

[~, x] = ode45(@(t,x) BT_dyn(t, x, Gamma_max, lamb0), [t0 t1], x0);
x1_ode = x(end,1:6);

F = x1_ode.' - x1;

end