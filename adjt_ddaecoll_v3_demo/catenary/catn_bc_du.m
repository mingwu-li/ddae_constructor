function [data, J] = catn_bc_du(prob, data, u)

T  = u(1);
x0 = u(2:3);
x1 = u(4:5);
p  = u(6);
T0 = u(7);

J = zeros(4,7);
J(1,1) = 1;
J(2,2) = 1;
J(3,4) = 1;
J(3,6) = -1;
J(4,7) = 1;

end