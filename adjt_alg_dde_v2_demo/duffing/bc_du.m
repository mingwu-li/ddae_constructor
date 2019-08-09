function [data J] = bc_du(prob, data, u)

% x0 = u(1:2);
% x1 = u(3:4);
% T0 = u(5);
% T  = u(6);
Omega = u(7);
% 
% y = [x0-x1; x0(2); T0; T-2*pi/Omega];

J = zeros(5,7);
J(1,1) = 1;
J(1,3) = -1;
J(2,2) = 1;
J(2,4) = -1;
J(3,2) = 1;
J(4,5) = 1;
J(5,6) = 1;
J(5,7) = 2*pi/Omega^2;

end