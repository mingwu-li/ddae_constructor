function [data y] = bc(prob, data, u)

x0 = u(1:2);
x1 = u(3:4);
T0 = u(5);
T  = u(6);
Omega = u(7);

y = [x0-x1; x0(2); T0; T-2*pi/Omega];

end