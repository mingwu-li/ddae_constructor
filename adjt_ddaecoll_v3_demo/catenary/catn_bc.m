function [data, y] = catn_bc(prob, data, u)

T  = u(1);
x0 = u(2:3);
x1 = u(4:5);
p  = u(6);
T0 = u(7);

y = [T-1; x0(1)-1; x1(1)-p; T0];
  
end