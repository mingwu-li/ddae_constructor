function out = duffing(t,x,y,p)
%CATENARY   'coll'-compatible encoding of catenary vector field

x1 = x(1,:);
x2 = x(2,:);
y1 = y(1,:);
y2 = y(2,:);

zeta  = p(1,:);
mu    = p(2,:);
a     = p(3,:);
b     = p(4,:);
alpha = p(5,:);
f     = p(6,:);
Omega = p(7,:);
phi   = p(8,:);

out(1,:) = x2;
out(2,:) = -2*zeta.*x2-x1-mu.*x1.^3+2*a.*y1+2*b.*y2+f.*cos(Omega.*t+phi);

end
