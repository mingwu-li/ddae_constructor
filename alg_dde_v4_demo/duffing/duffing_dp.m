function J = duffing_dp(t,x,y,p)
%CATENARY   'coll'-compatible encoding of catenary vector field

x1 = x(1,:);
x2 = x(2,:);
y1 = y(1,:);
y2 = y(2,:);

% zeta  = p(1,:);
% mu    = p(2,:);
% a     = p(3,:);
% b     = p(4,:);
% alpha = p(5,:);
f     = p(6,:);
Omega = p(7,:);
phi   = p(8,:);

J = zeros(2, 8, numel(t));
J(2,1,:) = -2*x2;
J(2,2,:) = -x1.^3;
J(2,3,:) = 2*y1;
J(2,4,:) = 2*y2;
J(2,6,:) = cos(Omega.*t+phi);
J(2,7,:) = -t.*f.*sin(Omega.*t+phi);
J(2,8,:) = -f.*sin(Omega.*t+phi);


% out(1,:) = x2;
% out(2,:) = -2*zeta.*x2-x1-mu.*x1.^3+2*a.*y1+2*b.*y2+f.*cos(Omega.*t+phi);

end