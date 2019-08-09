function J = duffing_dx(t,x,y,p)
%CATENARY   'coll'-compatible encoding of catenary vector field

x1 = x(1,:);

zeta  = p(1,:);
mu    = p(2,:);

J = zeros(2, 2, numel(t));
J(1,2,:) = 1;
J(2,1,:) = -1-3*mu.*x1.^2;
J(2,2,:) = -2*zeta;
% y(1,:) = x2;
% y(2,:) = -2*zeta.*x2-x1-mu.*x1.^3+2*a.*y1+2*b.*y2+f.*cos(Omega.*t+phi);

end