function J = duffing_dt(t,x,y,p)
%CATENARY   'coll'-compatible encoding of catenary vector field

f     = p(6,:);
Omega = p(7,:);
phi   = p(8,:);

J = zeros(2,numel(t));
J(2,:) = -f.*Omega.*sin(Omega.*t+phi);
% y(1,:) = x2;
% y(2,:) = -2*zeta.*x2-x1-mu.*x1.^3+2*a.*y1+2*b.*y2+f.*cos(Omega.*t+phi);

end