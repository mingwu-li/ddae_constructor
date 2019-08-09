function J = duffing_dy(t,x,y,p)
%CATENARY   'coll'-compatible encoding of catenary vector field

a     = p(3,:);
b     = p(4,:);

J = zeros(2, 2, numel(t));
J(2,1,:) = 2*a;
J(2,2,:) = 2*b;
% out(1,:) = x2;
% out(2,:) = -2*zeta.*x2-x1-mu.*x1.^3+2*a.*y1+2*b.*y2+f.*cos(Omega.*t+phi);

end