function J = catenary_dx(t,x,y,p)
%CATENARY   'coll'-compatible encoding of catenary vector field

x1 = x(1,:);
x2 = x(2,:);

J = zeros(2,2,numel(t));
J(2,1,:) = -(1+x2.^2)./x1.^2;
J(2,2,:) = 2*x2./x1;

end