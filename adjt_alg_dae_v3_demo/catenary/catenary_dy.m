function J = catenary_dy(t,x,y,p)
%CATENARY   'coll'-compatible encoding of catenary vector field

J = zeros(2,1,numel(t));
J(1,1,:) = 1;

end