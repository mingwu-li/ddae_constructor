function y = catenary(t,x,y,p)
%CATENARY   'coll'-compatible encoding of catenary vector field

x1 = x(1,:);
x2 = x(2,:);

y(1,:) = y;
y(2,:) = (1+x2.^2)./x1;

end
