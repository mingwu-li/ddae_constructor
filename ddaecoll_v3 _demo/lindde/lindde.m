function y = lindde(t,x,y,p)
%CATENARY   'coll'-compatible encoding of catenary vector field

x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);
y1 = y(1,:);
y2 = y(2,:);

y(1,:) = y1;
y(2,:) = y1+y2;
y(3,:) = x2;

end
