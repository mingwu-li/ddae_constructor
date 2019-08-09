function f = gfunc(t,x,y,p)
%CATENARY   'coll'-compatible encoding of catenary vector field

x1 = x(1,:);
x2 = x(2,:);

f = x1+x2+y-1;

end