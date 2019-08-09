function f = model1(t,x,y,p)
%CATENARY   'coll'-compatible encoding of catenary vector field

x1 = x(1,:);
y1 = y(1,:);
y2 = y(2,:);
y3 = y(3,:);

alpha = p(1,:);
% tau1  = p(2,:);
% tau2  = p(3,:);

f = x1.^2-0.5*y1.*y2+alpha.*sin(y1)-8*x1.*y3+1;

% x1(t)' =
% x1(t)^2-0.5*x1(t-tau1)*x2(t)+alpha*sin(x1(t-tau))-8*x1(t)*x2(t-tau2)+1
% with algebraic equations
% -2*x1(t)^3+x2(t)*(1+x1(t-tau1))
% In our framework, we have x=x1, y1=x1(t-tau1), y2=x2, y3=y2(t-tau2)


end
