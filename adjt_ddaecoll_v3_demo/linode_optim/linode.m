function f = linode(t, x, y, p)

x1 = x(1,:);
x2 = x(2,:);
k  = p(1,:);
th = p(2,:);

f(1,:) = y;
f(2,:) = -y-k.*x1+cos(t+th);

end
