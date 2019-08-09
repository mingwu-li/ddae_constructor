function f = linode_v0(t, x, p)

x1 = x(1,:);
x2 = x(2,:);
k  = p(1,:);
th = p(2,:);

f(1,:) = x2;
f(2,:) = -x2-k.*x1+cos(t+th);

end