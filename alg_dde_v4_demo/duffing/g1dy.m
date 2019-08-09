function J = g1dy(t,y,x,p)

J = zeros(2,2,numel(t));
J(1,1,:) = 1;
J(2,2,:) = 1;
% f = y-x;

end