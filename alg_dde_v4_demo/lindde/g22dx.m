function J = g22dx(t,y,x,p)

J = zeros(1,3,numel(t));
J(1,2,:) = -1;
% f = y-x(2,:);

end