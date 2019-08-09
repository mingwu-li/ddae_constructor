function J = linode_dy(t, x, y, p) %#ok<INUSL>

k = p(1,:);

J = zeros(2,1,numel(t));
J(1,1,:) = 1;
J(2,1,:) = -1;

end
