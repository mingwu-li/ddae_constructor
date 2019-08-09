function J = linode_dx(t, x, y, p) %#ok<INUSL>

k = p(1,:);

J = zeros(2,2,numel(t));
J(2,1,:) = -k;

end
