function [data, y] = adj_coup(prob, data, u) %#ok<INUSL>


seg  = data.ddaecoll_seg;
NTST = seg.ddaecoll.NTST;
NCOL = seg.ddaecoll.NCOL;
g_x  = [0 1];
J_xcn = kron(eye(NTST*NCOL), g_x);
% J_xcn = eye(NTST*NCOL);
J_ycn = -eye(NTST*NCOL);
y = [seg.Wda'*J_xcn, seg.Waa'*J_ycn];
% y = rand(size(y));

end