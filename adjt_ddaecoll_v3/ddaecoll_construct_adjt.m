function [prob, data] = ddaecoll_construct_adjt(prob, tbid, data, sol)
%COLL_CONSTRUCT_ADJT   Add COLL adjoint problem.

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: coll_construct_opt.m 2872 2015-08-06 20:12:06Z hdankowicz $

% [data, sol] = init_data(data, sol);
opt  = data.ddaecoll_opt;
seg  = data.ddaecoll_seg;
if ~isempty(sol.l0)
  sol.l0  = interp1(sol.tbp', sol.l0', seg.tbpd', 'pchip')';
  sol.l0  = sol.l0(:);
  if ~isempty(sol.tl0)
    sol.tl0 = interp1(sol.tbp', sol.tl0', seg.tbpd', 'pchip')';
    sol.tl0 = sol.tl0(:);
  end
end
prob = coco_add_adjt(prob, tbid, @adj, @adj_DU, data, 'l0', sol.l0, ...
  'tl0', sol.tl0, 'adim', opt.adim);
if ~isempty(seg.pnames)
  pfid   = coco_get_id(tbid, 'pars');
  dnames = coco_get_id('d', seg.pnames);
  axidx  = coco_get_adjt_data(prob, tbid, 'axidx');
  prob   = coco_add_adjt(prob, pfid, dnames, 'aidx', axidx(opt.p_idx), ...
    'l0', sol.pars_l0, 'tl0', sol.pars_tl0);
end

end

function [data, sol] = init_data(data, sol)
%INIT_DATA   Initialize data for COLL adjoint problem.

seg  = data.ddaecoll_seg;
int  = seg.int;
maps = seg.maps;
mesh = seg.mesh;

NCOL = int.NCOL;
NTST = maps.NTST;
xdim = int.dim;
pdim = maps.pdim;

cndim  = NCOL*xdim;
bpdim  = xdim*(NCOL+1);
xbpdim = NTST*(NCOL+1)*xdim;
xcnnum = NTST*NCOL;
xcndim = NTST*NCOL*xdim;
cntdim = (NTST-1)*xdim;
addim  = xcndim+cntdim+2*xdim+2+pdim;

opt.adim    = [xbpdim, addim];
opt.xcn_idx = (1:xcndim)';
opt.x0_idx  = xcndim + cntdim + (1:xdim)';
opt.x1_idx  = opt.x0_idx(end) + (1:xdim)';
opt.T0_idx  = opt.x1_idx(end) + 1;
opt.T_idx   = opt.T0_idx + 1;
opt.p_idx   = opt.T_idx + (1:pdim)';

opt.fbp_idx = 1:xbpdim;

opt.id0 = [eye(xdim); zeros(xbpdim-xdim,xdim)];
opt.id1 = [zeros(xbpdim-xdim,xdim); eye(xdim)];

opt.xtr    = (1:addim)';
opt.xtr(1:xcndim+cntdim) = 0;
opt.xtrend = opt.xtr(end-2*xdim-pdim-1:end);
opt.ftr    = (1:xbpdim)';
opt.ftr(xdim+1:end-xdim) = 0;
opt.ftrend = opt.ftr(end-xdim+1:end);

opt.dJrows = xcndim;
opt.dJcols = addim*(xbpdim+2+pdim);

opt.dpdTtcn = repmat(permute(mesh.tcn, [2 3 1]), [xdim, pdim]);
opt.dTtcn   = repmat(mesh.tcn', [xdim 1]);
opt.dxdTtcn = repmat(permute(mesh.tcn, [2 3 1]), [xdim, xdim]);

% for Jacobian of adjoint with respect to delta_x
rows = repmat(1:xdim, [1 xdim]);
rows = repmat(rows(:), [1 xdim]) + xcndim*repmat(0:xdim-1, [xdim^2 1]);
rows = repmat(rows(:), [1 xcnnum]) + xdim*repmat(0:xcnnum-1, [xdim^3 1]);
opt.dxdxrows1 = rows;
cols = kron(1:xcndim, ones(1,xdim));
opt.dxdxcols1 = repmat(reshape(cols, [xdim^2 xcnnum]), [xdim 1]);

dxdxrows = repmat(reshape(1:xcndim, [cndim NTST]), [xdim*bpdim, 1]);
cols = 1 + xdim*(0:NCOL-1);
cols = repmat(cols(:), [1 xdim]) + repmat(0:xdim-1, [NCOL 1]);
cols = repmat(cols(:), [1 bpdim]) + addim*repmat(0:bpdim-1, [cndim 1]);
cols = repmat(cols(:), [1 NTST]) + ...
  (cndim+addim*bpdim)*repmat(0:NTST-1, [cndim*bpdim 1]);
dxdxcols = kron(cols, ones(xdim,1));

idx = 1:cndim;
idx = repmat(idx(:), [1 xdim]) + xcndim*repmat(0:xdim-1, [cndim 1]);
idx = repmat(idx(:), [1 bpdim]) + ...
  xdim*xcndim*repmat(0:bpdim-1, [xdim*cndim 1]);
idx = repmat(idx(:), [1 NTST]) + ...
  xdim*(NCOL+cndim*xbpdim)*repmat(0:NTST-1, [xdim*cndim*bpdim 1]);
opt.dxdxidx = idx(:);

dxdT0rows = repmat(reshape(1:xcndim, [xdim xcnnum]), [xdim 1]);
dxdT0cols = xbpdim*addim + repmat(1:xcndim, [xdim 1]);
dxdTrows  = repmat(reshape(1:xcndim, [xdim xcnnum]), [xdim 1]);
dxdTcols  = (xbpdim+1)*addim + repmat(1:xcndim, [xdim 1]);

dxdprows = repmat(reshape(1:xcndim, [xdim xcnnum]), [xdim*pdim 1]);
cols = (xbpdim+2)*addim + repmat(1:xdim, [xdim 1]);
cols = repmat(cols(:), [1 pdim]) + addim*repmat(0:pdim-1, [xdim^2 1]);
cols = repmat(cols(:), [1 xcnnum]) + ...
  xdim*repmat(0:xcnnum-1, [xdim^2*pdim 1]);
dxdpcols = cols;

opt.dxrows = [dxdxrows(:); dxdT0rows(:); dxdTrows(:); dxdprows(:)];
opt.dxcols = [dxdxcols(:); dxdT0cols(:); dxdTcols(:); dxdpcols(:)];

% for Jacobian of adjoint with respect to delta_T0 & delta_T
dT0dxrows = repmat(reshape(1:xcndim, [cndim NTST]), [bpdim 1]);
dT0dxcols = addim-pdim-1+addim*repmat(0:xbpdim-1, [cndim 1]);

idx = 1:cndim;
idx = repmat(idx(:), [1 bpdim]) + xcndim*repmat(0:bpdim-1, [cndim 1]);
idx = repmat(idx(:), [1 NTST]) + ...
  (cndim+bpdim*xcndim)*repmat(0:NTST-1, [cndim*bpdim 1]);
opt.dT0dxidx = idx(:);
opt.dTdxidx  = opt.dT0dxidx;

dT0dT0rows = 1:xcndim;
dT0dTrows  = 1:xcndim;
dT0dT0cols = (xcndim+cntdim+2*xdim+1+xbpdim*addim)*ones(xcndim,1);
dT0dTcols  = (xcndim+cntdim+2*xdim+1+(xbpdim+1)*addim)*ones(xcndim,1);

dT0dprows = repmat(reshape(1:xcndim, [xdim xcnnum]), [pdim 1]);
cols = xcndim+cntdim+2*xdim+(xbpdim+2)*addim + ones(xdim,1);
cols = repmat(cols(:), [1 pdim]) + addim*repmat(0:pdim-1, [xdim 1]);
dT0dpcols = repmat(cols(:), [1 xcnnum]);

dT0rows = [dT0dxrows(:); dT0dT0rows(:); dT0dTrows(:); dT0dprows(:)];
dT0cols = [dT0dxcols(:); dT0dT0cols(:); dT0dTcols(:); dT0dpcols(:)];

% for Jacobian of adjoint with respect to delta_p
rows = reshape(1:xcndim*pdim, [xdim xcnnum pdim]);
rows = permute(repmat(rows, [xdim 1 1]), [1 3 2]);
opt.dpdxrows1 = rows(:,:);
cols = kron(1:xcndim, ones(1,xdim));
opt.dpdxcols1 = repmat(reshape(cols, [xdim^2 xcnnum]), [pdim 1]);

dpdxrows2 = repmat(reshape(1:xcndim, [cndim NTST]), [pdim*bpdim 1]);
cols = repmat(1:pdim, [cndim, 1]);
dpdxcols2 = repmat(cols(:), [1, xbpdim]) + ...
  xcndim+cntdim+2*xdim+2+addim*repmat(0:xbpdim-1, [cndim*pdim 1]);

idx = 1:cndim;
idx = repmat(idx(:), [1 pdim*bpdim]) + ...
  xcndim*repmat(0:pdim*bpdim-1, [cndim 1]);
idx = repmat(idx(:), [1 NTST]) + ...
  (cndim+pdim*bpdim*xcndim)*repmat(0:NTST-1, [cndim*pdim*bpdim 1]);
opt.dpdxidx = idx(:);

dpdT0rows = maps.fdprows;
dpdT0cols = maps.fdpcols + xcndim+cntdim+2*xdim+2+xbpdim*addim;
dpdTrows  = maps.fdprows;
dpdTcols  = maps.fdpcols + xcndim+cntdim+2*xdim+2+(xbpdim+1)*addim;

dpdprows = repmat(reshape(1:xcndim, [xdim, xcnnum]), [pdim^2, 1]);
cols = repmat(1:pdim, [xdim 1]);
cols = repmat(cols(:), [1 pdim]) + xcndim+cntdim+2*xdim+2 + ...
  (xbpdim+2)*addim+addim*repmat(0:pdim-1, [xdim*pdim 1]);
cols = repmat(cols(:), [1 xcnnum]);
dpdpcols = cols;

dprows = [dpdxrows2(:); dpdT0rows(:); dpdTrows(:); dpdprows(:)];
dpcols = [dpdxcols2(:); dpdT0cols(:); dpdTcols(:); dpdpcols(:)];

opt.dT0Tprows = [dT0rows; dT0rows; dprows];
opt.dT0Tpcols = [dT0cols; 1+dT0cols; dpcols];

opt.fid   = coco_get_id(data.oid, 'coll');
    
data.coll_opt  = opt;

if nargout==2 && ~isempty(sol.l0)
  sol.l0  = interp1(sol.tbp', sol.l0', mesh.tbp, 'pchip')';
  sol.l0  = sol.l0(:);
  if ~isempty(sol.tl0)
    sol.tl0 = interp1(sol.tbp', sol.tl0', mesh.tbp, 'pchip')';
    sol.tl0 = sol.tl0(:);
  end
end

end

function [data, J] = adj(prob, data, u) %#ok<INUSL>

seg  = data.ddaecoll_seg;
opt  = data.ddaecoll_opt;
NTST = seg.ddaecoll.NTST;

x  = u(seg.xbp_idx);
y  = u(seg.ybp_idx);
T0 = u(seg.T0_idx);
T  = u(seg.T_idx);
p  = u(seg.p_idx);

xcn = reshape(seg.Wdd*x, seg.x_shp); % Values at collocation nodes
pcn = repmat(p, seg.p_rep);
ycn = reshape(seg.Wda*y, seg.y_shp); 
tcn = seg.Md*T+T0;

xcna = reshape(seg.Wad*x, seg.x_shp); % Values at collocation nodes
ycna = reshape(seg.Waa*y, seg.y_shp); 
tcna = seg.Ma*T+T0;

fcn   = seg.fhan(tcn', xcn, ycn, pcn);
fdtcn = seg.dfdthan(tcn', xcn, ycn, pcn);
fdxcn = seg.dfdxhan(tcn', xcn, ycn, pcn);
fdycn = seg.dfdyhan(tcna', xcna, ycna, pcn); % evaluated at collocation set a 
fdpcn = seg.dfdphan(tcn', xcn, ycn, pcn);
 
% adjoint with respect to delta_x
dxode = sparse(seg.dxrows, seg.dxcols, fdxcn(:));
J = -2*NTST*seg.Wddp'-T*seg.Wdd'*dxode;

% adjoint with respect to delta_x^(j+1)(-1), delta_x^(1)(-1), and delta_x^(N)(1)
J = [ J, seg.Q', -opt.id0, opt.id1 ];

% adjoint with respect to delta_y
dyode = sparse(seg.dyrows, seg.dycols, fdycn(:));
if isempty(dyode)
    J = [ J, [] ];
else
    J = [ J, -T*seg.Wad'*dyode ];
end


% adjoint with respect to T0 and T
dT0ode = T*fdtcn;
dTode  = fcn + T*fdtcn.*opt.dTtcn;
J = [ J, -(0.5/NTST)*seg.Wdd'*seg.wts2*[dT0ode(:) dTode(:)] ];

% adjoint with respect to p
dpode = fdpcn;
dpode = sparse(seg.dprows, seg.dpcols, dpode(:));
if isempty(dpode)
    J = [ J, [] ];
else
    J = [ J, -(0.5*T/NTST)*seg.Wdd'*seg.wts2*dpode ];
end

end

function [data, dJ] = adj_DU(prob, data, u) %#ok<INUSL>

[data, dJ] = coco_ezDFDX('f(o,d,x)', prob, data, @adj, u);

% [data, DadF] = coco_ezDFDX('f(o,d,x)',  prob, data, @ddaecoll_F, u);

% pr   = data.pr;
% opt  = pr.coll_opt;
% seg  = pr.coll_seg;
% maps = seg.maps;
% mesh = seg.mesh;
% 
% x  = u(maps.xbp_idx);
% T0 = u(maps.T0_idx);
% T  = u(maps.T_idx);
% p  = u(maps.p_idx);
% 
% xcn = reshape(maps.W*x, maps.x_shp);
% pcn = repmat(p, maps.p_rep);
% tcn = T0+T*mesh.tcn';
% 
% fdxcn   =   pr.ode_DFDX(pr, tcn, xcn, pcn);
% fdpcn   =   pr.ode_DFDP(pr, tcn, xcn, pcn);
% fdtcn   =   pr.ode_DFDT(pr, tcn, xcn, pcn);
% fdxdxcn = pr.ode_DFDXDX(pr, tcn, xcn, pcn);
% fdxdpcn = pr.ode_DFDXDP(pr, tcn, xcn, pcn);
% fdpdpcn = pr.ode_DFDPDP(pr, tcn, xcn, pcn);
% fdtdtcn = pr.ode_DFDTDT(pr, tcn, xcn, pcn);
% fdtdpcn = pr.ode_DFDTDP(pr, tcn, xcn, pcn);
% fdtdxcn = pr.ode_DFDTDX(pr, tcn, xcn, pcn);
% 
% dJrows = opt.dJrows;
% dJcols = opt.dJcols;
% 
% % Jacobians of adjoint with respect to delta_x
% dxdxode  = mesh.fdxdxka.*(T*fdxdxcn);
% dxdxode  = sparse(opt.dxdxrows1, opt.dxdxcols1, dxdxode(:))*maps.W;
% dxdT0ode = mesh.fdxka.*(T*fdtdxcn);
% dxdTode  = mesh.fdxka.*(fdxcn+T*fdtdxcn.*opt.dxdTtcn);
% dxdpode  = mesh.fdxdpka.*(T*fdxdpcn);
% 
% dxvals = [dxdxode(opt.dxdxidx); dxdT0ode(:); dxdTode(:); dxdpode(:)];
% 
% % Jacobians of adjoint with respect to T0, T, and p
% dT0dxode  = mesh.fdxka.*(T*fdtdxcn);
% dT0dxode  = sparse(maps.fdxrows, maps.fdxcols, dT0dxode(:))*maps.W;
% dT0dT0ode = mesh.fka.*(T*fdtdtcn);
% dT0dTode  = mesh.fka.*(fdtcn+T*fdtdtcn.*opt.dTtcn);
% dT0dpode  = mesh.fdpka.*(T*fdtdpcn);
% dTdxode   = mesh.fdxka.*(fdxcn+T*fdtdxcn.*opt.dxdTtcn);
% dTdxode   = sparse(maps.fdxrows, maps.fdxcols, dTdxode(:))*maps.W;
% dTdT0ode  = mesh.fka.*(fdtcn+T*fdtdtcn.*opt.dTtcn);
% dTdTode   = mesh.fka.*(2*fdtcn+T*fdtdtcn.*opt.dTtcn).*opt.dTtcn;
% dTdpode   = mesh.fdpka.*(fdpcn+T*fdtdpcn.*opt.dpdTtcn);
% dpdxode  = mesh.fdxdpka.*(T*fdxdpcn);
% dpdxode  = sparse(opt.dpdxrows1, opt.dpdxcols1, dpdxode(:))*maps.W;
% dpdT0ode = mesh.fdpka.*(T*fdtdpcn);
% dpdTode  = mesh.fdpka.*(fdpcn+T*fdtdpcn.*opt.dpdTtcn);
% dpdpode  = mesh.fdpdpka.*(T*fdpdpcn);
% 
% dT0vals = [dT0dxode(opt.dT0dxidx); dT0dT0ode(:); dT0dTode(:); dT0dpode(:)];
% dTvals  = [dTdxode(opt.dTdxidx); dTdT0ode(:); dTdTode(:); dTdpode(:)];
% dpvals  = [dpdxode(opt.dpdxidx); dpdT0ode(:); dpdTode(:); dpdpode(:)];
% 
% dT0Tpvals = [dT0vals; dTvals; dpvals];
% 
% dJ = sparse(opt.dT0Tprows, opt.dT0Tpcols, dT0Tpvals, dJrows, dJcols);
% dJ = mesh.wts2*dJ;
% dJ = dJ + sparse(opt.dxrows, opt.dxcols, dxvals, dJrows, dJcols);
% dJ = -(0.5/maps.NTST)*maps.W'*dJ;

end
