function data = ddaecoll_adjt_init_data(prob, src_data) %#ok<INUSL>
%COLL_ADJT_INIT_DATA   Initialize data structure for a 'coll' adjoint problem.
%
% See also ODE_INIT_DATA, COCO_SAVE_DATA, COCO_FUNC_DATA.

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: ode_init_data.m 2902 2015-10-09 18:06:32Z hdankowicz $

data.ddaecoll_seg = src_data;

seg  = src_data.ddaecoll;
NCOL = seg.NCOL;
NTST = seg.NTST;
dim  = src_data.dim;
ydim = src_data.ydim;
pdim = src_data.pdim;


bpdim  = dim*(NCOL+1);      % Number of basepoint values per interval
xbpdim = dim*(NCOL+1)*NTST; % Number of basepoint values
ybpdim = ydim*NCOL*NTST;    % Number of basepoint values for y
cndim  = dim*NCOL;          % Number of collocation node values per interval
xcnnum = NCOL*NTST;         % Number of collocation nodes
xcndim = dim*NCOL*NTST;     % Number of collocation conditions
ycndim = ybpdim;            % Number of collocation conditions for y 
cntnum = NTST-1;            % Number of internal boundaries
cntdim = dim*(NTST-1);      % Number of continuity conditions


addim  = xcndim+cntdim+2*dim+ycndim+2+pdim;

opt.adim    = [xbpdim, addim];
opt.xcn_idx = (1:xcndim)';
opt.x0_idx  = xcndim + cntdim + (1:dim)';
opt.x1_idx  = opt.x0_idx(end) + (1:dim)';
opt.ycn_idx = opt.x1_idx(end) + (1:ycndim)';
if ycndim~=0
opt.T0_idx  = opt.ycn_idx(end) + 1;
else
opt.T0_idx  = opt.x1_idx(end) + 1;
end
opt.T_idx   = opt.T0_idx + 1;
opt.p_idx   = opt.T_idx + (1:pdim)';

opt.fbp_idx = 1:xbpdim;

opt.id0 = [eye(dim); zeros(xbpdim-dim,dim)];
opt.id1 = [zeros(xbpdim-dim,dim); eye(dim)];

opt.xtr    = (1:addim)';
opt.xtr(1:xcndim+cntdim) = 0;
opt.xtrend = opt.xtr(end-2*dim-pdim-1:end);
opt.ftr    = (1:xbpdim)';
opt.ftr(dim+1:end-dim) = 0;
opt.ftrend = opt.ftr(end-dim+1:end);

opt.dJrows = xcndim;
opt.dJcols = addim*(xbpdim+2+pdim);

% opt.dpdTtcn = repmat(permute(mesh.tcn, [2 3 1]), [dim, pdim]);
tcn         = src_data.Md;
opt.dTtcn   = repmat(tcn', [dim 1]);
% opt.dxdTtcn = repmat(permute(mesh.tcn, [2 3 1]), [dim, dim]);

% for Jacobian of adjoint with respect to delta_x
% rows = repmat(1:dim, [1 dim]);
% rows = repmat(rows(:), [1 dim]) + xcndim*repmat(0:dim-1, [dim^2 1]);
% rows = repmat(rows(:), [1 xcnnum]) + dim*repmat(0:xcnnum-1, [dim^3 1]);
% opt.dxdxrows1 = rows;
% cols = kron(1:xcndim, ones(1,dim));
% opt.dxdxcols1 = repmat(reshape(cols, [dim^2 xcnnum]), [dim 1]);
% 
% dxdxrows = repmat(reshape(1:xcndim, [cndim NTST]), [dim*bpdim, 1]);
% cols = 1 + dim*(0:NCOL-1);
% cols = repmat(cols(:), [1 dim]) + repmat(0:dim-1, [NCOL 1]);
% cols = repmat(cols(:), [1 bpdim]) + addim*repmat(0:bpdim-1, [cndim 1]);
% cols = repmat(cols(:), [1 NTST]) + ...
%   (cndim+addim*bpdim)*repmat(0:NTST-1, [cndim*bpdim 1]);
% dxdxcols = kron(cols, ones(dim,1));
% 
% idx = 1:cndim;
% idx = repmat(idx(:), [1 dim]) + xcndim*repmat(0:dim-1, [cndim 1]);
% idx = repmat(idx(:), [1 bpdim]) + ...
%   dim*xcndim*repmat(0:bpdim-1, [dim*cndim 1]);
% idx = repmat(idx(:), [1 NTST]) + ...
%   dim*(NCOL+cndim*xbpdim)*repmat(0:NTST-1, [dim*cndim*bpdim 1]);
% opt.dxdxidx = idx(:);
% 
% dxdT0rows = repmat(reshape(1:xcndim, [dim xcnnum]), [dim 1]);
% dxdT0cols = xbpdim*addim + repmat(1:xcndim, [dim 1]);
% dxdTrows  = repmat(reshape(1:xcndim, [dim xcnnum]), [dim 1]);
% dxdTcols  = (xbpdim+1)*addim + repmat(1:xcndim, [dim 1]);
% 
% dxdprows = repmat(reshape(1:xcndim, [dim xcnnum]), [dim*pdim 1]);
% cols = (xbpdim+2)*addim + repmat(1:dim, [dim 1]);
% cols = repmat(cols(:), [1 pdim]) + addim*repmat(0:pdim-1, [dim^2 1]);
% cols = repmat(cols(:), [1 xcnnum]) + ...
%   dim*repmat(0:xcnnum-1, [dim^2*pdim 1]);
% dxdpcols = cols;
% 
% opt.dxrows = [dxdxrows(:); dxdT0rows(:); dxdTrows(:); dxdprows(:)];
% opt.dxcols = [dxdxcols(:); dxdT0cols(:); dxdTcols(:); dxdpcols(:)];
% 
% % for Jacobian of adjoint with respect to delta_T0 & delta_T
% dT0dxrows = repmat(reshape(1:xcndim, [cndim NTST]), [bpdim 1]);
% dT0dxcols = addim-pdim-1+addim*repmat(0:xbpdim-1, [cndim 1]);
% 
% idx = 1:cndim;
% idx = repmat(idx(:), [1 bpdim]) + xcndim*repmat(0:bpdim-1, [cndim 1]);
% idx = repmat(idx(:), [1 NTST]) + ...
%   (cndim+bpdim*xcndim)*repmat(0:NTST-1, [cndim*bpdim 1]);
% opt.dT0dxidx = idx(:);
% opt.dTdxidx  = opt.dT0dxidx;
% 
% dT0dT0rows = 1:xcndim;
% dT0dTrows  = 1:xcndim;
% dT0dT0cols = (xcndim+cntdim+2*dim+1+xbpdim*addim)*ones(xcndim,1);
% dT0dTcols  = (xcndim+cntdim+2*dim+1+(xbpdim+1)*addim)*ones(xcndim,1);
% 
% dT0dprows = repmat(reshape(1:xcndim, [dim xcnnum]), [pdim 1]);
% cols = xcndim+cntdim+2*dim+(xbpdim+2)*addim + ones(dim,1);
% cols = repmat(cols(:), [1 pdim]) + addim*repmat(0:pdim-1, [dim 1]);
% dT0dpcols = repmat(cols(:), [1 xcnnum]);
% 
% dT0rows = [dT0dxrows(:); dT0dT0rows(:); dT0dTrows(:); dT0dprows(:)];
% dT0cols = [dT0dxcols(:); dT0dT0cols(:); dT0dTcols(:); dT0dpcols(:)];
% 
% % for Jacobian of adjoint with respect to delta_p
% rows = reshape(1:xcndim*pdim, [dim xcnnum pdim]);
% rows = permute(repmat(rows, [dim 1 1]), [1 3 2]);
% opt.dpdxrows1 = rows(:,:);
% cols = kron(1:xcndim, ones(1,dim));
% opt.dpdxcols1 = repmat(reshape(cols, [dim^2 xcnnum]), [pdim 1]);
% 
% dpdxrows2 = repmat(reshape(1:xcndim, [cndim NTST]), [pdim*bpdim 1]);
% cols = repmat(1:pdim, [cndim, 1]);
% dpdxcols2 = repmat(cols(:), [1, xbpdim]) + ...
%   xcndim+cntdim+2*dim+2+addim*repmat(0:xbpdim-1, [cndim*pdim 1]);
% 
% idx = 1:cndim;
% idx = repmat(idx(:), [1 pdim*bpdim]) + ...
%   xcndim*repmat(0:pdim*bpdim-1, [cndim 1]);
% idx = repmat(idx(:), [1 NTST]) + ...
%   (cndim+pdim*bpdim*xcndim)*repmat(0:NTST-1, [cndim*pdim*bpdim 1]);
% opt.dpdxidx = idx(:);
% 
% dpdT0rows = maps.fdprows;
% dpdT0cols = maps.fdpcols + xcndim+cntdim+2*dim+2+xbpdim*addim;
% dpdTrows  = maps.fdprows;
% dpdTcols  = maps.fdpcols + xcndim+cntdim+2*dim+2+(xbpdim+1)*addim;
% 
% dpdprows = repmat(reshape(1:xcndim, [dim, xcnnum]), [pdim^2, 1]);
% cols = repmat(1:pdim, [dim 1]);
% cols = repmat(cols(:), [1 pdim]) + xcndim+cntdim+2*dim+2 + ...
%   (xbpdim+2)*addim+addim*repmat(0:pdim-1, [dim*pdim 1]);
% cols = repmat(cols(:), [1 xcnnum]);
% dpdpcols = cols;
% 
% dprows = [dpdxrows2(:); dpdT0rows(:); dpdTrows(:); dpdprows(:)];
% dpcols = [dpdxcols2(:); dpdT0cols(:); dpdTcols(:); dpdpcols(:)];
% 
% opt.dT0Tprows = [dT0rows; dT0rows; dprows];
% opt.dT0Tpcols = [dT0cols; 1+dT0cols; dpcols];

% opt.fid   = coco_get_id(data.oid, 'coll');
    
data.ddaecoll_opt  = opt;

end
