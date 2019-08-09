function data = alg_dde_adjt_init_data(prob, data) %#ok<INUSL>
%COLL_ADJT_INIT_DATA   Initialize data structure for a 'coll' adjoint problem.
%
% See also ODE_INIT_DATA, COCO_SAVE_DATA, COCO_FUNC_DATA.

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: ode_init_data.m 2902 2015-10-09 18:06:32Z hdankowicz $

seg  = data.ddaecoll_seg;
NTST = seg.ddaecoll.NTST; % Number of mesh intervals
NCOL = seg.ddaecoll.NCOL; % Degree of polynomial interpolants
dim  = seg.dim;           % State-space dimension
ydim = seg.ydim;
pdim = seg.pdim;          % Number of problem parameters

yidim   = data.yidim;
yicndim = data.yicndim;
% ycndim  = ydim*NCOL*NTST;
xcnnum  = NCOL*NTST;
xcndim  = dim*NCOL*NTST;     % Number of collocation conditions
addim   = xcndim+yicndim+2+pdim;
yiidx   = data.yidx(:);
yiidx   = repmat(yiidx, [1,xcnnum])+repmat((0:xcnnum-1)*yidim, [numel(yiidx),1]);
yiidx   = yiidx(:);

% yiidx   = (data.yidx(1)-1)*xcnnum+1:data.yidx(end)*xcnnum;

data.adim   = [yicndim addim];
data.yiidx  = yiidx;
data.xcnidx = (1:xcndim);
data.ycnidx = xcndim+(1:yicndim);
data.T0idx  = xcndim+yicndim+1;
data.Tidx   = xcndim+yicndim+2;
data.pidx   = xcndim+yicndim+2+(1:pdim);

wts         = seg.wts;
wts         = repmat(wts, [yidim NTST]);
data.wtsy   = spdiags(wts(:), 0, yicndim, yicndim);

data.sd     = seg.Md;   
tsa         = seg.tsa;
tzd         = seg.tzd;
Lda         = coll_L(tsa, tzd);
rows        = reshape(1:yicndim, [yidim*NCOL NTST]);
rows        = repmat(rows, [yidim*NCOL 1]);
cols        = repmat(1:yicndim, [yidim*NCOL 1]);
Wda         = repmat(kron(Lda, eye(yidim)), [1 NTST]); % Interpolation matrix
data.Wda    = sparse(rows, cols, Wda);   


end
