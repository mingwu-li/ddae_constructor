function data = alg_dae_adjt_init_data(prob, data) %#ok<INUSL>
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


ycndim = ydim*NCOL*NTST;    % Number of basepoint values for y
xcndim = dim*NCOL*NTST;     % Number of collocation conditions
addim  = xcndim+ycndim+2+pdim;

wts         = seg.wts;
wts         = repmat(wts, [ydim NTST]);
data.wtsy   = spdiags(wts(:), 0, ycndim, ycndim);
data.adim   = [ycndim addim];

tcn         = seg.Md;
data.dTytcn = repmat(tcn', [ydim 1]);


end
