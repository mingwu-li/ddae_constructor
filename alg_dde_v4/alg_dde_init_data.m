function data = alg_dde_init_data(data)
%COLL_INIT_DATA   Initialize toolbox data for an instance of 'coll'.
%
% Populate remaining fields of the toolbox data structure used by 'coll'
% function objects.
%
% DATA = COLL_INIT_DATA(DATA, X0, P0)
%
% DATA - Toolbox data structure.
% X0   - Initial solution guess for discretized trajectory.
% P0   - Initial solution guess for problem parameters.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_init_data.m 2839 2015-03-05 17:09:01Z fschild $

seg  = data.ddaecoll_seg;
NTST = seg.ddaecoll.NTST; % Number of mesh intervals
NCOL = seg.ddaecoll.NCOL; % Degree of polynomial interpolants
dim  = seg.dim;    % State-space dimension
ydim = seg.ydim;
pdim = seg.pdim;      % Number of problem parameters

% ny = ydim/dim;

% xcndim = dim*NCOL*NTST;     % Number of collocation conditions
xcnnum = NCOL*NTST;
xbpdim = dim*(NCOL+1)*NTST;

yidim  = numel(data.yidx);
ycndim = yidim*NCOL*NTST;

data.xbpdim  = xbpdim;
data.mg      = (0:NTST)/NTST;
data.yidim   = yidim;
data.yicndim = ycndim;
data.y_shp   = [yidim, xcnnum];
data.xbp_idx = seg.xbp_idx;
data.ybp_idx = xbpdim+(1:ycndim)';
data.T0_idx  = xbpdim+ycndim+1;
data.T_idx   = xbpdim+ycndim+2;
data.p_idx   = xbpdim+ycndim+2+(1:pdim)';
yiidx        = repmat(data.yidx(:), [1, xcnnum]);
yiidx        = yiidx+repmat((0:xcnnum-1)*ydim, [yidim, 1]);
yiidx        = seg.ybp_idx(yiidx(:));
guidx        = [data.xbp_idx; yiidx; seg.T0_idx; seg.T_idx; seg.p_idx];
data.guidx   = data.ddae_uidx(guidx);

data.sa     = seg.Ma;   
tsa         = seg.tsa;
tza         = seg.tza;
Laa         = coll_L(tsa, tza);
rows        = reshape(1:ycndim, [yidim*NCOL NTST]);
rows        = repmat(rows, [yidim*NCOL 1]);
cols        = repmat(1:ycndim, [yidim*NCOL 1]);
Waa         = repmat(kron(Laa, eye(yidim)), [1 NTST]); % Interpolation matrix
data.Waa    = sparse(rows, cols, Waa);   


% (data.yidx-1)*xcndim+data.ybp_idx
% ycndim = ydim*NCOL*NTST;    % Number of basepoint values for y
% ycnnum = NCOL*NTST;         % Number of collocation nodes
% xcndim = dim*NCOL*NTST;     % Number of collocation conditions
% 
% data.ycnnum = ycnnum;
% 
% 
% 
% data.xbpnum  = xbpnum;
% data.xcnnum  = xcnnum;
% data.xbp_idx = (1:xbpdim)'; % Index array for basepoint values
% data.ybp_idx = xbpdim+(1:ybpdim)';
% data.T0_idx  = xbpdim+ybpdim+1;
% data.T_idx   = xbpdim+ybpdim+2;    % Index for interval length


% data.gdxrows = repmat(reshape(1:ycndim, [ydim ycnnum]), [dim 1]);   % Index array for vectorization
% data.gdxcols = repmat(1:xcndim, [ydim 1]);                          % Index array for vectorization
% data.gdprows = repmat(reshape(1:ycndim, [ydim ycnnum]), [pdim 1]);  % Index array for vectorization
% data.gdpcols = repmat(1:pdim, [ydim ycnnum]);                       % Index array for vectorization
% data.gdyrows = repmat(reshape(1:ycndim, [ydim ycnnum]), [ydim 1]);  % Index array for vectorization
% data.gdycols = repmat(1:ycndim, [ydim 1]);                          % Index array for vectorization
% data.gdtrows = reshape(1:ycndim, [ydim ycnnum]);                    % Index array for vectorization
% data.gdtcols = repmat(1:ycnnum, [ydim 1]);                          % Index array for vectorization

end


