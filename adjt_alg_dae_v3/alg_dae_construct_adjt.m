function [prob, data] = alg_dae_construct_adjt(prob, tbid, data, sol)
%COLL_CONSTRUCT_ADJT   Add COLL adjoint problem.

% Copyright (C) Frank Schilder, Harry Dankowicz, Mingwu Li
% $Id: coll_construct_opt.m 2872 2015-08-06 20:12:06Z hdankowicz $

% [data, sol] = init_data(data, sol);
opt  = data.ddaecoll_opt;
seg  = data.ddaecoll_seg;
if ~isempty(sol.l0)
    if strcmpi(seg.ddaecoll.Apoints, 'Uniform')
        sol.tbp = sol.tbp(seg.taL_idx);  % delete replicated ones before interpolation
        sol.l0  = sol.l0(:, seg.taL_idx);
    end
    sol.l0  = interp1(sol.tbp', sol.l0', seg.tbpa', 'pchip', 'extrap')';
    sol.l0  = sol.l0(:);
    if ~isempty(sol.tl0)
        if strcmpi(seg.ddaecoll.Apoints, 'Uniform')
            sol.tl0 = sol.tl0(:, seg.taL_idx);
        end        
        sol.tl0 = interp1(sol.tbp', sol.tl0', seg.tbpa', 'pchip', 'extrap')';
        sol.tl0 = sol.tl0(:);
    end
end

% [data, axidx] = coco_get_adjt_data(prob, 'ddaecoll', 'data', 'axidx');

prob = coco_add_adjt(prob, tbid, @adj, @adj_DU, data, 'aidx', ...
    data.axidx([opt.xcn_idx; opt.ycn_idx; opt.T0_idx; opt.T_idx; opt.p_idx]),...
    'l0', sol.l0, 'tl0', sol.tl0, 'adim', data.adim);

end


function [data, J] = adj(prob, data, u) %#ok<INUSL>

seg  = data.ddaecoll_seg;
NTST = seg.ddaecoll.NTST;

x  = u(seg.xbp_idx);
y  = u(seg.ybp_idx);
T0 = u(seg.T0_idx);
T  = u(seg.T_idx);
p  = u(seg.p_idx);

xcnd = reshape(seg.Wdd*x, seg.x_shp); % Values at collocation nodes
ycnd = reshape(seg.Wda*y, seg.y_shp); 
tcnd = seg.Md*T+T0;
xcna = reshape(seg.Wad*x, seg.x_shp); % Values at collocation nodes
ycna = reshape(seg.Waa*y, seg.y_shp); 
tcna = seg.Ma*T+T0;
pcn  = repmat(p, seg.p_rep);


gdtcn = data.dgdthan(tcnd', xcnd, ycnd, pcn);
gdxcn = data.dgdxhan(tcnd', xcnd, ycnd, pcn);
gdycn = data.dgdyhan(tcna', xcna, ycna, pcn);
gdpcn = data.dgdphan(tcnd', xcnd, ycnd, pcn);
 
% adjoint with respect to delta_x
dxode = sparse(data.gdxrows, data.gdxcols, gdxcn(:));
J = seg.Wda'*dxode;

% adjoint with respect to delta_y
dyode = sparse(data.gdyrows, data.gdycols, gdycn(:));
if isempty(dyode)
    J = [ J, [] ];
else
    J = [ J, seg.Waa'*dyode ];
end


% adjoint with respect to T0 and T
dT0ode = gdtcn;
dTode  = gdtcn.*data.dTytcn;
J = [ J, (0.5/NTST)*seg.Wda'*data.wtsy*[dT0ode(:) dTode(:)] ];

% adjoint with respect to p
dpode = gdpcn;
dpode = sparse(data.gdprows, data.gdpcols, dpode(:));
if isempty(dpode)
    J = [ J, [] ];
else
    J = [ J, (0.5/NTST)*seg.Wda'*data.wtsy*dpode ];
end

end

function [data, dJ] = adj_DU(prob, data, u) 

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
