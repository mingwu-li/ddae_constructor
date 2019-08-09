function [prob, data] = alg_dde_construct_adjt(prob, tbid, data, sol)
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

fid  = sprintf('gseg.%s',data.gidx);
fid  = coco_get_id(tbid, fid);
prob = coco_add_adjt(prob, fid, @adj, @adj_DU, data, 'aidx', ...
    data.axidx([opt.xcn_idx; opt.ycn_idx(data.yiidx); opt.T0_idx; opt.T_idx; opt.p_idx]),...
    'l0', sol.l0, 'tl0', sol.tl0, 'adim', data.adim);

end


function [data, J] = adj(prob, data, u) %#ok<INUSL>

seg  = data.ddaecoll_seg;

x  = u(seg.xbp_idx);
y  = u(seg.ybp_idx);
T0 = u(seg.T0_idx);
T  = u(seg.T_idx);
p  = u(seg.p_idx);

NTST  = seg.ddaecoll.NTST;
NCOL  = seg.ddaecoll.NCOL;
tsa   = seg.tsa;
tsd   = seg.tsd;
dim   = seg.dim;
pdim  = seg.pdim;
yidim = data.yidim;
J     = zeros(data.adim);
s     = 0;
for j=1:data.nseg
    gseg = data.gseg{j}; 
    si   = gseg.si(T0, T, p);
    se   = gseg.se(T0, T, p);
    % Find subset of nodes in the interval [sa, sb) or [sa, sb]; we restrict
    % Gauss and Gauss-Radau with left end included
    idx  = intersect(find(data.sa>=si),find(data.sa<se));
    sa   = data.sa(idx);
    if isempty(sa)
        continue;
    end

    Nsa  = numel(sa);
    % A vectorization based implementation
    tj   = T0+T*sa;                                % time at collocation nodes sa
    rowp = (1:yidim*Nsa);
    colp = rowp+(idx(1)-1)*yidim;
    Pj   = sparse(rowp, colp, ones(yidim*Nsa,1), yidim*Nsa, yidim*NTST*NCOL);
    yj   = Pj*data.Waa*y;
    yj   = reshape(yj, [yidim,Nsa]);              % ycn at collocation nodes sa
    pj   = repmat(p, [1,Nsa]);
    Omegaj = Pj*data.wtsy*Pj';
    
    if ~isempty(gseg.taux)
        % shifted node sets and the sets of intervals
        sax   = sa-gseg.taux(T0,T,p)/T;
        idsax = floor(seg.ddaecoll.NTST*sax)+1;    % index for subinterval belongs to
        idsax(idsax==NTST+1) = NTST;
        numid = numel(unique(idsax));
        dimb = numid*dim*(NCOL+1);
        rowt = (1:dimb);
        colt = rowt+(idsax(1)-1)*dim*(NCOL+1);
        Tj   = sparse(rowt, colt, ones(dimb,1), dimb, data.xbpdim);
        xbpj = Tj*x;
        jund = find(data.mg<=sax(1), 1, 'last');
        jbar = find(data.mg>sax(end), 1, 'first');
        if isempty(jbar) % at the right end of whole interval, i.e., 1
            jbar = NTST+1;
        end
        if jund==jbar-1                  % A single subinterval
            tx  = 2*(NTST*sax-(jund-1))-1;     % go to [-1,1]
            Lx  = kron(coll_L(tsd, tx),eye(dim));
            DLx = kron(coll_Lp(tsd, tx),eye(dim)); 
            xj  = Lx*xbpj;
            xj  = reshape(xj, [dim, Nsa]);
        elseif jund==jbar-2              % Two subintervals
            id1 = numel(find(idsax==idsax(1)));
            tx1 = 2*(NTST*sax(1:id1)-(jund-1))-1;
            Lx1 = coll_L(tsd, tx1);           % 1st subinterval
            DLx1 = coll_Lp(tsd, tx1);
            tx2  = 2*(NTST*sax(id1+1:end)-jund)-1;
            Lx2  = coll_L(tsd, tx2);           % 2nd subinterval
            DLx2 = coll_Lp(tsd, tx2);
            Lx  = blkdiag(kron(Lx1,eye(dim)), kron(Lx2,eye(dim)));
            DLx = blkdiag(kron(DLx1,eye(dim)), kron(DLx2,eye(dim)));
            xj  = Lx*xbpj;
            xj  = reshape(xj, [dim, Nsa]);
        elseif jund<jbar-2               % More than two subintervals
            % We use periodicity of collocation nodes to simplify
            % calculation
            nxF = numel(find(idsax==idsax(1)));
            nxL = numel(find(idsax==idsax(end)));
            tss = 2*(NTST*sax(nxF+1:nxF+NCOL)-jund)-1;
            Lxs = coll_L(tsd, tss);           % 2nd - (Last-1) subintervals            
            LxF = Lxs(NCOL-nxF+1:end,:);      % First subintervals
            DLxs = coll_Lp(tsd, tss);  
            DLxF = DLxs(NCOL-nxF+1:end,:); 
            if nxL>NCOL
                tse = 2*(NTST*sax(end-nxL+1:end)-(jbar-2))-1;
                LxL = coll_L(tsd, tse);       % Last subintervals
                DLxL = coll_Lp(tsd, tse);
            else
                LxL = Lxs(1:nxL,:);           % Last subintervals
                DLxL = DLxs(1:nxL,:);
            end
            Lx1 = kron(LxF,eye(dim));
            Lx2 = kron(eye(jbar-jund-2), kron(Lxs,eye(dim)));
            Lx3 = kron(LxL,eye(dim));
            Lx  = blkdiag(Lx1,Lx2,Lx3);
            DLx1 = kron(DLxF,eye(dim));
            DLx2 = kron(eye(jbar-jund-2), kron(DLxs,eye(dim)));
            DLx3 = kron(DLxL,eye(dim));
            DLx  = blkdiag(DLx1,DLx2,DLx3);           
            xj  = Lx*xbpj;
            xj  = reshape(xj, [dim, Nsa]);                
        else
            disp('bar{j} and underline{j} is not matched well');
        end
        DLx = 2*NTST*DLx*Tj;
        
        
        % Find subset of nodes in the interval [sa, sb) or [sa, sb]; we restrict
        % Gauss and Gauss-Radau with left end included
        jdx  = intersect(find(data.sd>=si-gseg.taux(T0,T,p)/T),find(data.sd<se-gseg.taux(T0,T,p)/T));
        sd   = data.sd(jdx);
        if isempty(sd)
            continue;
        end
        Nsd = numel(sd);
        tjd = T0+T*(sd+gseg.taux(T0,T,p)/T);
        rowp = (1:dim*Nsd);
        colp = rowp+(jdx(1)-1)*dim;
        Pjd   = sparse(rowp, colp, ones(dim*Nsd,1), dim*Nsd, dim*NTST*NCOL);
        xjd   = Pjd*seg.Wdd*x;
        xjd   = reshape(xjd, [dim,Nsd]);              % ycn at collocation nodes sa
        pjd   = repmat(p, [1,Nsd]);        
       
        
        sdx   = sd+gseg.taux(T0,T,p)/T;
        idsdx = floor(seg.ddaecoll.NTST*sdx)+1;    % index for subinterval belongs to
        idsdx(idsdx==NTST+1) = NTST;
        numid = numel(unique(idsdx));
        dimb = numid*yidim*NCOL;
        rowt = (1:dimb);
        colt = rowt+(idsdx(1)-1)*yidim*NCOL;
        Tjd  = sparse(rowt, colt, ones(dimb,1), dimb, data.yicndim);
        ybpj = Tjd*y;
        jund = find(data.mg<=sdx(1), 1, 'last');
        jbar = find(data.mg>sdx(end), 1, 'first');
        if isempty(jbar) % at the right end of whole interval, i.e., 1
            jbar = NTST+1;
        end
        if jund==jbar-1                  % A single subinterval
            ty  = 2*(NTST*sdx-(jund-1))-1;     % go to [-1,1]
            Ly  = kron(coll_L(tsa, ty),eye(yidim));
            yjd = Ly*ybpj;
            yjd = reshape(yjd, [yidim, Nsd]);
        elseif jund==jbar-2              % Two subintervals
            id1 = numel(find(idsdx==idsdx(1)));
            ty1 = 2*(NTST*sdx(1:id1)-(jund-1))-1;
            Ly1 = coll_L(tsa, ty1);           % 1st subinterval
            ty2  = 2*(NTST*sdx(id1+1:end)-jund)-1;
            Ly2  = coll_L(tsa, ty2);           % 2nd subinterval
            Ly  = blkdiag(kron(Ly1,eye(yidim)), kron(Ly2,eye(yidim)));
            yjd = Ly*ybpj;
            yjd = reshape(yjd, [yidim, Nsd]);
        elseif jund<jbar-2               % More than two subintervals
            % We use periodicity of collocation nodes to simplify
            % calculation
            nyF = numel(find(idsdx==idsdx(1)));
            nyL = numel(find(idsdx==idsdx(end)));
            tss = 2*(NTST*sdx(nyF+1:nyF+NCOL)-jund)-1;
            Lys = coll_L(tsa, tss);           % 2nd - (Last-1) subintervals            
            LyF = Lys(NCOL-nyF+1:end,:);      % First subintervals
            if nyL>NCOL
                tse = 2*(NTST*sdx(end-nxL+1:end)-(jbar-2))-1;
                LyL = coll_L(tsa, tse);       % Last subintervals
            else
                LyL = Lys(1:nyL,:);           % Last subintervals
            end
            Ly1 = kron(LyF,eye(yidim));
            Ly2 = kron(eye(jbar-jund-2), kron(Lys,eye(yidim)));
            Ly3 = kron(LyL,eye(yidim));
            Ly  = blkdiag(Ly1,Ly2,Ly3);          
            yjd = Ly*ybpj;
            yjd = reshape(yjd, [yidim, Nsd]);                
        else
            disp('bar{j} and underline{j} is not matched well');
        end        
        
        
        % Differentiation of tau_{j}^x 
        DtauxdT0han = gseg.Dtaux{1};
        DtauxdThan  = gseg.Dtaux{2};
        Dtauxdphan  = gseg.Dtaux{3};
        % w.r.t. T0
        dtauxdT0 = DtauxdT0han(T0, T, p);
        % w.r.t. T
        dtauxdT  = DtauxdThan(T0, T, p);
        % w.r.t. p
        dtauxdp  = Dtauxdphan(T0,T, p);
        
        if ~isempty(gseg.phi)
            tphi = tj-gseg.taut(T0,T,p);
            phij = gseg.phi(tphi');
            phidim  = size(phij,1);
            
            tphid = tjd-gseg.taut(T0,T,p);
            phijd = gseg.phi(tphid');
            
            % Differentiation of tau_{j}^t 
            DtautdT0han = gseg.Dtaut{1};
            DtautdThan  = gseg.Dtaut{2};
            Dtautdphan  = gseg.Dtaut{3};
            % w.r.t. T0
            dtautdT0 = DtautdT0han(T0, T, p);
            % w.r.t. T
            dtautdT  = DtautdThan(T0, T, p);
            % w.r.t. p
            dtautdp  = Dtautdphan(T0,T, p);
            
            % Differenation of phi
            Dphij = gseg.Dphi(tphi');
            
            % Differention of g
            dgdthan = gseg.Dg{1};
            dgdyhan = gseg.Dg{2};
            dgdxhan = gseg.Dg{3};
            dgdphihan = gseg.Dg{4};
            dgdphan = gseg.Dg{5};
            
            % w.r.t t
            dtode  = dgdthan(tj', yj, xj, phij, pj);
            % w.r.t. yi & yibp
            dyode  = dgdyhan(tj', yj, xj, phij, pj);
            dyrows = repmat(reshape(1:Nsa*yidim, [yidim Nsa]), [yidim 1]); 
            dycols = repmat(1:Nsa*yidim, [yidim 1]);
            dyode  = sparse(dyrows, dycols, dyode(:));
             % w.r.t x & xbp
            dxode  = dgdxhan(tj', yj, xj, phij, pj);
            dxrows = repmat(reshape(1:Nsa*yidim, [yidim Nsa]), [dim 1]); 
            dxcols = repmat(1:Nsa*dim, [yidim 1]);  
            dxode  = sparse(dxrows, dxcols, dxode(:));         
            % w.r.t. phi           
            dphiode  = dgdphihan(tj', yj, xj, phij, pj);
            dphirows = repmat(reshape(1:Nsa*yidim, [yidim Nsa]), [phidim 1]); 
            dphicols = repmat(1:Nsa*phidim, [yidim 1]);  
            dphiode  = sparse(dphirows, dphicols, dphiode(:));             
            % w.r.t. p
            dpode  = dgdphan(tj', yj, xj, phij, pj);
            dprows = repmat(reshape(1:Nsa*yidim, [yidim Nsa]), [pdim 1]); 
            dpcols = repmat(1:pdim, [yidim Nsa]);  
            dpode  = sparse(dprows, dpcols, dpode(:));             
 
            % w.r.t x & xbp
            dxdode = dgdxhan(tjd', yjd, xjd, phijd, pjd);
            dxrows = repmat(reshape(1:Nsd*yidim, [yidim Nsd]), [yidim 1]); 
            dxcols = repmat(1:Nsd*yidim, [yidim 1]);  
            dxdode = sparse(dxrows, dxcols, dxdode(:)); 
            
            % variation w.r.t. x & xbp
            J_xcn = (Ly*Tjd)'*dxdode;
            % variation w.r.t. y & ybp
            J_ycn = (Pj*data.Waa)'*dyode;
            % variation w.r.t. T0
            J_T0  = dtode(:)-dxode*DLx*x*dtauxdT0/T+dphiode*Dphij(:)*(1-dtautdT0);
            J_T0  = (Pj*data.Waa)'*Omegaj*J_T0/(2*NTST);
            % variation w.r.t. T
            dtode = dtode.*repmat(sa', [yidim,1]);
            duode = Dphij.*repmat(sa'-dtautdT, [phidim,1]);
            J_T   = dtode(:)-dxode*DLx*x*(dtauxdT*T-gseg.taux(T0,T,p))/T^2+dphiode*duode(:);
            J_T   = (Pj*data.Waa)'*Omegaj*J_T/(2*NTST);
            % variation w.r.t. p
            J_p   = dpode-dxode*DLx*x*dtauxdp/T-dphiode*Dphij(:)*dtautdp;
            J_p   = (Pj*data.Waa)'*Omegaj*J_p/(2*NTST);
            Jj    = [J_xcn, J_ycn, J_T0, J_T, J_p];
            
            xcnj_idx = (jdx(1)-1)*yidim+1:jdx(end)*yidim;
            ycnj_idx = (idx(1)-1)*dim+1:idx(end)*dim;
            colidx   = [data.xcnidx(xcnj_idx); data.ycnidx(ycnj_idx); data.T0idx; data.Tidx; data.pidx];
            J(:,colidx) = J(:,colidx) + Jj;
            
%             J_xbp = dxode*Lx*Tj;
%             J_ybp = dyode*Pj*data.Waa;
%             J_T0  = dtode(:)-dxode*DLx*x*dtauxdT0/T+dphiode*Dphij(:)*(1-dtautdT0);
%             dtode = dtode.*repmat(sa, [yidim,1]);
%             duode = Dphij.*repmat(sa-dtautdT, [phidim,1]);
%             J_T   = dtode(:)-dxode*DLx*x*(dtauxdT*T-gseg.taux(T0,T,p))/T^2+dphiode*duode(:);
%             J_p   = dpode-dxode*DLx*x*dtauxdp/T-dphiode*Dphij(:)*dtautdp;
%             
%             J(s+1:s+Nsa*yidim,:) = [J_xbp J_ybp J_T0, J_T, J_p];
        else            
            % Differention of g
            dgdthan = gseg.Dg{1};
            dgdyhan = gseg.Dg{2};
            dgdxhan = gseg.Dg{3};
            dgdphan = gseg.Dg{4};
            
            % w.r.t t
            dtode  = dgdthan(tj', yj, xj, pj);
            % w.r.t. yi & yibp
            dyode  = dgdyhan(tj', yj, xj, pj);
            dyrows = repmat(reshape(1:Nsa*yidim, [yidim Nsa]), [yidim 1]); 
            dycols = repmat(1:Nsa*yidim, [yidim 1]);
            dyode  = sparse(dyrows, dycols, dyode(:));
             % w.r.t x & xbp
            dxode  = dgdxhan(tj', yj, xj, pj);
            dxrows = repmat(reshape(1:Nsa*yidim, [yidim Nsa]), [dim 1]); 
            dxcols = repmat(1:Nsa*dim, [yidim 1]);  
            dxode  = sparse(dxrows, dxcols, dxode(:));         
            % w.r.t. p
            dpode  = dgdphan(tj', yj, xj, pj);
            dprows = repmat(reshape(1:Nsa*yidim, [yidim Nsa]), [pdim 1]); 
            dpcols = repmat(1:pdim, [yidim Nsa]);  
            dpode  = sparse(dprows, dpcols, dpode(:));             

            % w.r.t x & xbp
            dxdode = dgdxhan(tjd', yjd, xjd, pjd);
            dxrows = repmat(reshape(1:Nsd*yidim, [yidim Nsd]), [yidim 1]); 
            dxcols = repmat(1:Nsd*yidim, [yidim 1]);  
            dxdode = sparse(dxrows, dxcols, dxdode(:));             
            
            % variation w.r.t. x & xbp
            J_xcn = (Ly*Tjd)'*dxdode;
            % variation w.r.t. y & ybp
            J_ycn = (Pj*data.Waa)'*dyode;
            % variation w.r.t. T0
            J_T0  = dtode(:)-dxode*DLx*x*dtauxdT0/T;
            J_T0  = (Pj*data.Waa)'*Omegaj*J_T0/(2*NTST);
            % variation w.r.t. T
            dtode = dtode.*repmat(sa', [yidim,1]);
            J_T   = dtode(:)-dxode*DLx*x*(dtauxdT*T-gseg.taux(T0,T,p))/T^2;
            J_T   = (Pj*data.Waa)'*Omegaj*J_T/(2*NTST);
            % variation w.r.t. p
            J_p   = dpode-dxode*DLx*x*dtauxdp/T;
            J_p   = (Pj*data.Waa)'*Omegaj*J_p/(2*NTST);
            
            Jj    = [J_xcn, J_ycn, J_T0, J_T, J_p];
            
            xcnj_idx = (jdx(1)-1)*dim+1:jdx(end)*dim;
            ycnj_idx = (idx(1)-1)*yidim+1:idx(end)*yidim;
            colidx   = [data.xcnidx(xcnj_idx), data.ycnidx(ycnj_idx), data.T0idx, data.Tidx, data.pidx];
            J(:,colidx) = J(:,colidx) + Jj;
            
%             J(s+1:s+Nsa*yidim,:) = [J_xbp J_ybp J_T0, J_T, J_p];
            
%             f(s+1:s+Nsa*yidim) = gseg.gfunc(tj',yj,xj,pj);
        end 
    else
        if ~isempty(gseg.phi)         
            tphi = tj-gseg.taut(T0,T,p);
            phij = gseg.phi(tphi');
            phidim  = size(phij,1);
            
            % Differentiation of tau_{j}^t 
            DtautdT0han = gseg.Dtaut{1};
            DtautdThan  = gseg.Dtaut{2};
            Dtautdphan  = gseg.Dtaut{3};
            % w.r.t. T0
            dtautdT0 = DtautdT0han(T0, T, p); 
            % w.r.t. T
            dtautdT = DtautdThan(T0, T, p);
            % w.r.t. p
            dtautdp = Dtautdphan(T0,T, p);
            
            % Differenation of phi
            Dphij = gseg.Dphi(tphi');
            
            % Differention of g
            dgdthan = gseg.Dg{1};
            dgdyhan = gseg.Dg{2};
            dgdphihan = gseg.Dg{3};
            dgdphan = gseg.Dg{4};
            
            % w.r.t t
            dtode = dgdthan(tj', yj, phij, pj);
            % w.r.t. yi & yibp
            dyode  = dgdyhan(tj', yj, phij, pj);
            dyrows = repmat(reshape(1:Nsa*yidim, [yidim Nsa]), [yidim 1]); 
            dycols = repmat(1:Nsa*yidim, [yidim 1]);
            dyode  = sparse(dyrows, dycols, dyode(:));     
            % w.r.t. phi           
            dphiode  = dgdphihan(tj', yj, phij, pj);
            dphirows = repmat(reshape(1:Nsa*yidim, [yidim Nsa]), [phidim 1]); 
            dphicols = repmat(1:Nsa*phidim, [yidim 1]);  
            dphiode  = sparse(dphirows, dphicols, dphiode(:));             
            % w.r.t. p
            dpode  = dgdphan(tj', yj, phij, pj);
            dprows = repmat(reshape(1:Nsa*yidim, [yidim Nsa]), [pdim 1]); 
            dpcols = repmat(1:pdim, [yidim Nsa]);  
            dpode  = sparse(dprows, dpcols, dpode(:));             
            
            % variation w.r.t. x & xbp
            dxode = sparse(Nsa*yidim, numel(x));
            J_xcn = (Ly*Tjd)'*dxode;
            % variation w.r.t. y & ybp
            J_ycn = (Pj*data.Waa)'*dyode;
            % variation w.r.t. T0
            J_T0  = dtode(:)+dphiode*Dphij(:)*(1-dtautdT0);
            J_T0  = (Pj*data.Waa)'*Omegaj*J_T0/(2*NTST);
            % variation w.r.t. T
            dtode = dtode.*repmat(sa', [yidim,1]);
            duode = Dphij.*repmat(sa'-dtautdT, [phidim,1]);
            J_T   = dtode(:)+dphiode*duode(:);
            J_T   = (Pj*data.Waa)'*Omegaj*J_T/(2*NTST);
            % variation w.r.t. p
            J_p   = dpode-dphiode*Dphij(:)*dtautdp;
            J_p   = (Pj*data.Waa)'*Omegaj*J_p/(2*NTST);
            Jj    = [J_xcn, J_ycn, J_T0, J_T, J_p];
            
            xcnj_idx = (jdx(1)-1)*dim+1:jdx(end)*dim;
            ycnj_idx = (idx(1)-1)*yidim+1:idx(end)*yidim;
            colidx   = [data.xcnidx(xcnj_idx), data.ycnidx(ycnj_idx), data.T0idx, data.Tidx, data.pidx];
            J(:,colidx) = J(:,colidx) + Jj;            
            
            
%             J_xbp = sparse(Nsa*yidim, numel(x));
%             J_ybp = dyode*Pj*data.Waa;
%             J_T0  = dtode(:)+dphiode*Dphij(:)*(1-dtautdT0);
%             dt1   = repmat(sa', [yidim,1]);
% %             dtode = dtode.*repmat(sa, [yidim,1]);
%             duode = Dphij.*repmat(sa'-dtautdT, [phidim,1]);
%             J_T   = dtode(:).*dt1(:)+dphiode*duode(:);
%             J_p   = dpode-dphiode*Dphij(:)*dtautdp;
%             
%             J(s+1:s+Nsa*yidim,:) = [J_xbp J_ybp J_T0, J_T, J_p];
                      
%            f(s+1:s+Nsa*yidim) = gseg.gfunc(tj',yj,phij,pj);
        else       
            % Differention of g
            dgdthan = gseg.Dg{1};
            dgdyhan = gseg.Dg{2};
            dgdphan = gseg.Dg{3};
            
            % w.r.t t
            dtode = dgdthan(tj', yj, pj);
            % w.r.t. yi & yibp
            dyode = dgdyhan(tj', yj, pj);
            dyrows = repmat(reshape(1:Nsa*yidim, [yidim Nsa]), [yidim 1]); 
            dycols = repmat(1:Nsa*yidim, [yidim 1]);
            dyode  = sparse(dyrows, dycols, dyode(:));                 
            % w.r.t. p
            dpode = dgdphan(tj', yj, pj);
            dprows = repmat(reshape(1:Nsa*yidim, [yidim Nsa]), [pdim 1]); 
            dpcols = repmat(1:pdim, [yidim Nsa]);  
            dpode  = sparse(dprows, dpcols, dpode(:));             

            % variation w.r.t. x & xbp
            dxode = sparse(Nsa*yidim, numel(x));
            J_xcn = (Ly*Tjd)'*dxode;
            % variation w.r.t. y & ybp
            J_ycn = (Pj*data.Waa)'*dyode;
            % variation w.r.t. T0
            J_T0  = dtode(:);
            J_T0  = (Pj*data.Waa)'*Omegaj*J_T0/(2*NTST);
            % variation w.r.t. T
            dtode = dtode.*repmat(sa', [yidim,1]);
            J_T   = dtode(:);
            J_T   = (Pj*data.Waa)'*Omegaj*J_T/(2*NTST);
            % variation w.r.t. p
            J_p   = dpode;
            J_p   = (Pj*data.Waa)'*Omegaj*J_p/(2*NTST);
            Jj    = [J_xcn, J_ycn, J_T0, J_T, J_p];
            
            xcnj_idx = (jdx(1)-1)*dim+1:jdx(end)*dim;
            ycnj_idx = (idx(1)-1)*yidim+1:idx(end)*yidim;
            colidx   = [data.xcnidx(xcnj_idx), data.ycnidx(ycnj_idx), data.T0idx, data.Tidx, data.pidx];
            J(:,colidx) = J(:,colidx) + Jj;                
            
%             J_xbp = sparse(Nsa*yidim, numel(x));
%             J_ybp = dyode*Pj*data.Waa;
%             J_T0  = dtode(:);
%             dtode = dtode.*repmat(sa, [yidim,1]);
%             J_T   = dtode(:);
%             J_p   = dpode;
%             
%             J(s+1:s+Nsa*yidim,:) = [J_xbp J_ybp J_T0, J_T, J_p];
%             
%             f(s+1:s+Nsa*yidim) = gseg.gfunc(tj',yj,pj);
        end
    end
    s = s+Nsa*yidim;
%     toc
end

% % adjoint with respect to delta_x
% dxode = sparse(data.gdxrows, data.gdxcols, gdxcn(:));
% J = seg.Wda'*dxode;
% 
% % adjoint with respect to delta_y
% dyode = sparse(data.gdyrows, data.gdycols, gdycn(:));
% if isempty(dyode)
%     J = [ J, [] ];
% else
%     J = [ J, seg.Waa'*dyode ];
% end
% 
% 
% % adjoint with respect to T0 and T
% dT0ode = gdtcn;
% dTode  = gdtcn.*data.dTytcn;
% J = [ J, (0.5/NTST)*seg.Wda'*data.wtsy*[dT0ode(:) dTode(:)] ];
% 
% % adjoint with respect to p
% dpode = gdpcn;
% dpode = sparse(data.gdprows, data.gdpcols, dpode(:));
% if isempty(dpode)
%     J = [ J, [] ];
% else
%     J = [ J, (0.5/NTST)*seg.Wda'*data.wtsy*dpode ];
% end

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
