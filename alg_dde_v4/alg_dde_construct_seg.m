function prob = alg_dde_construct_seg(prob, tbid, data)
%COLL_CONSTRUCT_SEG   Append an instance of 'coll' to problem.
%
% Add collocation and continuity conditions, monitor functions that
% evaluate to the problem parameters, and corresponding inactive
% continuation parameters.
%
% PROB = COLL_CONSTRUCT_SEG(PROB, TBID, DATA, SOL)
%
% PROB - Continuation problem structure.
% TBID - Toolbox instance identifier.
% DATA - Toolbox data structure.
% SOL  - Initial solution guess.

% Copyright (C) Frank Schilder, Harry Dankowicz
% $Id: coll_construct_seg.m 2839 2015-03-05 17:09:01Z fschild $

assert(numel(data.ddae_uidx)==data.ddaecoll_seg.T_idx+numel(data.ddaecoll_seg.p_idx), ...
    'size of uidx do not match');

data.gseg = cell(data.nseg,1);
for j=1:data.nseg
    gseg.si = data.sj{j};
    gseg.se = data.sj{j+1};
    gseg.id = j;
    gj      = data.g{j};
    gseg.gfunc = gj{1};
    gseg.taux  = [];
    gseg.taut  = [];
    gseg.phi   = [];
    gseg.Dg    = [];
    gseg.Dtaux = [];
    gseg.Dtaut = [];
    gseg.Dphi  = [];
    ngarg      = numel(gj);
    if ngarg>=2
        gseg.taux = gj{2};
        if ngarg>=3
            gseg.taut = gj{3};
            if ngarg>=4
                gseg.phi = gj{4};
                if ngarg>=5
                    gseg.Dg = gj{5};
%                     Dg      = gj{5};
%                     nDg     = numel(Dg);
%                     gseg.Dg = cell(5,1);
%                     switch nDg
%                         case 1
%                             gseg.Dg{1}   = Dg;
%                         case 2
%                             gseg.Dg(1:2) = Dg;
%                         case 3
%                             gseg.Dg(1:3) = Dg;
%                         case 4
%                             gseg.Dg(1:4) = Dg;
%                         case 5
%                             gseg.Dg(1:5) = Dg;
%                     end
                    if ngarg>=6
                        gseg.Dtaux = gj{6};
%                         Dtaux      = gj{6};
%                         nDtaux     = numel(Dtaux);
%                         gseg.Dtaux = cell(3,1);
%                         switch nDtaux
%                             case 1
%                                 gseg.Dtaux{1}   = Dtaux;
%                             case 2
%                                 gseg.Dtaux(1:2) = Dtaux;
%                             case 3
%                                 gseg.Dtaux(1:3) = Dtaux;
%                         end                        
                        if ngarg>=7
                            gseg.Dtaut = gj{7};
%                             Dtaut      = gj{7};
%                             nDtaut     = numel(Dtaut);
%                             gseg.Dtaut = cell(3,1);
%                             switch nDtaut
%                                 case 1
%                                     gseg.Dtaut{1}   = Dtaut;
%                                 case 2
%                                     gseg.Dtaut(1:2) = Dtaut;
%                                 case 3
%                                     gseg.Dtaut(1:3) = Dtaut;
%                             end  
                            if ngarg>=8
                                gseg.Dphi = gj{8};
                            end
                        end
                    end
                end
            end
        end 
    end
    data.gseg{j} = gseg;
end

fid  = sprintf('gseg.%s',data.gidx);
fid  = coco_get_id(tbid, fid);
prob = coco_add_func(prob, fid, @alg_dde_F, @alg_dde_DFDU, data, 'zero', 'uidx', data.guidx);
prob = coco_add_slot(prob, fid, @coco_save_data, data, 'save_full');

end

function [data f] = alg_dde_F(prob, data, u)
%COLL_F   Collocation zero problem.
%
% Expects vectorized encoding of vector field as function of two arguments.

seg = data.ddaecoll_seg;
x   = u(data.xbp_idx); % Extract basepoint values
y   = u(data.ybp_idx);
T0  = u(data.T0_idx);
T   = u(data.T_idx);   % Extract interval length
p   = u(data.p_idx);   % Extract problem parameters

xx  = reshape(x, seg.xbp_shp); % Values at base nodes
% pp  = repmat(p, seg.p_rep);
yy  = reshape(y, data.y_shp); 
% tcn = seg.Wta*(seg.Ma*T+T0);

ycn   = reshape(data.Waa*y, data.y_shp);

NTST  = seg.ddaecoll.NTST;
NCOL  = seg.ddaecoll.NCOL;
% tsa   = seg.tsa;
tsd   = seg.tsd;
dim   = seg.dim;
yidim = data.yidim;
f     = zeros(data.yicndim,1);
s     = 0;
% yref  = zeros(data.yicndim,1);
% sref  = 0;
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
    % Find the index of collocation nodes (1,...,Nm) of sa(1)
%     idsa = floor(seg.ddaecoll.NTST*sa)+1; % index for subinterval belongs to
%     idj1 = numel(find(idsa==idsa(1)));
%     idj1 = idsa(1)*NCOL+1-idj1;           % index of sa(1)
%     idjend = idj1+numel(idx)-1;           % index of sa(end)

    Nsa  = numel(sa);
    % A vectorization based implementation
    tj   = T0+T*sa;                                % time at collocation nodes sa
%     yj   = ycn(:,idx);                             % ycn at collocation nodes sa
    rowp = (1:yidim*Nsa);
    colp = rowp+(idx(1)-1)*yidim;
    Pj   = sparse(rowp, colp, ones(yidim*Nsa,1), yidim*Nsa, yidim*NTST*NCOL);
    yj   = Pj*data.Waa*y;
    yj   = reshape(yj, [yidim,Nsa]);
    pj   = repmat(p, [1,Nsa]);
    
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
%         xbpj = xx(:,(idsax(1)-1)*(NCOL+1)+1:idsax(end)*(NCOL+1));
        jund = find(data.mg<=sax(1), 1, 'last');
        jbar = find(data.mg>sax(end), 1, 'first');
        if isempty(jbar) % at the right end of whole interval, i.e., 1
            jbar = NTST+1;
        end
        if jund==jbar-1                  % A single subinterval
            tx = 2*(NTST*sax-(jund-1))-1;     % go to [-1,1]
            Lx = kron(coll_L(tsd, tx),eye(dim));
            xj = Lx*xbpj;
            xj = reshape(xj, [dim, Nsa]);
%             xj = coll_L(tsd, tx)*xbpj';
%             xj = xj';
        elseif jund==jbar-2              % Two subintervals
            id1 = numel(find(idsax==idsax(1)));
            tx1 = 2*(NTST*sax(1:id1)-(jund-1))-1;
            Lx1 = coll_L(tsd, tx1);           % 1st subinterval
            tx2 = 2*(NTST*sax(id1+1:end)-jund)-1;
            Lx2 = coll_L(tsd, tx2);           % 2nd subinterval
            Lx  = blkdiag(kron(Lx1,eye(dim)), kron(Lx2,eye(dim)));
            xj  = Lx*xbpj;
            xj  = reshape(xj, [dim, Nsa]);
%             Lx  = blkdiag(Lx1, Lx2);
%             xj  = Lx*xbpj';
%             xj  = xj';
        elseif jund<jbar-2               % More than two subintervals
            % We use periodicity of collocation nodes to simplify
            % calculation
            nxF = numel(find(idsax==idsax(1)));
            nxL = numel(find(idsax==idsax(end)));
            tss = 2*(NTST*sax(nxF+1:nxF+NCOL)-jund)-1;
            Lxs = coll_L(tsd, tss);           % 2nd - (Last-1) subintervals            
            LxF = Lxs(NCOL-nxF+1:end,:);      % First subintervals
            if nxL>NCOL
                tse = 2*(NTST*sax(end-nxL+1:end)-(jbar-2))-1;
%                 LxL = [Lxs; Lextra];          % Last subintervals
                LxL = coll_L(tsd, tse);       % Last subintervals
            else
                LxL = Lxs(1:nxL,:);           % Last subintervals
            end
            Lx1 = kron(LxF,eye(dim));
            Lx2 = kron(eye(jbar-jund-2), kron(Lxs,eye(dim)));
            Lx3 = kron(LxL,eye(dim));
            Lx  = blkdiag(Lx1,Lx2,Lx3);
            xj  = Lx*xbpj;
            xj  = reshape(xj, [dim, Nsa]);                
            
%             Lx  = blkdiag(LxF, kron(eye(jbar-jund-2), Lxs), LxL);
%             xj  = Lx*xbpj';
%             xj  = xj';
        else
            disp('bar{j} and underline{j} is not matched well');
        end
        
        if ~isempty(gseg.phi)
            tphi = tj-gseg.taut(T0,T,p);
            phij = gseg.phi(tphi');
            f(s+1:s+Nsa*yidim) = gseg.gfunc(tj',yj,xj,phij,pj);
        else
            f(s+1:s+Nsa*yidim) = gseg.gfunc(tj',yj,xj,pj);
        end 
    else
        if ~isempty(gseg.phi)
           tphi = tj-gseg.taut(T0,T,p);
           phij = gseg.phi(tphi');
           f(s+1:s+Nsa*yidim) = gseg.gfunc(tj',yj,phij,pj);
        else
           f(s+1:s+Nsa*yidim) = gseg.gfunc(tj',yj,pj);
        end
    end
    s = s+Nsa*yidim;
%     toc
end

assert(s==data.yicndim, 'Number of collocation equations is incorret in %s', data.gidx);

end

function [data J] = alg_dde_DFDU(prob, data, u)
%COLL_DFDU   Linearization of collocation zero problem.
%
% Expects vectorized encoding of Jacobians of vector field with respect to
% its arguments. When dfdxhan and/or dfdphan are empty, approximate
% Jacobians are obtained using numerical differentiation.
% tic
seg = data.ddaecoll_seg;
x   = u(data.xbp_idx); % Extract basepoint values
y   = u(data.ybp_idx);
T0  = u(data.T0_idx);
T   = u(data.T_idx);   % Extract interval length
p   = u(data.p_idx);   % Extract problem parameters

% xx  = reshape(x, seg.xbp_shp); % Values at base nodes
% pp  = repmat(p, seg.p_rep);
% yy  = reshape(y, data.y_shp); 
% tcn = seg.Wta*(seg.Ma*T+T0);

% ycn   = reshape(data.Waa*y, data.y_shp);

NTST  = seg.ddaecoll.NTST;
NCOL  = seg.ddaecoll.NCOL;
% tsa   = seg.tsa;
tsd   = seg.tsd;
dim   = seg.dim;
pdim  = seg.pdim;
yidim = data.yidim;
J     = zeros(data.yicndim,numel(u));
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
            tx = 2*(NTST*sax-(jund-1))-1;     % go to [-1,1]
            Lx = kron(coll_L(tsd, tx),eye(dim));
            DLx = kron(coll_Lp(tsd, tx),eye(dim)); 
            xj = Lx*xbpj;
            xj = reshape(xj, [dim, Nsa]);
        elseif jund==jbar-2              % Two subintervals
            id1 = numel(find(idsax==idsax(1)));
            tx1 = 2*(NTST*sax(1:id1)-(jund-1))-1;
            Lx1 = coll_L(tsd, tx1);           % 1st subinterval
            DLx1 = coll_Lp(tsd, tx1);
            tx2 = 2*(NTST*sax(id1+1:end)-jund)-1;
            Lx2 = coll_L(tsd, tx2);           % 2nd subinterval
            DLx2 = coll_Lp(tsd, tx2);
            Lx  = blkdiag(kron(Lx1,eye(dim)), kron(Lx2,eye(dim)));
            DLx  = blkdiag(kron(DLx1,eye(dim)), kron(DLx2,eye(dim)));
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
        
        % Differentiation of tau_{j}^x 
        DtauxdT0han = [];
        DtauxdThan  = [];
        Dtauxdphan  = [];
        switch numel(gseg.Dtaux)
            case 1
                DtauxdT0han = gseg.Dtaux{1};
            case 2
                DtauxdT0han = gseg.Dtaux{1};
                DtauxdThan  = gseg.Dtaux{2};
            case 3
                DtauxdT0han = gseg.Dtaux{1};
                DtauxdThan  = gseg.Dtaux{2};
                Dtauxdphan  = gseg.Dtaux{3};
        end
        % w.r.t. T0
        if isempty(DtauxdT0han)
          f   = @(x,p) gseg.taux(x, p(1), p(2:end));
          Tp  = [T; p];
          dtauxdT0 = coco_ezDFDX('f(x,p)', f, T0, Tp);
        else
          dtauxdT0 = DtauxdT0han(T0, T, p);
        end 
        % w.r.t. T
        if isempty(DtauxdThan)
          f   = @(x,p) gseg.taux(p(1), x, p(2:end));
          T0p = [T0; p];
          dtauxdT = coco_ezDFDX('f(x,p)', f, T, T0p);
        else
          dtauxdT = DtauxdThan(T0, T, p);
        end 
        % w.r.t. p
        if isempty(Dtauxdphan)
          f   = @(x,p) gseg.taux(p(1), p(2), x);
          T0T = [T0; T];
          dtauxdp = coco_ezDFDX('f(x,p)', f, p, T0T);
        else
          dtauxdp = Dtauxdphan(T0,T, p);
        end 
        
        
        if ~isempty(gseg.phi)
            tphi = tj-gseg.taut(T0,T,p);
            phij = gseg.phi(tphi');
            phidim  = size(phij,1);
            
            % Differentiation of tau_{j}^t 
            DtautdT0han = [];
            DtautdThan  = [];
            Dtautdphan  = [];
            switch numel(gseg.Dtaut)
                case 1
                    DtautdT0han = gseg.Dtaut{1};
                case 2
                    DtautdT0han = gseg.Dtaut{1};
                    DtautdThan  = gseg.Dtaut{2};
                case 3
                    DtautdT0han = gseg.Dtaut{1};
                    DtautdThan  = gseg.Dtaut{2};
                    Dtautdphan  = gseg.Dtaut{3};
            end
            % w.r.t. T0
            if isempty(DtautdT0han)
              f   = @(x,p) gseg.taut(x, p(1), p(2:end));
              Tp  = [T; p];
              dtautdT0 = coco_ezDFDX('f(x,p)', f, T0, Tp);
            else
              dtautdT0 = DtautdT0han(T0, T, p);
            end 
            % w.r.t. T
            if isempty(DtautdThan)
              f   = @(x,p) gseg.taut(p(1), x, p(2:end));
              T0p = [T0; p];
              dtautdT = coco_ezDFDX('f(x,p)', f, T, T0p);
            else
              dtautdT = DtautdThan(T0, T, p);
            end 
            % w.r.t. p
            if isempty(Dtautdphan)
              f   = @(x,p) gseg.taut(p(1), p(2), x);
              T0T = [T0; T];
              dtautdp = coco_ezDFDX('f(x,p)', f, p, T0T);
            else
              dtautdp = Dtautdphan(T0,T, p);
            end 
            
            % Differenation of phi
            if isempty(gseg.Dphi)
              f   = @(x,p) gseg.phi(x);
              Dphij = coco_ezDFDX('f(x,p)v', f, tphi', []);                
            else
               Dphij = gseg.Dphi(tphi');
            end
            
            % Differention of g
            dgdthan = [];
            dgdyhan = [];
            dgdxhan = [];
            dgdphihan = [];
            dgdphan = [];            
            switch numel(gseg.Dg)
                case 1
                    dgdthan = gseg.Dg{1};
                case 2
                    dgdthan = gseg.Dg{1};
                    dgdyhan = gseg.Dg{2};
                case 3
                    dgdthan = gseg.Dg{1};
                    dgdyhan = gseg.Dg{2};
                    dgdxhan = gseg.Dg{3};
                case 4
                    dgdthan = gseg.Dg{1};
                    dgdyhan = gseg.Dg{2};
                    dgdxhan = gseg.Dg{3};
                    dgdphihan = gseg.Dg{4};
                case 5
                    dgdthan = gseg.Dg{1};
                    dgdyhan = gseg.Dg{2};
                    dgdxhan = gseg.Dg{3};
                    dgdphihan = gseg.Dg{4};
                    dgdphan = gseg.Dg{5};
            end                
            
            % w.r.t t
            if isempty(dgdthan)
              f      = @(x,p) gseg.gfunc(x, p(1:yidim,:), p(yidim+1:yidim+dim,:), p(yidim+dim+1:yidim+dim+phidim,:), p(yidim+dim+phidim+1:end,:));
              yxphip = [yj; xj; phij; pj];
              dtode = coco_ezDFDX('f(x,p)v', f, tj', yxphip);
            else
              dtode = dgdthan(tj', yj, xj, phij, pj);
            end             
            % w.r.t. yi & yibp
            if isempty(dgdyhan)
              f      = @(x,p) gseg.gfunc(p(1,:), x, p(2:1+dim,:), p(2+dim:1+dim+phidim,:), p(2+dim+phidim:end,:));
              txphip = [tj'; xj; phij; pj];
              dyode = coco_ezDFDX('f(x,p)v', f, yj, txphip);
            else
              dyode = dgdyhan(tj', yj, xj, phij, pj);
            end
            dyrows = repmat(reshape(1:Nsa*yidim, [yidim Nsa]), [yidim 1]); 
            dycols = repmat(1:Nsa*yidim, [yidim 1]);
            dyode  = sparse(dyrows, dycols, dyode(:));
             % w.r.t x & xbp
            if isempty(dgdxhan)
              f      = @(x,p) gseg.gfunc(p(1,:), p(2:1+yidim,:), x, p(2+yidim:1+yidim+phidim,:), p(2+yidim+phidim:end,:));
              typhip = [tj'; yj; phij; pj];
              dxode = coco_ezDFDX('f(x,p)v', f, xj, typhip);
            else
              dxode = dgdxhan(tj', yj, xj, phij, pj);
            end
            dxrows = repmat(reshape(1:Nsa*yidim, [yidim Nsa]), [dim 1]); 
            dxcols = repmat(1:Nsa*dim, [yidim 1]);  
            dxode  = sparse(dxrows, dxcols, dxode(:));         
            % w.r.t. phi           
            if isempty(dgdphihan)
              f       = @(x,p) gseg.gfunc(p(1,:), p(2:1+yidim,:), p(2+yidim:1+yidim+dim,:), x, p(2+yidim+dim:end,:));
              tyxp    = [tj'; yj; xj; pj];
              dphiode = coco_ezDFDX('f(x,p)v', f, phij, tyxp);
            else
              dphiode = dgdphihan(tj', yj, xj, phij, pj);
            end  
            dphirows = repmat(reshape(1:Nsa*yidim, [yidim Nsa]), [phidim 1]); 
            dphicols = repmat(1:Nsa*phidim, [yidim 1]);  
            dphiode  = sparse(dphirows, dphicols, dphiode(:));             
            % w.r.t. p
            if isempty(dgdphan)
              f      = @(x,p) gseg.gfunc(p(1,:), p(2:1+yidim,:), p(2+yidim:1+yidim+dim,:), p(2+yidim+dim:end,:), x);
              tyxphi = [tj'; yj; xj; phij];
              dpode = coco_ezDFDX('f(x,p)v', f, pj, tyxphi);
            else
              dpode = dgdphan(tj', yj, xj, phij, pj);
            end
            dprows = repmat(reshape(1:Nsa*yidim, [yidim Nsa]), [pdim 1]); 
            dpcols = repmat(1:pdim, [yidim Nsa]);  
            dpode  = sparse(dprows, dpcols, dpode(:));             
            
            J_xbp = dxode*Lx*Tj;
            J_ybp = dyode*Pj*data.Waa;
            J_T0  = dtode(:)-dxode*DLx*x*dtauxdT0/T+dphiode*Dphij(:)*(1-dtautdT0);
            dtode = dtode.*repmat(sa', [yidim,1]);
            duode = Dphij.*repmat(sa'-dtautdT, [phidim,1]);
            J_T   = dtode(:)-dxode*DLx*x*(dtauxdT*T-gseg.taux(T0,T,p))/T^2+dphiode*duode(:);
            J_p   = dpode-dxode*DLx*x*dtauxdp/T-dphiode*Dphij(:)*dtautdp;
            
            J(s+1:s+Nsa*yidim,:) = [J_xbp J_ybp J_T0, J_T, J_p];
        else            
            % Differention of g
            dgdthan = [];
            dgdyhan = [];
            dgdxhan = [];
            dgdphan = [];            
            switch numel(gseg.Dg)
                case 1
                    dgdthan = gseg.Dg{1};
                case 2
                    dgdthan = gseg.Dg{1};
                    dgdyhan = gseg.Dg{2};
                case 3
                    dgdthan = gseg.Dg{1};
                    dgdyhan = gseg.Dg{2};
                    dgdxhan = gseg.Dg{3};
                case 4
                    dgdthan = gseg.Dg{1};
                    dgdyhan = gseg.Dg{2};
                    dgdxhan = gseg.Dg{3};
                    dgdphan = gseg.Dg{4};
            end                
            
            % w.r.t t
            if isempty(dgdthan)
              f     = @(x,p) gseg.gfunc(x, p(1:yidim,:), p(yidim+1:yidim+dim,:), p(yidim+dim+1:end,:));
              yxp   = [yj; xj; pj];
              dtode = coco_ezDFDX('f(x,p)v', f, tj', yxp);
            else
              dtode = dgdthan(tj', yj, xj, pj);
            end             
            % w.r.t. yi & yibp
            if isempty(dgdyhan)
              f     = @(x,p) gseg.gfunc(p(1,:), x, p(2:1+dim,:), p(2+dim:end,:));
              txp   = [tj'; xj; pj];
              dyode = coco_ezDFDX('f(x,p)v', f, yj, txp);
            else
              dyode = dgdyhan(tj', yj, xj, pj);
            end
            dyrows = repmat(reshape(1:Nsa*yidim, [yidim Nsa]), [yidim 1]); 
            dycols = repmat(1:Nsa*yidim, [yidim 1]);
            dyode  = sparse(dyrows, dycols, dyode(:));
             % w.r.t x & xbp
            if isempty(dgdxhan)
              f     = @(x,p) gseg.gfunc(p(1,:), p(2:1+yidim,:), x, p(2+yidim:end,:));
              typ   = [tj'; yj; pj];
              dxode = coco_ezDFDX('f(x,p)v', f, xj, typ);
            else
              dxode = dgdxhan(tj', yj, xj, pj);
            end
            dxrows = repmat(reshape(1:Nsa*yidim, [yidim Nsa]), [dim 1]); 
            dxcols = repmat(1:Nsa*dim, [yidim 1]);  
            dxode  = sparse(dxrows, dxcols, dxode(:));         
            % w.r.t. p
            if isempty(dgdphan)
              f     = @(x,p) gseg.gfunc(p(1,:), p(2:1+yidim,:), p(2+yidim:1+yidim+dim,:), x);
              tyx   = [tj'; yj; xj];
              dpode = coco_ezDFDX('f(x,p)v', f, pj, tyx);
            else
              dpode = dgdphan(tj', yj, xj, pj);
            end
            dprows = repmat(reshape(1:Nsa*yidim, [yidim Nsa]), [pdim 1]); 
            dpcols = repmat(1:pdim, [yidim Nsa]);  
            dpode  = sparse(dprows, dpcols, dpode(:));             
            
            J_xbp = dxode*Lx*Tj;
            J_ybp = dyode*Pj*data.Waa;
            J_T0  = dtode(:)-dxode*DLx*x*dtauxdT0/T;
            
            dt1   = repmat(sa', [yidim,1]);            
%             dtode = dtode.*repmat(sa', [yidim,1]);
            J_T   = dtode(:).*dt1(:)-dxode*DLx*x*(dtauxdT*T-gseg.taux(T0,T,p))/T^2;
            J_p   = dpode-dxode*DLx*x*dtauxdp/T;
            
            J(s+1:s+Nsa*yidim,:) = [J_xbp J_ybp J_T0, J_T, J_p];
            
%             f(s+1:s+Nsa*yidim) = gseg.gfunc(tj',yj,xj,pj);
        end 
    else
        if ~isempty(gseg.phi)         
            tphi = tj-gseg.taut(T0,T,p);
            phij = gseg.phi(tphi');
            phidim  = size(phij,1);
            
            % Differentiation of tau_{j}^t 
            DtautdT0han = [];
            DtautdThan  = [];
            Dtautdphan  = [];
            switch numel(gseg.Dtaut)
                case 1
                    DtautdT0han = gseg.Dtaut{1};
                case 2
                    DtautdT0han = gseg.Dtaut{1};
                    DtautdThan  = gseg.Dtaut{2};
                case 3
                    DtautdT0han = gseg.Dtaut{1};
                    DtautdThan  = gseg.Dtaut{2};
                    Dtautdphan  = gseg.Dtaut{3};
            end
            % w.r.t. T0
            if isempty(DtautdT0han)
              f   = @(x,p) gseg.taut(x, p(1), p(2:end));
              Tp  = [T; p];
              dtautdT0 = coco_ezDFDX('f(x,p)', f, T0, Tp);
            else
              dtautdT0 = DtautdT0han(T0, T, p);
            end 
            % w.r.t. T
            if isempty(DtautdThan)
              f   = @(x,p) gseg.taut(p(1), x, p(2:end));
              T0p = [T0; p];
              dtautdT = coco_ezDFDX('f(x,p)', f, T, T0p);
            else
              dtautdT = DtautdThan(T0, T, p);
            end 
            % w.r.t. p
            if isempty(Dtautdphan)
              f   = @(x,p) gseg.taut(p(1), p(2), x);
              T0T = [T0; T];
              dtautdp = coco_ezDFDX('f(x,p)', f, p, T0T);
            else
              dtautdp = Dtautdphan(T0,T, p);
            end 
            
            % Differenation of phi
            if isempty(gseg.Dphi)
              f   = @(x,p) gseg.phi(x);
              Dphij = coco_ezDFDX('f(x,p)v', f, tphi', tphi'); % may use 'f(o,x)v' option  
              Dphij = reshape(Dphij, [phidim,Nsa]);
            else
              Dphij = gseg.Dphi(tphi');
            end
            
            % Differention of g
            dgdthan = [];
            dgdyhan = [];
            dgdphihan = [];
            dgdphan = [];            
            switch numel(gseg.Dg)
                case 1
                    dgdthan = gseg.Dg{1};
                case 2
                    dgdthan = gseg.Dg{1};
                    dgdyhan = gseg.Dg{2};
                case 3
                    dgdthan = gseg.Dg{1};
                    dgdyhan = gseg.Dg{2};
                    dgdphihan = gseg.Dg{3};
                case 4
                    dgdthan = gseg.Dg{1};
                    dgdyhan = gseg.Dg{2};
                    dgdphihan = gseg.Dg{3};
                    dgdphan = gseg.Dg{4};
            end                
            
            % w.r.t t
            if isempty(dgdthan)
              f     = @(x,p) gseg.gfunc(x, p(1:yidim,:), p(yidim+1:yidim+phidim,:), p(yidim+phidim+1:end,:));
              yphip = [yj; phij; pj];
              dtode = coco_ezDFDX('f(x,p)v', f, tj', yphip);
            else
              dtode = dgdthan(tj', yj, phij, pj);
            end             
            % w.r.t. yi & yibp
            if isempty(dgdyhan)
              f     = @(x,p) gseg.gfunc(p(1,:), x, p(2:1+phidim,:), p(2+phidim:end,:));
              tphip = [tj'; phij; pj];
              dyode = coco_ezDFDX('f(x,p)v', f, yj, tphip);
            else
              dyode = dgdyhan(tj', yj, phij, pj);
            end
            dyrows = repmat(reshape(1:Nsa*yidim, [yidim Nsa]), [yidim 1]); 
            dycols = repmat(1:Nsa*yidim, [yidim 1]);
            dyode  = sparse(dyrows, dycols, dyode(:));     
            % w.r.t. phi           
            if isempty(dgdphihan)
              f       = @(x,p) gseg.gfunc(p(1,:), p(2:1+yidim,:), x, p(2+yidim:end,:));
              typ     = [tj'; yj; pj];
              dphiode = coco_ezDFDX('f(x,p)v', f, phij, typ);
            else
              dphiode = dgdphihan(tj', yj, phij, pj);
            end  
            dphirows = repmat(reshape(1:Nsa*yidim, [yidim Nsa]), [phidim 1]); 
            dphicols = repmat(1:Nsa*phidim, [yidim 1]);  
            dphiode  = sparse(dphirows, dphicols, dphiode(:));             
            % w.r.t. p
            if isempty(dgdphan)
              f     = @(x,p) gseg.gfunc(p(1,:), p(2:1+yidim,:), p(2+yidim:end,:), x);
              typhi = [tj'; yj; phij];
              dpode = coco_ezDFDX('f(x,p)v', f, pj, typhi);
            else
              dpode = dgdphan(tj', yj, phij, pj);
            end
            dprows = repmat(reshape(1:Nsa*yidim, [yidim Nsa]), [pdim 1]); 
            dpcols = repmat(1:pdim, [yidim Nsa]);  
            dpode  = sparse(dprows, dpcols, dpode(:));             
            
            J_xbp = sparse(Nsa*yidim, numel(x));
            J_ybp = dyode*Pj*data.Waa;
            J_T0  = dtode(:)+dphiode*Dphij(:)*(1-dtautdT0);
            dt1   = repmat(sa', [yidim,1]);
%             dtode = dtode.*repmat(sa, [yidim,1]);
            duode = Dphij.*repmat(sa'-dtautdT, [phidim,1]);
            J_T   = dtode(:).*dt1(:)+dphiode*duode(:);
            J_p   = dpode-dphiode*Dphij(:)*dtautdp;
            
            J(s+1:s+Nsa*yidim,:) = [J_xbp J_ybp J_T0, J_T, J_p];
                      
%            f(s+1:s+Nsa*yidim) = gseg.gfunc(tj',yj,phij,pj);
        else       
            % Differention of g
            dgdthan = [];
            dgdyhan = [];
            dgdphan = [];            
            switch numel(gseg.Dg)
                case 1
                    dgdthan = gseg.Dg{1};
                case 2
                    dgdthan = gseg.Dg{1};
                    dgdyhan = gseg.Dg{2};
                case 3
                    dgdthan = gseg.Dg{1};
                    dgdyhan = gseg.Dg{2};
                    dgdphan = gseg.Dg{3};
            end                
            
            % w.r.t t
            if isempty(dgdthan)
              f     = @(x,p) gseg.gfunc(x, p(1:yidim,:), p(yidim+1:end,:));
              yp    = [yj; pj];
              dtode = coco_ezDFDX('f(x,p)v', f, tj', yp);
            else
              dtode = dgdthan(tj', yj, pj);
            end             
            % w.r.t. yi & yibp
            if isempty(dgdyhan)
              f     = @(x,p) gseg.gfunc(p(1,:), x, p(2:end,:));
              tp    = [tj'; pj];
              dyode = coco_ezDFDX('f(x,p)v', f, yj, tp);
            else
              dyode = dgdyhan(tj', yj, pj);
            end
            dyrows = repmat(reshape(1:Nsa*yidim, [yidim Nsa]), [yidim 1]); 
            dycols = repmat(1:Nsa*yidim, [yidim 1]);
            dyode  = sparse(dyrows, dycols, dyode(:));                 
            % w.r.t. p
            if isempty(dgdphan)
              f     = @(x,p) gseg.gfunc(p(1,:), p(2:end,:), x);
              ty    = [tj'; yj];
              dpode = coco_ezDFDX('f(x,p)v', f, pj, ty);
            else
              dpode = dgdphan(tj', yj, pj);
            end
            dprows = repmat(reshape(1:Nsa*yidim, [yidim Nsa]), [pdim 1]); 
            dpcols = repmat(1:pdim, [yidim Nsa]);  
            dpode  = sparse(dprows, dpcols, dpode(:));             
            
            J_xbp = sparse(Nsa*yidim, numel(x));
            J_ybp = dyode*Pj*data.Waa;
            J_T0  = dtode(:);
            dtode = dtode.*repmat(sa', [yidim,1]);
            J_T   = dtode(:);
            J_p   = dpode;
            
            J(s+1:s+Nsa*yidim,:) = [J_xbp J_ybp J_T0, J_T, J_p];
            
%             f(s+1:s+Nsa*yidim) = gseg.gfunc(tj',yj,pj);
        end
    end
    s = s+Nsa*yidim;
%     toc
end

% toc

% tic
% [data, DadF] = coco_ezDFDX('f(o,d,x)',  prob, data, @alg_dde_F, u);
% % toc
% err = max(max(abs(J-DadF)));
% if err>1e-4
%     disp('Analytic Jacobian is incorrect for some entries');
% end
% J = DadF;
end



function A = coll_L(ts, tz)
%COLL_L   Evaluation of Lagrange polynomials.
%
% Use high-dimensional arrays for vectorized evaluation.
%
% A = COLL_L(TS, TZ)
% 
% A  - Array of interpolated values.
% TS - Array of basepoints.
% TZ - Array of interpolation points.

q = numel(ts);
p = numel(tz);

zi = repmat(reshape(tz, [p 1 1]), [1 q q]);
sj = repmat(reshape(ts, [1 q 1]), [p 1 q]);
sk = repmat(reshape(ts, [1 1 q]), [p q 1]);

t1 = zi-sk;
t2 = sj-sk;
idx = find(abs(t2)<=eps);
t1(idx) = 1;
t2(idx) = 1;

A = prod(t1./t2, 3);

end


function A = coll_Lp(ts, tz)
%COLL_LP   Evaluation of derivative of Lagrange polynomials.
%
% Use high-dimensional arrays for vectorized evaluation.
%
% A = COLL_LP(TS, TZ)
% 
% A  - Array of interpolated values
% TS - Array of basepoints
% TZ - Array of interpolation points

q = numel(ts);
p = numel(tz);

zi = repmat(reshape(tz, [p 1 1 1]), [1 q q q]);
sj = repmat(reshape(ts, [1 q 1 1]), [p 1 q q]);
sk = repmat(reshape(ts, [1 1 q 1]), [p q 1 q]);
sl = repmat(reshape(ts, [1 1 1 q]), [p q q 1]);

t3 = sj(:,:,:,1)-sk(:,:,:,1);
t4 = zi-sl;
t5 = sj-sl;

idx1 = find(abs(t5)<=eps);
idx2 = find(abs(t3)<=eps);
idx3 = find(abs(sk-sl)<=eps);
t5(union(idx1, idx3)) = 1;
t4(union(idx1, idx3)) = 1;
t3(idx2) = 1;
t3       = 1.0./t3;
t3(idx2) = 0;

A = sum(t3.*prod(t4./t5, 4), 3);

end
