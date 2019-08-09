function prob = alg_ddae_construct_seg(prob, tbid, data)
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
    gseg.tauy  = [];
    gseg.taut  = [];
    gseg.phi   = [];
    gseg.Dg    = [];
    gseg.Dtaux = [];
    gseg.Dtauy = [];
    gseg.Dtaut = [];
    gseg.Dphi  = [];
    ngarg      = numel(gj);
    if ngarg>=2
        gseg.taux = gj{2};
        if ngarg>=3
            gseg.tauy = gj{3};
            if ngarg>=4
                gseg.taut = gj{4};
                if ngarg>=5
                    gseg.phi = gj{5};
                    if ngarg>=6
                        gseg.Dg = gj{6};
                        if ngarg>=7
                            gseg.Dtaux = gj{7};
                            if ngarg>=8
                                gseg.Dtauy = gj{8};
                                if ngarg>=9
                                    gseg.Dtaut = gj{9};
                                    if ngarg>=10
                                        gseg.Dphi = gj{10};
                                    end
                                end
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
prob = coco_add_func(prob, fid, @alg_ddae_F, @alg_ddae_DFDU, data, 'zero', 'uidx', data.ddae_uidx);
prob = coco_add_slot(prob, fid, @coco_save_data, data, 'save_full');

end

function [data f] = alg_ddae_F(prob, data, u)
%COLL_F   Collocation zero problem.
%
% Expects vectorized encoding of vector field as function of two arguments.

seg = data.ddaecoll_seg;
x   = u(seg.xbp_idx); % Extract basepoint values
y   = u(seg.ybp_idx);
T0  = u(seg.T0_idx);
T   = u(seg.T_idx);   % Extract interval length
p   = u(seg.p_idx);   % Extract problem parameters

xx  = reshape(x, seg.xbp_shp); % Values at base nodes
% pp  = repmat(p, seg.p_rep);
yy  = reshape(y, seg.y_shp); 
% tcn = seg.Wta*(seg.Ma*T+T0);
xcn   = reshape(seg.Wad*x, seg.x_shp);
ycn   = reshape(seg.Waa*y, seg.y_shp);

NTST  = seg.ddaecoll.NTST;
NCOL  = seg.ddaecoll.NCOL;
tsa   = seg.tsa;
tsd   = seg.tsd;
gdim  = data.gdim;
ydim  = seg.ydim;
dim   = seg.dim;
f     = zeros(data.gcndim,1);
s     = 0;
yref  = zeros(data.gcndim,1);
sref  = 0;
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
    
    tj   = T0+T*sa;                                % time at collocation nodes sa
    xj   = xcn(:,idx);                             % xcn at collocation nodes sa
    yj   = ycn(:,idx);                             % ycn at collocation nodes sa
    pj   = repmat(p, [1,Nsa]);  
    
    %% For-loop based implementation
    flag = false;
    if flag
    if ~isempty(gseg.taux)
        % shifted node sets and the sets of intervals
        sax   = sa-gseg.taux(T0,T,p)/T;
        idsax = floor(seg.ddaecoll.NTST*sax)+1;    % index for subinterval belongs to
        idsax(idsax==NTST+1) = NTST;
    end 
    if ~isempty(gseg.tauy)
        % shifted node sets and the sets of intervals
        say   = sa-gseg.tauy(T0,T,p)/T;
        idsay = floor(seg.ddaecoll.NTST*say)+1;    % index for subinterval belongs to
        idsay(idsay==NTST+1) = NTST;
    end     
    for k=1:Nsa
        Tk = tj(k);                               % go to [T0,T0+T]
        kd = (k-1)*gdim+1:k*gdim;                 % collocation index
        if ~isempty(gseg.taux)
            tx = 2*(NTST*sax(k)-(idsax(k)-1))-1;                  % go to [-1,1]
            xk = xx(:,(idsax(k)-1)*(NCOL+1)+1:idsax(k)*(NCOL+1)); % base pts in that subinterval
            xk = coll_L(tsd, tx)*xk';                             % interpolation
            xk = xk';
            if ~isempty(gseg.tauy)
                ty = 2*(NTST*say(k)-(idsay(k)-1))-1;                  % go to [-1,1]
                yk = yy(:,(idsay(k)-1)*NCOL+1:idsay(k)*NCOL);         % base pts in that subinterval
                yk = coll_L(tsa, ty)*yk';                             % interpolation
                yk = yk';
                if ~isempty(gseg.phi)
                    tphi = Tk-gseg.taut(T0,T,p);
                    phik = gseg.phi(tphi);
                    yref(s+kd) = gseg.gfunc(Tk,xj(:,k),yj(:,k),xk,yk,phik,p);
                else
                    yref(s+kd) = gseg.gfunc(Tk,xj(:,k),yj(:,k),xk,yk,p);
                end
            else
                if ~isempty(gseg.phi)
                    tphi = Tk-gseg.taut(T0,T,p);
                    phik = gseg.phi(tphi);
                    yref(s+kd) = gseg.gfunc(Tk,xj(:,k),yj(:,k),xk,phik,p);
                else
                    yref(s+kd) = gseg.gfunc(Tk,xj(:,k),yj(:,k),xk,p);
                end                
            end
        else
            if ~isempty(gseg.tauy)
                ty = 2*(NTST*say(k)-(idsay(k)-1))-1;                  % go to [-1,1]
                yk = yy(:,(idsay(k)-1)*NCOL+1:idsay(k)*NCOL);         % base pts in that subinterval
                yk = coll_L(tsa, ty)*yk';                             % interpolation
                yk = yk';
                if ~isempty(gseg.phi)
                    tphi = Tk-gseg.taut(T0,T,p);
                    phik = gseg.phi(tphi);
                    yref(s+kd) = gseg.gfunc(Tk,xj(:,k),yj(:,k),yk,phik,p);
                else
                    yref(s+kd) = gseg.gfunc(Tk,xj(:,k),yj(:,k),yk,p);
                end
            else
                if ~isempty(gseg.phi)
                    tphi = Tk-gseg.taut(T0,T,p);
                    phik = gseg.phi(tphi);
                    yref(s+kd) = gseg.gfunc(Tk,xj(:,k),yj(:,k),phik,p);
                else
                    yref(s+kd) = gseg.gfunc(Tk,xj(:,k),yj(:,k),p);
                end                
            end            
        end
    end
    end
%     s   = s+Nsa*gdim; 
    
    %% vectorization based implementation 
    rowp = (1:dim*Nsa);
    colp = rowp+(idx(1)-1)*dim;
    Pxj  = sparse(rowp, colp, ones(dim*Nsa,1), dim*Nsa, dim*NTST*NCOL);
    xj   = Pxj*seg.Wad*x;
    xj   = reshape(xj, [dim,Nsa]);
    
    rowp = (1:ydim*Nsa);
    colp = rowp+(idx(1)-1)*ydim;
    Pyj  = sparse(rowp, colp, ones(ydim*Nsa,1), ydim*Nsa, ydim*NTST*NCOL);
    yj   = Pyj*seg.Waa*y;
    yj   = reshape(yj, [ydim,Nsa]); 
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
            tx  = 2*(NTST*sax-(jund-1))-1;     % go to [-1,1]
            Lx  = kron(coll_L(tsd, tx),eye(dim));
            xsj = Lx*xbpj;
            xsj = reshape(xsj, [dim, Nsa]);
%             xsj = coll_L(tsd, tx)*xbpj';
%             xsj = xsj';
        elseif jund==jbar-2              % Two subintervals
            id1 = numel(find(idsax==idsax(1)));
            tx1 = 2*(NTST*sax(1:id1)-(jund-1))-1;
            Lx1 = coll_L(tsd, tx1);           % 1st subinterval
            tx2 = 2*(NTST*sax(id1+1:end)-jund)-1;
            Lx2 = coll_L(tsd, tx2);           % 2nd subinterval
            Lx  = blkdiag(kron(Lx1,eye(dim)), kron(Lx2,eye(dim)));
            xsj = Lx*xbpj;
            xsj = reshape(xsj, [dim, Nsa]);            
%             Lx  = blkdiag(Lx1, Lx2);
%             xsj  = Lx*xbpj';
%             xsj  = xsj';
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
            xsj = Lx*xbpj;
            xsj = reshape(xsj, [dim, Nsa]);              
%             Lx  = blkdiag(LxF, kron(eye(jbar-jund-2), Lxs), LxL);
%             xsj  = Lx*xbpj';
%             xsj  = xsj';
        else
            disp('bar{j} and underline{j} is not matched well');
        end
        
        if ~isempty(gseg.tauy)
            % shifted node sets and the sets of intervals
            say   = sa-gseg.tauy(T0,T,p)/T;
            idsay = floor(seg.ddaecoll.NTST*say)+1;    % index for subinterval belongs to
            idsay(idsay==NTST+1) = NTST;
            numid = numel(unique(idsay));
            dimb = numid*ydim*NCOL;
            rowt = (1:dimb);
            colt = rowt+(idsay(1)-1)*ydim*NCOL;
            Tj   = sparse(rowt, colt, ones(dimb,1), dimb, data.ycndim);
            ybpj = Tj*y;             
%             ybpj = yy(:,(idsay(1)-1)*NCOL+1:idsay(end)*NCOL);
            jund = find(data.mg<=say(1), 1, 'last');
            jbar = find(data.mg>say(end), 1, 'first');
            if isempty(jbar) % at the right end of whole interval, i.e., 1
                jbar = NTST+1;
            end
            if jund==jbar-1                  % A single subinterval
                ty  = 2*(NTST*say-(jund-1))-1;     % go to [-1,1]
                Ly  = kron(coll_L(tsa, ty),eye(ydim));
                ysj = Ly*ybpj;
                ysj = reshape(ysj, [ydim, Nsa]);             
%                 ysj = coll_L(tsa, ty)*ybpj';
%                 ysj = ysj';
            elseif jund==jbar-2              % Two subintervals
                id1 = numel(find(idsay==idsay(1)));
                ty1 = 2*(NTST*say(1:id1)-(jund-1))-1;
                Ly1 = coll_L(tsa, ty1);           % 1st subinterval
                ty2 = 2*(NTST*say(id1+1:end)-jund)-1;
                Ly2 = coll_L(tsa, ty2);           % 2nd subinterval
                Ly  = blkdiag(kron(Ly1,eye(ydim)), kron(Ly2,eye(ydim)));
                ysj = Ly*ybpj;
                ysj = reshape(ysj, [ydim, Nsa]);                  
%                 Ly  = blkdiag(Ly1, Ly2);
%                 ysj  = Ly*ybpj';
%                 ysj  = ysj';
            elseif jund<jbar-2               % More than two subintervals
                % We use periodicity of collocation nodes to simplify
                % calculation
                nyF = numel(find(idsay==idsay(1)));
                nyL = numel(find(idsay==idsay(end)));
                tss = 2*(NTST*say(nyF+1:nyF+NCOL)-jund)-1;
                Lys = coll_L(tsa, tss);           % 2nd - (Last-1) subintervals            
                LyF = Lys(NCOL-nyF+1:end,:);      % First subintervals
                if nyL>NCOL
                    tse = 2*(NTST*say(end-nyL+1:end)-(jbar-2))-1;
    %                 LxL = [Lxs; Lextra];          % Last subintervals
                    LyL = coll_L(tsa, tse);       % Last subintervals
                else
                    LyL = Lys(1:nyL,:);           % Last subintervals
                end  
                Ly1 = kron(LyF,eye(ydim));
                Ly2 = kron(eye(jbar-jund-2), kron(Lys,eye(ydim)));
                Ly3 = kron(LyL,eye(ydim));
                Ly  = blkdiag(Ly1,Ly2,Ly3);
                ysj = Ly*ybpj;
                ysj = reshape(ysj, [ydim, Nsa]);                  
%                 Ly   = blkdiag(LyF, kron(eye(jbar-jund-2), Lys), LyL);
%                 ysj  = Ly*ybpj';
%                 ysj  = ysj';
            else
                disp('bar{j} and underline{j} is not matched well');
            end
        
            if ~isempty(gseg.phi)
                tphi = tj-gseg.taut(T0,T,p);
                phij = gseg.phi(tphi);
                f(s+1:s+Nsa*gdim) = gseg.gfunc(tj',xj,yj,xsj,ysj,phij,pj);
            else
                f(s+1:s+Nsa*gdim) = gseg.gfunc(tj',xj,yj,xsj,ysj,pj);
            end 
        else
            if ~isempty(gseg.phi)
                tphi = tj-gseg.taut(T0,T,p);
                phij = gseg.phi(tphi);
                f(s+1:s+Nsa*gdim) = gseg.gfunc(tj',xj,yj,xsj,phij,pj);
            else
                f(s+1:s+Nsa*gdim) = gseg.gfunc(tj',xj,yj,xsj,pj);
            end 
        end
    else
        if ~isempty(gseg.tauy)
            % shifted node sets and the sets of intervals
            say   = sa-gseg.tauy(T0,T,p)/T;
            idsay = floor(seg.ddaecoll.NTST*say)+1;    % index for subinterval belongs to
            idsay(idsay==NTST+1) = NTST;
            numid = numel(unique(idsay));
            dimb = numid*ydim*NCOL;
            rowt = (1:dimb);
            colt = rowt+(idsay(1)-1)*ydim*NCOL;
            Tj   = sparse(rowt, colt, ones(dimb,1), dimb, data.ycndim);
            ybpj = Tj*y;             
%             ybpj = yy(:,(idsay(1)-1)*NCOL+1:idsay(end)*NCOL);
            jund = find(data.mg<=say(1), 1, 'last');
            jbar = find(data.mg>say(end), 1, 'first');
            if isempty(jbar) % at the right end of whole interval, i.e., 1
                jbar = NTST+1;
            end
            if jund==jbar-1                  % A single subinterval
                ty  = 2*(NTST*say-(jund-1))-1;     % go to [-1,1]
                Ly  = kron(coll_L(tsa, ty),eye(ydim));
                ysj = Ly*ybpj;
                ysj = reshape(ysj, [ydim, Nsa]);             
%                 ysj = coll_L(tsa, ty)*ybpj';
%                 ysj = ysj';
            elseif jund==jbar-2              % Two subintervals
                id1 = numel(find(idsay==idsay(1)));
                ty1 = 2*(NTST*say(1:id1)-(jund-1))-1;
                Ly1 = coll_L(tsa, ty1);           % 1st subinterval
                ty2 = 2*(NTST*say(id1+1:end)-jund)-1;
                Ly2 = coll_L(tsa, ty2);           % 2nd subinterval
                Ly  = blkdiag(kron(Ly1,eye(ydim)), kron(Ly2,eye(ydim)));
                ysj = Ly*ybpj;
                ysj = reshape(ysj, [ydim, Nsa]);                  
%                 Ly  = blkdiag(Ly1, Ly2);
%                 ysj  = Ly*ybpj';
%                 ysj  = ysj';
            elseif jund<jbar-2               % More than two subintervals
                % We use periodicity of collocation nodes to simplify
                % calculation
                nyF = numel(find(idsay==idsay(1)));
                nyL = numel(find(idsay==idsay(end)));
                tss = 2*(NTST*say(nyF+1:nyF+NCOL)-jund)-1;
                Lys = coll_L(tsa, tss);           % 2nd - (Last-1) subintervals            
                LyF = Lys(NCOL-nyF+1:end,:);      % First subintervals
                if nyL>NCOL
                    tse = 2*(NTST*say(end-nyL+1:end)-(jbar-2))-1;
    %                 LxL = [Lxs; Lextra];          % Last subintervals
                    LyL = coll_L(tsa, tse);       % Last subintervals
                else
                    LyL = Lys(1:nyL,:);           % Last subintervals
                end  
                Ly1 = kron(LyF,eye(ydim));
                Ly2 = kron(eye(jbar-jund-2), kron(Lys,eye(ydim)));
                Ly3 = kron(LyL,eye(ydim));
                Ly  = blkdiag(Ly1,Ly2,Ly3);
                ysj = Ly*ybpj;
                ysj = reshape(ysj, [ydim, Nsa]);                  
%                 Ly   = blkdiag(LyF, kron(eye(jbar-jund-2), Lys), LyL);
%                 ysj  = Ly*ybpj';
%                 ysj  = ysj';
            else
                disp('bar{j} and underline{j} is not matched well');
            end        
            if ~isempty(gseg.phi)
                tphi = tj-gseg.taut(T0,T,p);
                phij = gseg.phi(tphi);
                f(s+1:s+Nsa*gdim) = gseg.gfunc(tj',xj,yj,ysj,phij,pj);
            else
                f(s+1:s+Nsa*gdim) = gseg.gfunc(tj',xj,yj,ysj,pj);
            end 
        else
            if ~isempty(gseg.phi)
                tphi = tj-gseg.taut(T0,T,p);
                phij = gseg.phi(tphi);
                f(s+1:s+Nsa*gdim) = gseg.gfunc(tj',xj,yj,phij,pj);
            else
                f(s+1:s+Nsa*gdim) = gseg.gfunc(tj',xj,yj,pj);
            end 
        end
    end
    s = s+Nsa*gdim;
%     toc
end

% assert(max(abs(f-yref))<1e-8, 'the vectorization is incorrect');


assert(s==data.gcndim, 'Number of collocation equations is incorret in %s', data.gidx);

end

function [data J] = alg_ddae_DFDU(prob, data, u)
%COLL_DFDU   Linearization of collocation zero problem.
%
% Expects vectorized encoding of Jacobians of vector field with respect to
% its arguments. When dfdxhan and/or dfdphan are empty, approximate
% Jacobians are obtained using numerical differentiation.
seg = data.ddaecoll_seg;
x   = u(seg.xbp_idx); % Extract basepoint values
y   = u(seg.ybp_idx);
T0  = u(seg.T0_idx);
T   = u(seg.T_idx);   % Extract interval length
p   = u(seg.p_idx);   % Extract problem parameters

% xx  = reshape(x, seg.xbp_shp); % Values at base nodes
% pp  = repmat(p, seg.p_rep);
% yy  = reshape(y, data.y_shp); 
% tcn = seg.Wta*(seg.Ma*T+T0);

% ycn   = reshape(data.Waa*y, data.y_shp);

NTST  = seg.ddaecoll.NTST;
NCOL  = seg.ddaecoll.NCOL;
tsa   = seg.tsa;
tsd   = seg.tsd;
dim   = seg.dim;
pdim  = seg.pdim;
ydim  = seg.ydim;
gdim  = data.gdim;
J     = zeros(data.gcndim,numel(u));
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
    
    
    rowp = (1:dim*Nsa);
    colp = rowp+(idx(1)-1)*dim;
    Pxj  = sparse(rowp, colp, ones(dim*Nsa,1), dim*Nsa, dim*NTST*NCOL);
    xj   = Pxj*seg.Wad*x;
    xj   = reshape(xj, [dim,Nsa]);                 % xcn at collocation nodes sa
    
    rowp = (1:ydim*Nsa);
    colp = rowp+(idx(1)-1)*ydim;
    Pyj  = sparse(rowp, colp, ones(ydim*Nsa,1), ydim*Nsa, ydim*NTST*NCOL);
    yj   = Pyj*seg.Waa*y;             
    yj   = reshape(yj, [ydim,Nsa]);                % ycn at collocation nodes sa
    
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
        else
            disp('bar{j} and underline{j} is not matched well');
        end
        xsj = Lx*xbpj;
        xsj = reshape(xsj, [dim, Nsa]);
        Lx  = Lx*Tj;
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
    end
        
        
    if ~isempty(gseg.tauy)
        % shifted node sets and the sets of intervals
        say   = sa-gseg.tauy(T0,T,p)/T;
        idsay = floor(seg.ddaecoll.NTST*say)+1;    % index for subinterval belongs to
        idsay(idsay==NTST+1) = NTST;
        numid = numel(unique(idsay));
        dimb = numid*ydim*NCOL;
        rowt = (1:dimb);
        colt = rowt+(idsay(1)-1)*ydim*NCOL;
        Tj   = sparse(rowt, colt, ones(dimb,1), dimb, data.ycndim);
        ybpj = Tj*y;             
        jund = find(data.mg<=say(1), 1, 'last');
        jbar = find(data.mg>say(end), 1, 'first');
        if isempty(jbar) % at the right end of whole interval, i.e., 1
            jbar = NTST+1;
        end
        if jund==jbar-1                  % A single subinterval
            ty  = 2*(NTST*say-(jund-1))-1;     % go to [-1,1]
            Ly  = kron(coll_L(tsa, ty),eye(ydim));
            DLy = kron(coll_Lp(tsa, ty),eye(ydim)); 
        elseif jund==jbar-2              % Two subintervals
            id1 = numel(find(idsay==idsay(1)));
            ty1 = 2*(NTST*say(1:id1)-(jund-1))-1;
            Ly1 = coll_L(tsa, ty1);           % 1st subinterval
            DLy1 = coll_Lp(tsa, ty1);
            ty2 = 2*(NTST*say(id1+1:end)-jund)-1;
            Ly2 = coll_L(tsa, ty2);           % 2nd subinterval
            DLy2 = coll_Lp(tsa, ty2);
            Ly  = blkdiag(kron(Ly1,eye(ydim)), kron(Ly2,eye(ydim)));
            DLy  = blkdiag(kron(DLy1,eye(ydim)), kron(DLy2,eye(ydim)));
        elseif jund<jbar-2               % More than two subintervals
            % We use periodicity of collocation nodes to simplify
            % calculation
            nyF = numel(find(idsay==idsay(1)));
            nyL = numel(find(idsay==idsay(end)));
            tss = 2*(NTST*say(nyF+1:nyF+NCOL)-jund)-1;
            Lys = coll_L(tsa, tss);           % 2nd - (Last-1) subintervals            
            LyF = Lys(NCOL-nyF+1:end,:);      % First subintervals
            DLys = coll_Lp(tsa, tss);
            DLyF = DLys(NCOL-nyF+1:end,:);
            if nyL>NCOL
                tse = 2*(NTST*say(end-nyL+1:end)-(jbar-2))-1;
                LyL = coll_L(tsa, tse);       % Last subintervals
                DLyL = coll_Lp(tsa, tse);
            else
                LyL = Lys(1:nyL,:);           % Last subintervals
                DLyL = DLys(1:nyL,:);
            end  
            Ly1 = kron(LyF,eye(ydim));
            Ly2 = kron(eye(jbar-jund-2), kron(Lys,eye(ydim)));
            Ly3 = kron(LyL,eye(ydim));
            Ly  = blkdiag(Ly1,Ly2,Ly3);
            DLy1 = kron(DLyF,eye(ydim));
            DLy2 = kron(eye(jbar-jund-2), kron(DLys,eye(ydim)));
            DLy3 = kron(DLyL,eye(ydim));
            DLy  = blkdiag(DLy1,DLy2,DLy3);                
        else
            disp('bar{j} and underline{j} is not matched well');
        end
        ysj = Ly*ybpj;
        ysj = reshape(ysj, [ydim, Nsa]);  
        Ly  = Ly*Tj;
        DLy = 2*NTST*DLy*Tj;
        
         % Differentiation of tau_{j}^y 
        DtauydT0han = [];
        DtauydThan  = [];
        Dtauydphan  = [];
        switch numel(gseg.Dtauy)
            case 1
                DtauydT0han = gseg.Dtauy{1};
            case 2
                DtauydT0han = gseg.Dtauy{1};
                DtauydThan  = gseg.Dtauy{2};
            case 3
                DtauydT0han = gseg.Dtauy{1};
                DtauydThan  = gseg.Dtauy{2};
                Dtauydphan  = gseg.Dtauy{3};
        end
        % w.r.t. T0
        if isempty(DtauydT0han)
          f   = @(x,p) gseg.tauy(x, p(1), p(2:end));
          Tp  = [T; p];
          dtauydT0 = coco_ezDFDX('f(x,p)', f, T0, Tp);
        else
          dtauydT0 = DtauydT0han(T0, T, p);
        end 
        % w.r.t. T
        if isempty(DtauydThan)
          f   = @(x,p) gseg.tauy(p(1), x, p(2:end));
          T0p = [T0; p];
          dtauydT = coco_ezDFDX('f(x,p)', f, T, T0p);
        else
          dtauydT = DtauydThan(T0, T, p);
        end 
        % w.r.t. p
        if isempty(Dtauydphan)
          f   = @(x,p) gseg.tauy(p(1), p(2), x);
          T0T = [T0; T];
          dtauydp = coco_ezDFDX('f(x,p)', f, p, T0T);
        else
          dtauydp = Dtauydphan(T0,T, p);
        end       
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
%           Dphij = coco_ezDFDX('f(x,p)v', f, tphi', []); 
          Dphij = coco_ezDFDX('f(x,p)v', f, tphi', tphi'); % may use 'f(o,x)v' option  
          Dphij = reshape(Dphij, [phidim,Nsa]);         
        else
          Dphij = gseg.Dphi(tphi');
        end
    end

    
    if ~isempty(gseg.taux)
        if ~isempty(gseg.tauy)
            if ~isempty(gseg.phi)
                % Differention of g
                dgdthan  = [];
                dgdxhan  = [];
                dgdyhan  = [];
                dgdxshan = [];
                dgdyshan = [];
                dgdphihan = [];
                dgdphan  = [];            
                switch numel(gseg.Dg)
                    case 1
                        dgdthan = gseg.Dg{1};
                    case 2
                        dgdthan = gseg.Dg{1};
                        dgdxhan = gseg.Dg{2};
                    case 3
                        dgdthan = gseg.Dg{1};
                        dgdxhan = gseg.Dg{2};
                        dgdyhan = gseg.Dg{3};
                    case 4
                        dgdthan  = gseg.Dg{1};
                        dgdxhan  = gseg.Dg{2};
                        dgdyhan  = gseg.Dg{3};
                        dgdxshan = gseg.Dg{4};
                    case 5
                        dgdthan  = gseg.Dg{1};
                        dgdxhan  = gseg.Dg{2};
                        dgdyhan  = gseg.Dg{3};
                        dgdxshan = gseg.Dg{4};
                        dgdyshan = gseg.Dg{5};
                    case 6
                        dgdthan  = gseg.Dg{1};
                        dgdxhan  = gseg.Dg{2};
                        dgdyhan  = gseg.Dg{3};
                        dgdxshan = gseg.Dg{4};
                        dgdyshan = gseg.Dg{5};
                        dgdphihan = gseg.Dg{6};
                    case 7
                        dgdthan  = gseg.Dg{1};
                        dgdxhan  = gseg.Dg{2};
                        dgdyhan  = gseg.Dg{3};
                        dgdxshan = gseg.Dg{4};
                        dgdyshan = gseg.Dg{5};
                        dgdphihan = gseg.Dg{6};
                        dgdphan  = gseg.Dg{7};
                end
                
                % w.r.t t
                if isempty(dgdthan)
                  f          = @(x,p) gseg.gfunc(x, p(1:dim,:), p(dim+1:ydim+dim,:), p(ydim+dim+1:ydim+2*dim,:), p(ydim+2*dim+1:2*ydim+2*dim,:), p(2*ydim+2*dim+1:2*ydim+2*dim+phidim,:), p(2*ydim+2*dim+phidim+1:end,:));
                  xyxsysphip = [xj; yj; xsj; ysj; phij; pj];
                  dtode      = coco_ezDFDX('f(x,p)v', f, tj', xyxsysphip);
                  dtode      = reshape(dtode,[gdim,Nsa]);
                else
                  dtode      = dgdthan(tj', xj, yj, xsj, ysj, phij, pj);
                end
                % w.r.t x
                if isempty(dgdxhan)
                  f          = @(x,p) gseg.gfunc(p(1,:), x, p(2:ydim+1,:), p(ydim+2:ydim+dim+1,:), p(ydim+dim+2:2*ydim+dim+1,:), p(2*ydim+dim+2:2*ydim+dim+phidim+1,:), p(2*ydim+dim+phidim+2:end,:));
                  tyxsysphip = [tj'; yj; xsj; ysj; phij; pj];
                  dxode      = coco_ezDFDX('f(x,p)v', f, xj, tyxsysphip);
                else
                  dxode      = dgdxhan(tj', xj, yj, xsj, ysj, phij, pj);
                end
                dxrows = repmat(reshape(1:Nsa*gdim, [gdim Nsa]), [dim 1]); 
                dxcols = repmat(1:Nsa*dim, [gdim 1]);  
                dxode  = sparse(dxrows, dxcols, dxode(:));                
                % w.r.t. y
                if isempty(dgdyhan)
                  f          = @(x,p) gseg.gfunc(p(1,:), p(2:dim+1,:), x, p(dim+2:2*dim+1,:), p(2*dim+2:2*dim+ydim+1,:), p(2*dim+ydim+2:2*dim+ydim+phidim+1,:), p(2*dim+ydim+phidim+2:end,:));
                  txxsysphip = [tj'; xj; xsj; ysj; phij; pj];
                  dyode      = coco_ezDFDX('f(x,p)v', f, yj, txxsysphip);
                else
                  dyode      = dgdyhan(tj', xj, yj, xsj, ysj, phij, pj);
                end
                dyrows = repmat(reshape(1:Nsa*gdim, [gdim Nsa]), [ydim 1]); 
                dycols = repmat(1:Nsa*ydim, [gdim 1]);  
                dyode  = sparse(dyrows, dycols, dyode(:));
                % w.r.t xs
                if isempty(dgdxshan)
                  f          = @(x,p) gseg.gfunc(p(1,:), p(2:dim+1,:), p(dim+2:ydim+dim+1,:), x, p(ydim+dim+2:2*ydim+dim+1,:), p(2*ydim+dim+2:2*ydim+dim+phidim+1,:), p(2*ydim+dim+phidim+2:end,:));
                  txyysphip  = [tj'; xj; yj; ysj; phij; pj];
                  dxsode     = coco_ezDFDX('f(x,p)v', f, xsj, txyysphip);
                else
                  dxsode     = dgdxshan(tj', xj, yj, xsj, ysj, phij, pj);
                end
                dxsode   = sparse(dxrows, dxcols, dxsode(:)); 
                % w.r.t. ys
                if isempty(dgdyshan)
                  f          = @(x,p) gseg.gfunc(p(1,:), p(2:dim+1,:), p(dim+2:dim+ydim+1,:), p(dim+ydim+2:2*dim+ydim+1,:), x, p(2*dim+ydim+2:2*dim+ydim+phidim+1,:), p(2*dim+ydim+phidim+2:end,:));
                  txyxsphip  = [tj'; xj; yj; xsj; phij; pj];
                  dysode     = coco_ezDFDX('f(x,p)v', f, ysj, txyxsphip);
                else
                  dysode     = dgdyshan(tj', xj, yj, xsj, ysj, phij, pj);
                end                
                dysode   = sparse(dyrows, dycols, dysode(:)); 
                % w.r.t. phi           
                if isempty(dgdphihan)
                  f          = @(x,p) gseg.gfunc(p(1,:), p(2:dim+1,:), p(dim+2:dim+ydim+1,:), p(dim+ydim+2:2*dim+ydim+1,:), p(2*dim+ydim+2:2*dim+2*ydim+1,:), x, p(2*dim+2*ydim+2:end,:));
                  txyxsysp   = [tj'; xj; yj; xsj; ysj; pj];
                  dphiode    = coco_ezDFDX('f(x,p)v', f, phij, txyxsysp);
                else
                  dphiode    = dgdphihan(tj', xj, yj, xsj, ysj, phij, pj);
                end  
                dphirows = repmat(reshape(1:Nsa*gdim, [gdim Nsa]), [phidim 1]); 
                dphicols = repmat(1:Nsa*phidim, [gdim 1]);  
                dphiode  = sparse(dphirows, dphicols, dphiode(:));             
                % w.r.t. p
                if isempty(dgdphan)
                  f          = @(x,p) gseg.gfunc(p(1,:), p(2:dim+1,:), p(dim+2:dim+ydim+1,:), p(dim+ydim+2:2*dim+ydim+1,:), p(2*dim+ydim+2:2*dim+2*ydim+1,:), p(2*dim+2*ydim+2:end,:), x);
                  txyxsysphi = [tj'; xj; yj; xsj; ysj; phij];
                  dpode      = coco_ezDFDX('f(x,p)v', f, pj, txyxsysphi);
                else
                  dpode = dgdphan(tj', xj, yj, xsj, ysj, phij, pj);
                end
                dprows = repmat(reshape(1:Nsa*gdim, [gdim Nsa]), [pdim 1]); 
                dpcols = repmat(1:pdim, [gdim Nsa]);  
                dpode  = sparse(dprows, dpcols, dpode(:)); 
                
                J_xbp = dxode*Pxj*seg.Wad+dxsode*Lx;
                J_ybp = dyode*Pyj*seg.Waa+dysode*Ly;
                J_T0  = dtode(:)-dxsode*DLx*x*dtauxdT0/T-dysode*DLy*y*dtauydT0/T+dphiode*Dphij(:)*(1-dtautdT0);               
                dtode = dtode.*repmat(sa', [gdim,1]);
                duode = Dphij.*repmat(sa'-dtautdT, [phidim,1]);
                J_T   = dtode(:)-dxsode*DLx*x*(dtauxdT*T-gseg.taux(T0,T,p))/T^2 ...
                        -dysode*DLy*y*(dtauydT*T-gseg.tauy(T0,T,p))/T^2+dphiode*duode(:);
                J_p   = dpode-dxsode*DLx*x*dtauxdp/T-dysode*DLy*y*dtauydp/T-dphiode*Dphij(:)*dtautdp;

                J(s+1:s+Nsa*gdim,:) = [J_xbp J_ybp J_T0, J_T, J_p];                                              
            else
                % Differention of g
                dgdthan  = [];
                dgdxhan  = [];
                dgdyhan  = [];
                dgdxshan = [];
                dgdyshan = [];
                dgdphan  = [];            
                switch numel(gseg.Dg)
                    case 1
                        dgdthan = gseg.Dg{1};
                    case 2
                        dgdthan = gseg.Dg{1};
                        dgdxhan = gseg.Dg{2};
                    case 3
                        dgdthan = gseg.Dg{1};
                        dgdxhan = gseg.Dg{2};
                        dgdyhan = gseg.Dg{3};
                    case 4
                        dgdthan  = gseg.Dg{1};
                        dgdxhan  = gseg.Dg{2};
                        dgdyhan  = gseg.Dg{3};
                        dgdxshan = gseg.Dg{4};
                    case 5
                        dgdthan  = gseg.Dg{1};
                        dgdxhan  = gseg.Dg{2};
                        dgdyhan  = gseg.Dg{3};
                        dgdxshan = gseg.Dg{4};
                        dgdyshan = gseg.Dg{5};
                    case 6
                        dgdthan  = gseg.Dg{1};
                        dgdxhan  = gseg.Dg{2};
                        dgdyhan  = gseg.Dg{3};
                        dgdxshan = gseg.Dg{4};
                        dgdyshan = gseg.Dg{5};
                        dgdphan  = gseg.Dg{6};
                end
                
                % w.r.t t
                if isempty(dgdthan)
                  f       = @(x,p) gseg.gfunc(x, p(1:dim,:), p(dim+1:ydim+dim,:), p(ydim+dim+1:ydim+2*dim,:), p(ydim+2*dim+1:2*ydim+2*dim,:), p(2*ydim+2*dim+1:end,:));
                  xyxsysp = [xj; yj; xsj; ysj; pj];
                  dtode   = coco_ezDFDX('f(x,p)v', f, tj', xyxsysp);
                  dtode   = reshape(dtode,[gdim,Nsa]);
                else
                  dtode   = dgdthan(tj', xj, yj, xsj, ysj, pj);
                end
                % w.r.t x
                if isempty(dgdxhan)
                  f       = @(x,p) gseg.gfunc(p(1,:), x, p(2:ydim+1,:), p(ydim+2:ydim+dim+1,:), p(ydim+dim+2:2*ydim+dim+1,:), p(2*ydim+dim+2:end,:));
                  tyxsysp = [tj'; yj; xsj; ysj; pj];
                  dxode   = coco_ezDFDX('f(x,p)v', f, xj, tyxsysp);
                else
                  dxode   = dgdxhan(tj', xj, yj, xsj, ysj, pj);
                end
                dxrows = repmat(reshape(1:Nsa*gdim, [gdim Nsa]), [dim 1]); 
                dxcols = repmat(1:Nsa*dim, [gdim 1]);  
                dxode  = sparse(dxrows, dxcols, dxode(:));                
                % w.r.t. y
                if isempty(dgdyhan)
                  f       = @(x,p) gseg.gfunc(p(1,:), p(2:dim+1,:), x, p(dim+2:2*dim+1,:), p(2*dim+2:2*dim+ydim+1,:), p(2*dim+ydim+2:end,:));
                  txxsysp = [tj'; xj; xsj; ysj; pj];
                  dyode   = coco_ezDFDX('f(x,p)v', f, yj, txxsysp);
                else
                  dyode   = dgdyhan(tj', xj, yj, xsj, ysj, pj);
                end
                dyrows = repmat(reshape(1:Nsa*gdim, [gdim Nsa]), [ydim 1]); 
                dycols = repmat(1:Nsa*ydim, [gdim 1]);  
                dyode  = sparse(dyrows, dycols, dyode(:));
                % w.r.t xs
                if isempty(dgdxshan)
                  f       = @(x,p) gseg.gfunc(p(1,:), p(2:dim+1,:), p(dim+2:ydim+dim+1,:), x, p(ydim+dim+2:2*ydim+dim+1,:), p(2*ydim+dim+2:end,:));
                  txyysp  = [tj'; xj; yj; ysj; pj];
                  dxsode  = coco_ezDFDX('f(x,p)v', f, xsj, txyysp);
                else
                  dxsode  = dgdxshan(tj', xj, yj, xsj, ysj, pj);
                end
                dxsode    = sparse(dxrows, dxcols, dxsode(:)); 
                % w.r.t. ys
                if isempty(dgdyshan)
                  f       = @(x,p) gseg.gfunc(p(1,:), p(2:dim+1,:), p(dim+2:dim+ydim+1,:), p(dim+ydim+2:2*dim+ydim+1,:), x, p(2*dim+ydim+2:end,:));
                  txyxsp  = [tj'; xj; yj; xsj; pj];
                  dysode  = coco_ezDFDX('f(x,p)v', f, ysj, txyxsp);
                else
                  dysode  = dgdyshan(tj', xj, yj, xsj, ysj, pj);
                end                
                dysode    = sparse(dyrows, dycols, dysode(:)); 
                % w.r.t. phi                       
                % w.r.t. p
                if isempty(dgdphan)
                  f       = @(x,p) gseg.gfunc(p(1,:), p(2:dim+1,:), p(dim+2:dim+ydim+1,:), p(dim+ydim+2:2*dim+ydim+1,:), p(2*dim+ydim+2:end,:), x);
                  txyxsys = [tj'; xj; yj; xsj; ysj];
                  dpode   = coco_ezDFDX('f(x,p)v', f, pj, txyxsys);
                else
                  dpode   = dgdphan(tj', xj, yj, xsj, ysj, pj);
                end
                dprows = repmat(reshape(1:Nsa*gdim, [gdim Nsa]), [pdim 1]); 
                dpcols = repmat(1:pdim, [gdim Nsa]);  
                dpode  = sparse(dprows, dpcols, dpode(:)); 
                
                J_xbp = dxode*Pxj*seg.Wad+dxsode*Lx;
                J_ybp = dyode*Pyj*seg.Waa+dysode*Ly;
                J_T0  = dtode(:)-dxsode*DLx*x*dtauxdT0/T-dysode*DLy*y*dtauydT0/T;               
                dtode = dtode.*repmat(sa', [gdim,1]);
                J_T   = dtode(:)-dxsode*DLx*x*(dtauxdT*T-gseg.taux(T0,T,p))/T^2 ...
                        -dysode*DLy*y*(dtauydT*T-gseg.tauy(T0,T,p))/T^2;
                J_p   = dpode-dxsode*DLx*x*dtauxdp/T-dysode*DLy*y*dtauydp/T;

                J(s+1:s+Nsa*gdim,:) = [J_xbp J_ybp J_T0, J_T, J_p];                                               
            end
        else
            if ~isempty(gseg.phi)
                % Differention of g
                dgdthan  = [];
                dgdxhan  = [];
                dgdyhan  = [];
                dgdxshan = [];
                dgdphihan = [];
                dgdphan  = [];            
                switch numel(gseg.Dg)
                    case 1
                        dgdthan = gseg.Dg{1};
                    case 2
                        dgdthan = gseg.Dg{1};
                        dgdxhan = gseg.Dg{2};
                    case 3
                        dgdthan = gseg.Dg{1};
                        dgdxhan = gseg.Dg{2};
                        dgdyhan = gseg.Dg{3};
                    case 4
                        dgdthan  = gseg.Dg{1};
                        dgdxhan  = gseg.Dg{2};
                        dgdyhan  = gseg.Dg{3};
                        dgdxshan = gseg.Dg{4};
                    case 5
                        dgdthan  = gseg.Dg{1};
                        dgdxhan  = gseg.Dg{2};
                        dgdyhan  = gseg.Dg{3};
                        dgdxshan = gseg.Dg{4};
                        dgdphihan = gseg.Dg{5};
                    case 6
                        dgdthan  = gseg.Dg{1};
                        dgdxhan  = gseg.Dg{2};
                        dgdyhan  = gseg.Dg{3};
                        dgdxshan = gseg.Dg{4};
                        dgdphihan = gseg.Dg{5};
                        dgdphan  = gseg.Dg{6};
                end
                
                % w.r.t t
                if isempty(dgdthan)
                  f        = @(x,p) gseg.gfunc(x, p(1:dim,:), p(dim+1:ydim+dim,:), p(ydim+dim+1:ydim+2*dim,:), p(ydim+2*dim+1:ydim+2*dim+phidim,:), p(ydim+2*dim+phidim+1:end,:));
                  xyxsphip = [xj; yj; xsj; phij; pj];
                  dtode    = coco_ezDFDX('f(x,p)v', f, tj', xyxsphip);
                  dtode    = reshape(dtode,[gdim,Nsa]);
                else
                  dtode    = dgdthan(tj', xj, yj, xsj, phij, pj);
                end
                % w.r.t x
                if isempty(dgdxhan)
                  f        = @(x,p) gseg.gfunc(p(1,:), x, p(2:ydim+1,:), p(ydim+2:ydim+dim+1,:), p(ydim+dim+2:ydim+dim+phidim+1,:), p(ydim+dim+phidim+2:end,:));
                  tyxsphip = [tj'; yj; xsj; phij; pj];
                  dxode    = coco_ezDFDX('f(x,p)v', f, xj, tyxsphip);
                else
                  dxode    = dgdxhan(tj', xj, yj, xsj, phij, pj);
                end
                dxrows = repmat(reshape(1:Nsa*gdim, [gdim Nsa]), [dim 1]); 
                dxcols = repmat(1:Nsa*dim, [gdim 1]);  
                dxode  = sparse(dxrows, dxcols, dxode(:));                
                % w.r.t. y
                if isempty(dgdyhan)
                  f        = @(x,p) gseg.gfunc(p(1,:), p(2:dim+1,:), x, p(dim+2:2*dim+1,:), p(2*dim+2:2*dim+phidim+1,:), p(2*dim+phidim+2:end,:));
                  txxsphip = [tj'; xj; xsj; phij; pj];
                  dyode    = coco_ezDFDX('f(x,p)v', f, yj, txxsphip);
                else
                  dyode    = dgdyhan(tj', xj, yj, xsj, phij, pj);
                end
                dyrows = repmat(reshape(1:Nsa*gdim, [gdim Nsa]), [ydim 1]); 
                dycols = repmat(1:Nsa*ydim, [gdim 1]);  
                dyode  = sparse(dyrows, dycols, dyode(:));
                % w.r.t xs
                if isempty(dgdxshan)
                  f        = @(x,p) gseg.gfunc(p(1,:), p(2:dim+1,:), p(dim+2:ydim+dim+1,:), x, p(ydim+dim+2:ydim+dim+phidim+1,:), p(ydim+dim+phidim+2:end,:));
                  txyphip  = [tj'; xj; yj; phij; pj];
                  dxsode   = coco_ezDFDX('f(x,p)v', f, xsj, txyphip);
                else
                  dxsode   = dgdxshan(tj', xj, yj, xsj, phij, pj);
                end
                dxsode   = sparse(dxrows, dxcols, dxsode(:)); 
                % w.r.t. ys
                % w.r.t. phi           
                if isempty(dgdphihan)
                  f        = @(x,p) gseg.gfunc(p(1,:), p(2:dim+1,:), p(dim+2:dim+ydim+1,:), p(dim+ydim+2:2*dim+ydim+1,:), x, p(2*dim+ydim+2:end,:));
                  txyxsp   = [tj'; xj; yj; xsj; pj];
                  dphiode  = coco_ezDFDX('f(x,p)v', f, phij, txyxsp);
                else
                  dphiode  = dgdphihan(tj', xj, yj, xsj, phij, pj);
                end  
                dphirows = repmat(reshape(1:Nsa*gdim, [gdim Nsa]), [phidim 1]); 
                dphicols = repmat(1:Nsa*phidim, [gdim 1]);  
                dphiode  = sparse(dphirows, dphicols, dphiode(:));             
                % w.r.t. p
                if isempty(dgdphan)
                  f        = @(x,p) gseg.gfunc(p(1,:), p(2:dim+1,:), p(dim+2:dim+ydim+1,:), p(dim+ydim+2:2*dim+ydim+1,:), p(2*dim+ydim+2:end,:), x);
                  txyxsphi = [tj'; xj; yj; xsj; phij];
                  dpode    = coco_ezDFDX('f(x,p)v', f, pj, txyxsphi);
                else
                  dpode    = dgdphan(tj', xj, yj, xsj, phij, pj);
                end
                dprows = repmat(reshape(1:Nsa*gdim, [gdim Nsa]), [pdim 1]); 
                dpcols = repmat(1:pdim, [gdim Nsa]);  
                dpode  = sparse(dprows, dpcols, dpode(:)); 
                
                J_xbp = dxode*Pxj*seg.Wad+dxsode*Lx;
                J_ybp = dyode*Pyj*seg.Waa;
                J_T0  = dtode(:)-dxsode*DLx*x*dtauxdT0/T+dphiode*Dphij(:)*(1-dtautdT0);               
                dtode = dtode.*repmat(sa', [gdim,1]);
                duode = Dphij.*repmat(sa'-dtautdT, [phidim,1]);
                J_T   = dtode(:)-dxsode*DLx*x*(dtauxdT*T-gseg.taux(T0,T,p))/T^2+dphiode*duode(:);
                J_p   = dpode-dxsode*DLx*x*dtauxdp/T-dphiode*Dphij(:)*dtautdp;

                J(s+1:s+Nsa*gdim,:) = [J_xbp J_ybp J_T0, J_T, J_p];                                                               
            else
                % Differention of g
                dgdthan  = [];
                dgdxhan  = [];
                dgdyhan  = [];
                dgdxshan = [];
                dgdphan  = [];            
                switch numel(gseg.Dg)
                    case 1
                        dgdthan = gseg.Dg{1};
                    case 2
                        dgdthan = gseg.Dg{1};
                        dgdxhan = gseg.Dg{2};
                    case 3
                        dgdthan = gseg.Dg{1};
                        dgdxhan = gseg.Dg{2};
                        dgdyhan = gseg.Dg{3};
                    case 4
                        dgdthan  = gseg.Dg{1};
                        dgdxhan  = gseg.Dg{2};
                        dgdyhan  = gseg.Dg{3};
                        dgdxshan = gseg.Dg{4};
                    case 5
                        dgdthan  = gseg.Dg{1};
                        dgdxhan  = gseg.Dg{2};
                        dgdyhan  = gseg.Dg{3};
                        dgdxshan = gseg.Dg{4};
                        dgdphan  = gseg.Dg{5};
                end
                
                % w.r.t t
                if isempty(dgdthan)
                  f     = @(x,p) gseg.gfunc(x, p(1:dim,:), p(dim+1:ydim+dim,:), p(ydim+dim+1:ydim+2*dim,:), p(ydim+2*dim+1:end,:));
                  xyxsp = [xj; yj; xsj; pj];
                  dtode = coco_ezDFDX('f(x,p)v', f, tj', xyxsp);
                  dtode = reshape(dtode,[gdim,Nsa]);
                else
                  dtode = dgdthan(tj', xj, yj, xsj, pj);
                end
                % w.r.t x
                if isempty(dgdxhan)
                  f     = @(x,p) gseg.gfunc(p(1,:), x, p(2:ydim+1,:), p(ydim+2:ydim+dim+1,:), p(ydim+dim+2:end,:));
                  tyxsp = [tj'; yj; xsj; pj];
                  dxode = coco_ezDFDX('f(x,p)v', f, xj, tyxsp);
                else
                  dxode = dgdxhan(tj', xj, yj, xsj, pj);
                end
                dxrows = repmat(reshape(1:Nsa*gdim, [gdim Nsa]), [dim 1]); 
                dxcols = repmat(1:Nsa*dim, [gdim 1]);  
                dxode  = sparse(dxrows, dxcols, dxode(:));                
                % w.r.t. y
                if isempty(dgdyhan)
                  f     = @(x,p) gseg.gfunc(p(1,:), p(2:dim+1,:), x, p(dim+2:2*dim+1,:), p(2*dim+2:end,:));
                  txxsp = [tj'; xj; xsj; pj];
                  dyode = coco_ezDFDX('f(x,p)v', f, yj, txxsp);
                else
                  dyode = dgdyhan(tj', xj, yj, xsj, pj);
                end
                dyrows = repmat(reshape(1:Nsa*gdim, [gdim Nsa]), [ydim 1]); 
                dycols = repmat(1:Nsa*ydim, [gdim 1]);  
                dyode  = sparse(dyrows, dycols, dyode(:));
                % w.r.t xs
                if isempty(dgdxshan)
                  f      = @(x,p) gseg.gfunc(p(1,:), p(2:dim+1,:), p(dim+2:ydim+dim+1,:), x, p(ydim+dim+2:end,:));
                  txyp   = [tj'; xj; yj; pj];
                  dxsode = coco_ezDFDX('f(x,p)v', f, xsj, txyp);
                else
                  dxsode   = dgdxshan(tj', xj, yj, xsj, pj);
                end
                dxsode   = sparse(dxrows, dxcols, dxsode(:)); 
                % w.r.t. ys
                % w.r.t. phi                     
                % w.r.t. p
                if isempty(dgdphan)
                  f     = @(x,p) gseg.gfunc(p(1,:), p(2:dim+1,:), p(dim+2:dim+ydim+1,:), p(dim+ydim+2:end,:), x);
                  txyxs = [tj'; xj; yj; xsj];
                  dpode = coco_ezDFDX('f(x,p)v', f, pj, txyxs);
                else
                  dpode = dgdphan(tj', xj, yj, xsj, pj);
                end
                dprows = repmat(reshape(1:Nsa*gdim, [gdim Nsa]), [pdim 1]); 
                dpcols = repmat(1:pdim, [gdim Nsa]);  
                dpode  = sparse(dprows, dpcols, dpode(:)); 
                
                J_xbp = dxode*Pxj*seg.Wad+dxsode*Lx;
                J_ybp = dyode*Pyj*seg.Waa;
                J_T0  = dtode(:)-dxsode*DLx*x*dtauxdT0/T;               
                dtode = dtode.*repmat(sa', [gdim,1]);
                J_T   = dtode(:)-dxsode*DLx*x*(dtauxdT*T-gseg.taux(T0,T,p))/T^2;
                J_p   = dpode-dxsode*DLx*x*dtauxdp/T;

                J(s+1:s+Nsa*gdim,:) = [J_xbp J_ybp J_T0, J_T, J_p];                                                               
                
            end
        end
    else
        if ~isempty(gseg.tauy)
            if ~isempty(gseg.phi)
                % Differention of g
                dgdthan  = [];
                dgdxhan  = [];
                dgdyhan  = [];
                dgdyshan = [];
                dgdphihan = [];
                dgdphan  = [];            
                switch numel(gseg.Dg)
                    case 1
                        dgdthan = gseg.Dg{1};
                    case 2
                        dgdthan = gseg.Dg{1};
                        dgdxhan = gseg.Dg{2};
                    case 3
                        dgdthan = gseg.Dg{1};
                        dgdxhan = gseg.Dg{2};
                        dgdyhan = gseg.Dg{3};
                    case 4
                        dgdthan  = gseg.Dg{1};
                        dgdxhan  = gseg.Dg{2};
                        dgdyhan  = gseg.Dg{3};
                        dgdyshan = gseg.Dg{4};
                    case 5
                        dgdthan  = gseg.Dg{1};
                        dgdxhan  = gseg.Dg{2};
                        dgdyhan  = gseg.Dg{3};
                        dgdyshan = gseg.Dg{4};
                        dgdphihan = gseg.Dg{5};
                    case 6
                        dgdthan  = gseg.Dg{1};
                        dgdxhan  = gseg.Dg{2};
                        dgdyhan  = gseg.Dg{3};
                        dgdyshan = gseg.Dg{4};
                        dgdphihan = gseg.Dg{5};
                        dgdphan  = gseg.Dg{6};
                end
                
                % w.r.t t
                if isempty(dgdthan)
                  f        = @(x,p) gseg.gfunc(x, p(1:dim,:), p(dim+1:ydim+dim,:), p(ydim+dim+1:2*ydim+dim,:), p(2*ydim+dim+1:2*ydim+dim+phidim,:), p(2*ydim+dim+phidim+1:end,:));
                  xyysphip = [xj; yj; ysj; phij; pj];
                  dtode    = coco_ezDFDX('f(x,p)v', f, tj', xyysphip);
                  dtode    = reshape(dtode,[gdim,Nsa]);
                else
                  dtode    = dgdthan(tj', xj, yj, ysj, phij, pj);
                end
                % w.r.t x
                if isempty(dgdxhan)
                  f        = @(x,p) gseg.gfunc(p(1,:), x, p(2:ydim+1,:), p(ydim+2:2*ydim+1,:), p(2*ydim+2:2*ydim+phidim+1,:), p(2*ydim+phidim+2:end,:));
                  tyysphip = [tj'; yj; ysj; phij; pj];
                  dxode    = coco_ezDFDX('f(x,p)v', f, xj, tyysphip);
                else
                  dxode    = dgdxhan(tj', xj, yj, ysj, phij, pj);
                end
                dxrows = repmat(reshape(1:Nsa*gdim, [gdim Nsa]), [dim 1]); 
                dxcols = repmat(1:Nsa*dim, [gdim 1]);  
                dxode  = sparse(dxrows, dxcols, dxode(:));                
                % w.r.t. y
                if isempty(dgdyhan)
                  f          = @(x,p) gseg.gfunc(p(1,:), p(2:dim+1,:), x, p(dim+2:dim+ydim+1,:), p(dim+ydim+2:dim+ydim+phidim+1,:), p(dim+ydim+phidim+2:end,:));
                  txysphip = [tj'; xj; ysj; phij; pj];
                  dyode    = coco_ezDFDX('f(x,p)v', f, yj, txysphip);
                else
                  dyode    = dgdyhan(tj', xj, yj, ysj, phij, pj);
                end
                dyrows = repmat(reshape(1:Nsa*gdim, [gdim Nsa]), [ydim 1]); 
                dycols = repmat(1:Nsa*ydim, [gdim 1]);  
                dyode  = sparse(dyrows, dycols, dyode(:));
                % w.r.t. xs 
                % w.r.t. ys
                if isempty(dgdyshan)
                  f        = @(x,p) gseg.gfunc(p(1,:), p(2:dim+1,:), p(dim+2:dim+ydim+1,:), x, p(dim+ydim+2:dim+ydim+phidim+1,:), p(dim+ydim+phidim+2:end,:));
                  txyphip  = [tj'; xj; yj; phij; pj];
                  dysode   = coco_ezDFDX('f(x,p)v', f, ysj, txyphip);
                else
                  dysode   = dgdyshan(tj', xj, yj, ysj, phij, pj);
                end                
                dysode     = sparse(dyrows, dycols, dysode(:)); 
                % w.r.t. phi           
                if isempty(dgdphihan)
                  f        = @(x,p) gseg.gfunc(p(1,:), p(2:dim+1,:), p(dim+2:dim+ydim+1,:), p(dim+ydim+2:dim+2*ydim+1,:), x, p(dim+2*ydim+2:end,:));
                  txyysp   = [tj'; xj; yj; ysj; pj];
                  dphiode  = coco_ezDFDX('f(x,p)v', f, phij, txyysp);
                else
                  dphiode  = dgdphihan(tj', xj, yj, ysj, phij, pj);
                end  
                dphirows = repmat(reshape(1:Nsa*gdim, [gdim Nsa]), [phidim 1]); 
                dphicols = repmat(1:Nsa*phidim, [gdim 1]);  
                dphiode  = sparse(dphirows, dphicols, dphiode(:));             
                % w.r.t. p
                if isempty(dgdphan)
                  f        = @(x,p) gseg.gfunc(p(1,:), p(2:dim+1,:), p(dim+2:dim+ydim+1,:), p(dim+ydim+2:dim+2*ydim+1,:), p(dim+2*ydim+2:end,:), x);
                  txyysphi = [tj'; xj; yj; ysj; phij];
                  dpode    = coco_ezDFDX('f(x,p)v', f, pj, txyysphi);
                else
                  dpode = dgdphan(tj', xj, yj, ysj, phij, pj);
                end
                dprows = repmat(reshape(1:Nsa*gdim, [gdim Nsa]), [pdim 1]); 
                dpcols = repmat(1:pdim, [gdim Nsa]);  
                dpode  = sparse(dprows, dpcols, dpode(:)); 
                
                J_xbp = dxode*Pxj*seg.Wad;
                J_ybp = dyode*Pyj*seg.Waa+dysode*Ly;
                J_T0  = dtode(:)-dysode*DLy*y*dtauydT0/T+dphiode*Dphij(:)*(1-dtautdT0);               
                dtode = dtode.*repmat(sa', [gdim,1]);
                duode = Dphij.*repmat(sa'-dtautdT, [phidim,1]);
                J_T   = dtode(:)-dysode*DLy*y*(dtauydT*T-gseg.tauy(T0,T,p))/T^2+dphiode*duode(:);
                J_p   = dpode-dysode*DLy*y*dtauydp/T-dphiode*Dphij(:)*dtautdp;

                J(s+1:s+Nsa*gdim,:) = [J_xbp J_ybp J_T0, J_T, J_p];                                               
            else
                % Differention of g
                dgdthan  = [];
                dgdxhan  = [];
                dgdyhan  = [];
                dgdyshan = [];
                dgdphan  = [];            
                switch numel(gseg.Dg)
                    case 1
                        dgdthan = gseg.Dg{1};
                    case 2
                        dgdthan = gseg.Dg{1};
                        dgdxhan = gseg.Dg{2};
                    case 3
                        dgdthan = gseg.Dg{1};
                        dgdxhan = gseg.Dg{2};
                        dgdyhan = gseg.Dg{3};
                    case 4
                        dgdthan  = gseg.Dg{1};
                        dgdxhan  = gseg.Dg{2};
                        dgdyhan  = gseg.Dg{3};
                        dgdyshan = gseg.Dg{4};
                    case 5
                        dgdthan  = gseg.Dg{1};
                        dgdxhan  = gseg.Dg{2};
                        dgdyhan  = gseg.Dg{3};
                        dgdyshan = gseg.Dg{4};
                        dgdphan  = gseg.Dg{5};
                end
                
                % w.r.t t
                if isempty(dgdthan)
                  f     = @(x,p) gseg.gfunc(x, p(1:dim,:), p(dim+1:ydim+dim,:), p(ydim+dim+1:2*ydim+dim,:), p(2*ydim+dim+1:end,:));
                  xyysp = [xj; yj; ysj; pj];
                  dtode = coco_ezDFDX('f(x,p)v', f, tj', xyysp);
                  dtode = reshape(dtode,[gdim,Nsa]);
                else
                  dtode = dgdthan(tj', xj, yj, ysj, pj);
                end
                % w.r.t x
                if isempty(dgdxhan)
                  f     = @(x,p) gseg.gfunc(p(1,:), x, p(2:ydim+1,:), p(ydim+2:2*ydim+1,:), p(2*ydim+2:end,:));
                  tyysp = [tj'; yj; ysj; pj];
                  dxode = coco_ezDFDX('f(x,p)v', f, xj, tyysp);
                else
                  dxode = dgdxhan(tj', xj, yj, ysj, pj);
                end
                dxrows = repmat(reshape(1:Nsa*gdim, [gdim Nsa]), [dim 1]); 
                dxcols = repmat(1:Nsa*dim, [gdim 1]);  
                dxode  = sparse(dxrows, dxcols, dxode(:));                
                % w.r.t. y
                if isempty(dgdyhan)
                  f     = @(x,p) gseg.gfunc(p(1,:), p(2:dim+1,:), x, p(dim+2:dim+ydim+1,:), p(dim+ydim+2:end,:));
                  txysp = [tj'; xj; ysj; pj];
                  dyode = coco_ezDFDX('f(x,p)v', f, yj, txysp);
                else
                  dyode = dgdyhan(tj', xj, yj, ysj, pj);
                end
                dyrows = repmat(reshape(1:Nsa*gdim, [gdim Nsa]), [ydim 1]); 
                dycols = repmat(1:Nsa*ydim, [gdim 1]);  
                dyode  = sparse(dyrows, dycols, dyode(:));
                % w.r.t. xs 
                % w.r.t. ys
                if isempty(dgdyshan)
                  f      = @(x,p) gseg.gfunc(p(1,:), p(2:dim+1,:), p(dim+2:dim+ydim+1,:), x, p(dim+ydim+2:end,:));
                  txyp   = [tj'; xj; yj; pj];
                  dysode = coco_ezDFDX('f(x,p)v', f, ysj, txyp);
                else
                  dysode = dgdyshan(tj', xj, yj, ysj, pj);
                end                
                dysode   = sparse(dyrows, dycols, dysode(:)); 
                % w.r.t. phi                      
                % w.r.t. p
                if isempty(dgdphan)
                  f     = @(x,p) gseg.gfunc(p(1,:), p(2:dim+1,:), p(dim+2:dim+ydim+1,:), p(dim+ydim+2:end,:), x);
                  txyys = [tj'; xj; yj; ysj];
                  dpode = coco_ezDFDX('f(x,p)v', f, pj, txyys);
                else
                  dpode = dgdphan(tj', xj, yj, ysj, pj);
                end
                dprows = repmat(reshape(1:Nsa*gdim, [gdim Nsa]), [pdim 1]); 
                dpcols = repmat(1:pdim, [gdim Nsa]);  
                dpode  = sparse(dprows, dpcols, dpode(:)); 
                
                J_xbp = dxode*Pxj*seg.Wad;
                J_ybp = dyode*Pyj*seg.Waa+dysode*Ly;
                J_T0  = dtode(:)-dysode*DLy*y*dtauydT0/T;               
                dtode = dtode.*repmat(sa', [gdim,1]);
                J_T   = dtode(:)-dysode*DLy*y*(dtauydT*T-gseg.tauy(T0,T,p))/T^2;
                J_p   = dpode-dysode*DLy*y*dtauydp/T;

                J(s+1:s+Nsa*gdim,:) = [J_xbp J_ybp J_T0, J_T, J_p];                                               
                
            end
        else
            if ~isempty(gseg.phi)
                % Differention of g
                dgdthan  = [];
                dgdxhan  = [];
                dgdyhan  = [];
                dgdphihan = [];
                dgdphan  = [];            
                switch numel(gseg.Dg)
                    case 1
                        dgdthan = gseg.Dg{1};
                    case 2
                        dgdthan = gseg.Dg{1};
                        dgdxhan = gseg.Dg{2};
                    case 3
                        dgdthan = gseg.Dg{1};
                        dgdxhan = gseg.Dg{2};
                        dgdyhan = gseg.Dg{3};
                    case 4
                        dgdthan  = gseg.Dg{1};
                        dgdxhan  = gseg.Dg{2};
                        dgdyhan  = gseg.Dg{3};
                        dgdphihan = gseg.Dg{4};
                    case 5
                        dgdthan  = gseg.Dg{1};
                        dgdxhan  = gseg.Dg{2};
                        dgdyhan  = gseg.Dg{3};
                        dgdphihan = gseg.Dg{4};
                        dgdphan  = gseg.Dg{5};
                end
                
                % w.r.t t
                if isempty(dgdthan)
                  f      = @(x,p) gseg.gfunc(x, p(1:dim,:), p(dim+1:ydim+dim,:), p(ydim+dim+1:ydim+dim+phidim,:), p(ydim+dim+phidim+1:end,:));
                  xyphip = [xj; yj; phij; pj];
                  dtode  = coco_ezDFDX('f(x,p)v', f, tj', xyphip);
                  dtode  = reshape(dtode,[gdim,Nsa]);
                else
                  dtode  = dgdthan(tj', xj, yj, phij, pj);
                end
                % w.r.t x
                if isempty(dgdxhan)
                  f      = @(x,p) gseg.gfunc(p(1,:), x, p(2:ydim+1,:), p(ydim+2:ydim+phidim+1,:), p(ydim+phidim+2:end,:));
                  typhip = [tj'; yj; phij; pj];
                  dxode  = coco_ezDFDX('f(x,p)v', f, xj, typhip);
                else
                  dxode    = dgdxhan(tj', xj, yj, phij, pj);
                end
                dxrows = repmat(reshape(1:Nsa*gdim, [gdim Nsa]), [dim 1]); 
                dxcols = repmat(1:Nsa*dim, [gdim 1]);  
                dxode  = sparse(dxrows, dxcols, dxode(:));                
                % w.r.t. y
                if isempty(dgdyhan)
                  f      = @(x,p) gseg.gfunc(p(1,:), p(2:dim+1,:), x, p(dim+2:dim+phidim+1,:), p(dim+phidim+2:end,:));
                  txphip = [tj'; xj; phij; pj];
                  dyode  = coco_ezDFDX('f(x,p)v', f, yj, txphip);
                else
                  dyode    = dgdyhan(tj', xj, yj, phij, pj);
                end
                dyrows = repmat(reshape(1:Nsa*gdim, [gdim Nsa]), [ydim 1]); 
                dycols = repmat(1:Nsa*ydim, [gdim 1]);  
                dyode  = sparse(dyrows, dycols, dyode(:));
                % w.r.t. xs 
                % w.r.t. ys
                % w.r.t. phi           
                if isempty(dgdphihan)
                  f       = @(x,p) gseg.gfunc(p(1,:), p(2:dim+1,:), p(dim+2:dim+ydim+1,:), x, p(dim+ydim+2:end,:));
                  txyp    = [tj'; xj; yj; pj];
                  dphiode = coco_ezDFDX('f(x,p)v', f, phij, txyp);
                else
                  dphiode  = dgdphihan(tj', xj, yj, phij, pj);
                end  
                dphirows = repmat(reshape(1:Nsa*gdim, [gdim Nsa]), [phidim 1]); 
                dphicols = repmat(1:Nsa*phidim, [gdim 1]);  
                dphiode  = sparse(dphirows, dphicols, dphiode(:));             
                % w.r.t. p
                if isempty(dgdphan)
                  f      = @(x,p) gseg.gfunc(p(1,:), p(2:dim+1,:), p(dim+2:dim+ydim+1,:), p(dim+ydim+2:end,:), x);
                  txyphi = [tj'; xj; yj; phij];
                  dpode  = coco_ezDFDX('f(x,p)v', f, pj, txyphi);
                else
                  dpode = dgdphan(tj', xj, yj, phij, pj);
                end
                dprows = repmat(reshape(1:Nsa*gdim, [gdim Nsa]), [pdim 1]); 
                dpcols = repmat(1:pdim, [gdim Nsa]);  
                dpode  = sparse(dprows, dpcols, dpode(:)); 
                
                J_xbp = dxode*Pxj*seg.Wad;
                J_ybp = dyode*Pyj*seg.Waa;
                J_T0  = dtode(:)+dphiode*Dphij(:)*(1-dtautdT0);               
                dtode = dtode.*repmat(sa', [gdim,1]);
                duode = Dphij.*repmat(sa'-dtautdT, [phidim,1]);
                J_T   = dtode(:)+dphiode*duode(:);
                J_p   = dpode-dphiode*Dphij(:)*dtautdp;

                J(s+1:s+Nsa*gdim,:) = [J_xbp J_ybp J_T0, J_T, J_p];                                               
                
            else
                % Differention of g
                dgdthan  = [];
                dgdxhan  = [];
                dgdyhan  = [];
                dgdphan  = [];            
                switch numel(gseg.Dg)
                    case 1
                        dgdthan = gseg.Dg{1};
                    case 2
                        dgdthan = gseg.Dg{1};
                        dgdxhan = gseg.Dg{2};
                    case 3
                        dgdthan = gseg.Dg{1};
                        dgdxhan = gseg.Dg{2};
                        dgdyhan = gseg.Dg{3};
                    case 4
                        dgdthan  = gseg.Dg{1};
                        dgdxhan  = gseg.Dg{2};
                        dgdyhan  = gseg.Dg{3};
                        dgdphan  = gseg.Dg{4};
                end
                
                % w.r.t t
                if isempty(dgdthan)
                  f     = @(x,p) gseg.gfunc(x, p(1:dim,:), p(dim+1:ydim+dim,:), p(ydim+dim+1:end,:));
                  xyp   = [xj; yj; pj];
                  dtode = coco_ezDFDX('f(x,p)v', f, tj', xyp);
                  dtode = reshape(dtode,[gdim,Nsa]);
                else
                  dtode = dgdthan(tj', xj, yj, pj);
                end
                % w.r.t x
                if isempty(dgdxhan)
                  f     = @(x,p) gseg.gfunc(p(1,:), x, p(2:ydim+1,:), p(ydim+2:end,:));
                  typ   = [tj'; yj; pj];
                  dxode = coco_ezDFDX('f(x,p)v', f, xj, typ);
                else
                  dxode = dgdxhan(tj', xj, yj, pj);
                end
                dxrows = repmat(reshape(1:Nsa*gdim, [gdim Nsa]), [dim 1]); 
                dxcols = repmat(1:Nsa*dim, [gdim 1]);  
                dxode  = sparse(dxrows, dxcols, dxode(:));                
                % w.r.t. y
                if isempty(dgdyhan)
                  f     = @(x,p) gseg.gfunc(p(1,:), p(2:dim+1,:), x, p(dim+2:end,:));
                  txp   = [tj'; xj; pj];
                  dyode = coco_ezDFDX('f(x,p)v', f, yj, txp);
                else
                  dyode = dgdyhan(tj', xj, yj, pj);
                end
                dyrows = repmat(reshape(1:Nsa*gdim, [gdim Nsa]), [ydim 1]); 
                dycols = repmat(1:Nsa*ydim, [gdim 1]);  
                dyode  = sparse(dyrows, dycols, dyode(:));
                % w.r.t. xs 
                % w.r.t. ys
                % w.r.t. phi                      
                % w.r.t. p
                if isempty(dgdphan)
                  f     = @(x,p) gseg.gfunc(p(1,:), p(2:dim+1,:), p(dim+2:end,:), x);
                  txy   = [tj'; xj; yj];
                  dpode = coco_ezDFDX('f(x,p)v', f, pj, txy);
                else
                  dpode = dgdphan(tj', xj, yj, pj);
                end
                dprows = repmat(reshape(1:Nsa*gdim, [gdim Nsa]), [pdim 1]); 
                dpcols = repmat(1:pdim, [gdim Nsa]);  
                dpode  = sparse(dprows, dpcols, dpode(:)); 
                
                J_xbp = dxode*Pxj*seg.Wad;
                J_ybp = dyode*Pyj*seg.Waa;
                J_T0  = dtode(:);               
                dtode = dtode.*repmat(sa', [gdim,1]);
                J_T   = dtode(:);
                J_p   = dpode;

                J(s+1:s+Nsa*gdim,:) = [J_xbp J_ybp J_T0, J_T, J_p];                                               
                
            end
        end
    end
    s = s+Nsa*gdim;
%     toc
end

% [data, DadF] = coco_ezDFDX('f(o,d,x)',  prob, data, @alg_ddae_F, u);
% assert(max(max(abs(J-DadF)))<1e-5, 'The Jacobian is incorrect');
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
