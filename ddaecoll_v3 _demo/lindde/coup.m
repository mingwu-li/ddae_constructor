function [data y] = coup(prob, data, u)


xbp = u(data.xbp_idx); % Extract basepoint values
ybp = u(data.ybp_idx);
T0  = u(data.T0_idx);
T   = u(data.T_idx);
tcn = data.Md*T+T0;
xbp = reshape(xbp, data.xbp_shp); % Values at collocation nodes
% ycn = reshape(ycn,data.y_shp);
ycn = reshape(data.Wda*ybp, data.y_shp);


x1 = xbp(1,:);
x2 = xbp(2,:);

NTST = data.ddaecoll.NTST;
NCOL = data.ddaecoll.NCOL;
xcnnum = NTST*NCOL;
tm = linspace(-1, 1, NCOL+1)';

y1 = zeros(1,xcnnum);
for k=1:xcnnum
    tk = tcn(k);
    if tk<=1+T0
        y1(k) = 1;
    else
        tauk = (tk-1-T0)/T;
        % find the interval index
        id   = ceil(tauk*NTST);
        % extract x at base points of this interval
        xbpk = x1((id-1)*(NCOL+1)+1:id*(NCOL+1));
        % extract interpolation matrices
        Lk   = coll_L((1+tm)/(2*NTST)+(id-1)/NTST,tauk);
        y1(k)= Lk*xbpk';
    end
end

y2 = zeros(1,xcnnum);
for k=1:xcnnum
    tk = tcn(k);
    if tk<=0.2+T0
        y2(k) = 1;
    else
        tauk = (tk-0.2-T0)/T;
        % find the interval index
        id   = ceil(tauk*NTST);
        % extract x at base points of this interval
        xbpk = x2((id-1)*(NCOL+1)+1:id*(NCOL+1));
        % extract interpolation matrices
        Lk   = coll_L((1+tm)/(2*NTST)+(id-1)/NTST,tauk);
        y2(k)= Lk*xbpk';
    end
end

y  = [y1;y2]-ycn;
y  = y(:);

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
