function [data y] = coup(prob, data, u)


xbp = u(data.xbp_idx); % Extract basepoint values
ybp = u(data.ybp_idx);
xcn = reshape(data.Wdd*xbp, data.x_shp); % Values at collocation nodes
ycn = reshape(data.Wda*ybp, data.y_shp);
x1 = xcn(2,:);
y  = x1-ycn;
y  = y(:);



end