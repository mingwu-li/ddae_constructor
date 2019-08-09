function catenary_no_y_optimization

coco_use_recipes_toolbox ddaecoll_v2 adjt_ddaecoll_v2

prob = coco_prob;


%% first run to find initial fold
% zero problems
t0 = [0; 1];
x0 = [1 0; 1 0];
Y0 = 0.8;
coll_args = {@catn, @catn_dt, @catn_dx, @catn_dy, @catn_dp, ...
  t0, x0, [], 'Y', Y0};
prob1 = ddaecoll_isol2seg(prob, '', coll_args{:});
[data, uidx] = coco_get_func_data(prob1, 'ddaecoll', 'data', 'uidx');
maps = data;
prob1 = coco_add_func(prob1, 'bc', @catn_bc, @catn_bc_du, ...
  data, 'zero', 'uidx', ...
  uidx([maps.T_idx; maps.x0_idx; maps.x1_idx; maps.p_idx; maps.T0_idx]));
prob1 = coco_add_pars(prob1, 'rightend', uidx(maps.x1_idx(1)), 'y1');
% adjoints
prob1 = adjt_ddaecoll_isol2seg(prob1, '');
[data, axidx] = coco_get_adjt_data(prob1, 'ddaecoll', 'data', 'axidx');
opt = data.ddaecoll_opt;
prob1 = coco_add_adjt(prob1, 'bc', 'aidx', ...
  axidx([opt.T_idx; opt.x0_idx; opt.x1_idx; opt.p_idx; opt.T0_idx]));

prob1 = coco_add_adjt(prob1, 'rightend', 'd.y1', 'aidx', axidx(opt.x1_idx(1)));

% bd1 = coco(prob1, 'catn1', [], 1, {'Y' 'd.Y'}, [0.9 2]);
bd1 = coco(prob1, 'catn1', [], 1, {'y1', 'd.y1', 'Y'}, [0 1.5]);

%% branch switch from fold to grow nontrivial adjoint
BPlab = coco_bd_labs(bd1, 'BP');
% zero problems
prob2 = ode_BP2coll(prob, '', 'catn1', BPlab(1));    
[data, uidx] = coco_get_func_data(prob2, 'coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob2 = coco_add_func(prob2, 'bc', @catn_bc, @catn_bc_du, ...
  @catn_bc_dudu, data, 'zero', 'uidx', ...
  uidx([maps.T_idx; maps.x0_idx; maps.x1_idx; maps.p_idx]));    
prob2 = coco_add_pars(prob2, 'rightend', uidx(maps.x1_idx(1)), 'y1');
% branch switch data
chart = coco_read_solution('catn1', BPlab(1), 'chart');
cdata = coco_get_chart_data(chart, 'lsol');
% adjoints
prob2 = adjt_BP2coll(prob2, '', 'catn1', BPlab(1));
[chart, aidx] = coco_read_adjoint('bc', 'catn1', BPlab(1), 'chart', 'aidx');
[data, axidx] = coco_get_adjt_data(prob2, 'coll', 'data', 'axidx');
opt = data.coll_opt;
prob2   = coco_add_adjt(prob2, 'bc', 'aidx', ...
  axidx([opt.T_idx; opt.x0_idx; opt.x1_idx; opt.p_idx]), ...
  'l0', chart.x, 'tl0', cdata.v(aidx));
[chart, aidx] = coco_read_adjoint('rightend', 'catn1', BPlab(1), 'chart', 'aidx');
prob2 = coco_add_adjt(prob2, 'rightend', 'd.y1', 'aidx', ...
  axidx(opt.x1_idx(1)), 'l0', chart.x, 'tl0', cdata.v(aidx));
dv_int = [chart.x(1) 1];

bd2 = coco(prob2, 'catn2', [], {'d.y1', 'y1', 'd.coll.T0', 'Y'}, {dv_int});

profile off

profile viewer

%% check the result
%% plot optimal results
lab = coco_bd_labs(bd2, 'EP');
sol = coll_read_solution('', 'catn2', lab(2));
x1=sol.xbp(:,1);
x2=sol.xbp(:,2);
t=sol.tbp;
figure(2)
plot(sol.tbp,sol.xbp(:,1),'r*'); hold on
plot(sol.tbp,sol.xbp(:,2),'b*');
legend('x_1, test code','x_2, test code');

sol = coll_read_adjoint('', 'catn2', lab(2));
[data, uidx] = coco_get_func_data(prob2, 'coll', 'data', 'uidx');
maps = data.coll_seg.maps;
lbp     = reshape(sol.l0(maps.xbp_idx), maps.xbp_shp)';
sol.lbp = lbp(maps.tbp_idx,:);
l1 = sol.lbp(:,1);
l2 = sol.lbp(:,2);

figure(3)
plot(t,l1,'r*'); hold on
plot(t,l2,'b*');
legend('\lambda_1, my code','\lambda_2, my code');


sol.t=t;
yy=l1.*x2+l2.*(1+x2.^2)./x1;
[chart, aidx] = coco_read_adjoint('bc', 'catn2', lab(2), 'chart', 'aidx');
l3=trapz(sol.t,yy); % done see chart.x(1)
l4=l1(1);           % done see chart.x(2)
l5=0;               % done see chart.x(3)

l1dot=l2.*(1+x2.^2)./x1.^2;
l2dot=-(l1+2*l2.*x2./x1);
dl1=(l1(3:end)-l1(1:end-2))/(sol.t(2)-sol.t(1))/2; % central difference
dl2=(l2(3:end)-l2(1:end-2))/(sol.t(2)-sol.t(1))/2; % central difference
figure(4)
plot(sol.t,l1dot);hold on
plot(sol.t(2:end-1),dl1,'ro')
plot(sol.t,l2dot);
plot(sol.t(2:end-1),dl2,'ko')


 
end  



function y = catn(t, x, y, p)

x1 = x(1,:);
x2 = x(2,:);

y(1,:) = x2;
y(2,:) = (1+x2.^2)./x1;

end

function y = catn_dt(t, x, y, p)
%CATN   'coll'-compatible encoding of catenary vector field

y = zeros(2,numel(t));

end


function y = catn_dx(t, x, y, p)
%CATN   'coll'-compatible encoding of catenary vector field

x1 = x(1,:);
x2 = x(2,:);

y = zeros(2,2,numel(x1));
y(1,2,:) = 1;
y(2,1,:) = -(1+x2.^2)./x1.^2;
y(2,2,:) = 2*x2./x1;

end

function y = catn_dy(t, x, y, p)
%CATN   'coll'-compatible encoding of catenary vector field

y = [];

end


function y = catn_dp(t, x, y, p)

x1 = x(1,:);
y = zeros(2,1,numel(x1));

end


function y = catn_dxdx(x, p)
%CATN   'coll'-compatible encoding of catenary vector field

x1 = x(1,:);
x2 = x(2,:);

y = zeros(2,2,2,numel(x1));
y(2,1,1,:) = 2*(1+x2.^2)./x1.^3;
y(2,1,2,:) = -2*x2./x1.^2;
y(2,2,1,:) = -2*x2./x1.^2;
y(2,2,2,:) = 2./x1;

end



function y = catn_dxdp(x, p)

x1 = x(1,:);
x2 = x(2,:);

y = zeros(2,2,1,numel(x1));

end


function y = catn_dpdp(x, p)

x1 = x(1,:);
y = zeros(2,1,1,numel(x1));

end


function [data, y] = catn_bc(prob, data, u)

T  = u(1);
x0 = u(2:3);
x1 = u(4:5);
p  = u(6);
T0 = u(7);

y = [T-1; x0(1)-1; x1(1)-p; T0];
  
end


function [data, J] = catn_bc_du(prob, data, u)

T  = u(1);
x0 = u(2:3);
x1 = u(4:5);
p  = u(6);
T0 = u(7);

J = zeros(4,7);
J(1,1) = 1;
J(2,2) = 1;
J(3,4) = 1;
J(3,6) = -1;
J(4,7) = 1;

end


function [data, dJ] = catn_bc_dudu(prob, data, u)

dJ = zeros(3,6,6);

end


