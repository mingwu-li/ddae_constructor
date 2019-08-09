function simplifed_dodel_coll_demo_optimization()

prob = coco_prob;
prob = coco_set(prob, 'coll', 'TOL', 1e-1);
% prob = coco_set(prob, 'coll', 'NTST', 25);
% prob = coco_set(prob, 'cont', 'NAdapt', 1);

%% first run to find initial fold
% zero problems
t0 = [0; 1];
x0 = [0 0;0 0];

%%
coll_args = {@obv, @obv_dx, @obv_dp, ...
  @obv_dxdx, @obv_dxdp, @obv_dpdp, t0, x0, ...
  {'l1', 'l2', 'l3'}, [0.1;0.1;0.1]};
prob1 = ode_isol2coll(prob, '', coll_args{:});
[data, uidx] = coco_get_func_data(prob1, 'coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob1 = coco_add_func(prob1, 'bc', @obv_bc, @obv_bc_du, ...
  @obv_bc_dudu, data, 'zero', 'uidx', ...
  uidx([maps.T_idx; maps.x0_idx; maps.x1_idx; maps.p_idx]));
prob1 = coco_add_func(prob1, 'obj', @objhan, @objhan_du, @objhan_dudu...
    , data, 'inactive', 'obj', 'uidx', uidx(maps.p_idx));
% adjoints
prob1 = adjt_isol2coll(prob1, '');
[data, axidx] = coco_get_adjt_data(prob1, 'coll', 'data', 'axidx');
opt = data.coll_opt;
prob1 = coco_add_adjt(prob1, 'bc', 'aidx', ...
  axidx([opt.T_idx; opt.x0_idx; opt.x1_idx; opt.p_idx]));

prob1 = coco_add_adjt(prob1, 'obj', 'd.obj', 'aidx', axidx(opt.p_idx));

bd1 = coco(prob1, 'obv1', [], 1, {'obj', 'd.obj', 'd.coll.T0', 'd.l2', 'd.l3', 'l1', 'd.l1', 'l2', 'l3'}, [0 2]);

%% branch switch from fold to grow nontrivial adjoint
BPlab = coco_bd_labs(bd1, 'BP');
% zero problems
prob2 = ode_BP2coll(prob, '', 'obv1', BPlab(2));    
[data, uidx] = coco_get_func_data(prob2, 'coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob2 = coco_add_func(prob2, 'bc', @obv_bc, @obv_bc_du, ...
  @obv_bc_dudu, data, 'zero', 'uidx', ...
  uidx([maps.T_idx; maps.x0_idx; maps.x1_idx; maps.p_idx]));    
prob2 = coco_add_func(prob2, 'obj', @objhan, @objhan_du, @objhan_dudu...
    , data, 'inactive', 'obj', 'uidx', uidx(maps.p_idx));
% branch switch data
chart = coco_read_solution('obv1', BPlab(2), 'chart');
cdata = coco_get_chart_data(chart, 'lsol');
% adjoints
prob2 = adjt_BP2coll(prob2, '', 'obv1', BPlab(2));
[chart, aidx] = coco_read_adjoint('bc', 'obv1', BPlab(2), 'chart', 'lidx');
[data, axidx] = coco_get_adjt_data(prob2, 'coll', 'data', 'axidx');
opt = data.coll_opt;
prob2   = coco_add_adjt(prob2, 'bc', 'aidx', ...
  axidx([opt.T_idx; opt.x0_idx; opt.x1_idx; opt.p_idx]), ...
  'l0', chart.x, 'tl0', cdata.v(aidx));
[chart, aidx] = coco_read_adjoint('obj', 'obv1', BPlab(2), 'chart', 'lidx');
prob2 = coco_add_adjt(prob2, 'obj', 'd.obj', 'aidx', ...
  axidx(opt.p_idx),...
  'l0', chart.x, 'tl0', cdata.v(aidx));

prob2 = coco_add_event(prob2, 'opt', 'BP', 'd.obj', '>', 1);
bd2 = coco(prob2, 'obv2', [], {'obj', 'd.obj', 'd.coll.T0', 'd.l2', 'd.l3', 'l1', 'd.l1', 'l2', 'l3'}, {[0 3], [0 1.1]});

%% continue to let d.p2=0
lab = coco_bd_labs(bd2, 'opt');
% zero problems
prob3 = ode_coll2coll(prob, '', 'obv2', lab(1));    
[data, uidx] = coco_get_func_data(prob3, 'coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob3 = coco_add_func(prob3, 'bc', @obv_bc, @obv_bc_du, ...
  @obv_bc_dudu, data, 'zero', 'uidx', ...
  uidx([maps.T_idx; maps.x0_idx; maps.x1_idx; maps.p_idx]));    
prob3 = coco_add_func(prob3, 'obj', @objhan, @objhan_du, @objhan_dudu...
    , data, 'inactive', 'obj', 'uidx', uidx(maps.p_idx));
% adjoints
prob3 = adjt_BP2coll(prob3, '', 'obv2', lab(1));
chart = coco_read_adjoint('bc', 'obv2', lab(1), 'chart');
[data, axidx] = coco_get_adjt_data(prob3, 'coll', 'data', 'axidx');
opt = data.coll_opt;
prob3   = coco_add_adjt(prob3, 'bc', 'aidx', ...
  axidx([opt.T_idx; opt.x0_idx; opt.x1_idx; opt.p_idx]), ...
  'l0', chart.x);
chart = coco_read_adjoint('obj', 'obv2', lab(1), 'chart');
prob3 = coco_add_adjt(prob3, 'obj', 'd.obj', 'aidx', ...
  axidx(opt.p_idx),...
  'l0', chart.x);

prob3 = coco_add_event(prob3, 'opt', 'BP', 'd.l2', '=', 0);
bd3 = coco(prob3, 'obv3', [], {'obj', 'l2', 'd.coll.T0', 'd.l2', 'd.l3', 'l1', 'd.l1', 'd.obj', 'l3'}, [0 3]);

%% continue to let d.p3=0
lab = coco_bd_labs(bd3, 'opt');
% zero problems
prob4 = ode_coll2coll(prob, '', 'obv3', lab(1));    
[data, uidx] = coco_get_func_data(prob4, 'coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob4 = coco_add_func(prob4, 'bc', @obv_bc, @obv_bc_du, ...
  @obv_bc_dudu, data, 'zero', 'uidx', ...
  uidx([maps.T_idx; maps.x0_idx; maps.x1_idx; maps.p_idx]));    
prob4 = coco_add_func(prob4, 'obj', @objhan, @objhan_du, @objhan_dudu...
    , data, 'inactive', 'obj', 'uidx', uidx(maps.p_idx));
% adjoints
prob4 = adjt_BP2coll(prob4, '', 'obv3', lab(1));
chart = coco_read_adjoint('bc', 'obv3', lab(1), 'chart');
[data, axidx] = coco_get_adjt_data(prob4, 'coll', 'data', 'axidx');
opt = data.coll_opt;
prob4   = coco_add_adjt(prob4, 'bc', 'aidx', ...
  axidx([opt.T_idx; opt.x0_idx; opt.x1_idx; opt.p_idx]), ...
  'l0', chart.x);
chart = coco_read_adjoint('obj', 'obv3', lab(1), 'chart');
prob4 = coco_add_adjt(prob4, 'obj', 'd.obj', 'aidx', ...
  axidx(opt.p_idx), 'l0', chart.x);

prob4 = coco_add_event(prob4, 'opt', 'BP', 'd.l3', '=', 0);
bd4 = coco(prob4, 'obv4', [], {'obj', 'l2', 'd.coll.T0', 'l3', 'd.l3', 'l1', 'd.l1', 'd.obj', 'd.l2'}, [0 3]);

end



function y=obv(x,p)

x1 = x(1,:);
x2 = x(2,:);
p1 = p(1,:);
p2 = p(2,:);
p3 = p(3,:); % p_bar


y(1,:) = x2;
y(2,:) = -p1.*exp(x1+p2.*x1.^2+p3.*x1.^4);

end


function y=obv_dx(x,p)

x1 = x(1,:);
x2 = x(2,:);
p1 = p(1,:);
p2 = p(2,:);
p3 = p(3,:); % p_bar

y = zeros(2,2,numel(x1));
y(1,2,:)=1;
pp=x1+p2.*x1.^2+p3.*x1.^4;
y(2,1,:)=-p1.*exp(pp).*(1+2*p2.*x1+4*p3.*x1.^3);

end


function y=obv_dp(x,p)

x1 = x(1,:);
x2 = x(2,:);
p1 = p(1,:);
p2 = p(2,:);
p3 = p(3,:); % p_bar

y = zeros(2,3,numel(x1));
p=x1+p2.*x1.^2+p3.*x1.^4;
y(2,1,:)=-exp(p);
y(2,2,:)=-p1.*exp(p).*x1.^2;
y(2,3,:)=-p1.*exp(p).*x1.^4;

end


function y=obv_dxdx(x,p)

x1 = x(1,:);
x2 = x(2,:);
p1 = p(1,:);
p2 = p(2,:);
p3 = p(3,:); % p_bar

y = zeros(2,2,2,numel(x1));
pp=x1+p2.*x1.^2+p3.*x1.^4;
y(2,1,1,:)=-p1.*exp(pp).*(1+2*p2.*x1+4*p3.*x1.^3).^2-p1.*exp(pp).*(2*p2+12*p3.*x1.^2);

end


function y=obv_dxdp(x,p)

x1 = x(1,:);
x2 = x(2,:);
p1 = p(1,:);
p2 = p(2,:);
p3 = p(3,:); % p_bar

y = zeros(2,2,3,numel(x1));
pp=x1+p2.*x1.^2+p3.*x1.^4;
y(2,1,1,:) = -exp(pp).*(1+2*p2.*x1+4*p3.*x1.^3);
y(2,1,2,:) = -p1.*exp(pp).*(1+2*p2.*x1+4*p3.*x1.^3).*x1.^2-2*p1.*exp(pp).*x1;
y(2,1,3,:) = -p1.*exp(pp).*(1+2*p2.*x1+4*p3.*x1.^3).*x1.^4-4*p1.*exp(pp).*x1.^3;

end


function y=obv_dpdp(x,p)

x1 = x(1,:);
x2 = x(2,:);
p1 = p(1,:);
p2 = p(2,:);
p3 = p(3,:); % p_bar

y = zeros(2,3,3,numel(x1));
p=x1+p2.*x1.^2+p3.*x1.^4;
y(2,1,2,:)=-exp(p).*x1.^2;
y(2,1,3,:)=-exp(p).*x1.^4;

y(2,2,1,:)=-exp(p).*x1.^2;
y(2,2,2,:)=-p1.*exp(p).*x1.^4;
y(2,2,3,:)=-p1.*exp(p).*x1.^6;

y(2,3,1,:)=-exp(p).*x1.^4;
y(2,3,2,:)=-p1.*exp(p).*x1.^6;
y(2,3,3,:)=-p1.*exp(p).*x1.^8;

end


function [data, y] = obv_bc(prob, data, u)

T  = u(1);
x0 = u(2:3);
x1 = u(4:5);
p  = u(6:8);

y = [T-1; x0(1); x1(1)];
  
end


function [data, J] = obv_bc_du(prob, data, u)

J = zeros(3,8);
J(1,1) = 1;
J(2,2) = 1;
J(3,4) = 1;

end


function [data, dJ] = obv_bc_dudu(prob, data, u)

dJ = zeros(3,8,8);

end


function [data, y] = objhan(prob, data, u)


p = u;    % Extract problem parameters
y = (p(1).^2+p(2).^2+p(3).^2)/10;

end


function [data, J] = objhan_du(prob, data, u)

p = u;
J = zeros(1,3);
J(1,1) = p(1)/5;
J(1,2) = p(2)/5;
J(1,3) = p(3)/5;

end

function [data, dJ] = objhan_dudu(prob, data, u)

dJ = zeros(1,3,3);
dJ(1,1,1) = 1/5;
dJ(1,2,2) = 1/5;
dJ(1,3,3) = 1/5;

end