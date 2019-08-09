function ddae_monotonic()

coco_use_recipes_toolbox ddaecoll_v3 adjt_ddaecoll_v3 alg_dae_v3 adjt_alg_dae_v3


%% we reformulate the original problem such that x2=y
clear
prob = coco_prob;
prob = coco_set(prob, 'coll', 'TOL', 1e-1);
prob = coco_set(prob,'ddaecoll','Dpoints','Gauss_Radau');
% prob = coco_set(prob,'ddaecoll','Dnodes','Gauss_Radau');
prob = coco_set(prob,'ddaecoll','Apoints','Gauss_Radau');
% prob = coco_set(prob,'ddaecoll','ANodes','Gauss_Radau');
% prob = coco_set(prob,'coll','NTST',15);
% prob = coco_set(prob,'coll','NCOL',4);
% prob = coco_set(prob,'cont', 'ItMX', 550);
% prob = coco_set(prob, 'cont', 'NAdapt', 5);
p1 = [-0.1, 3.5];
p2 = [-0.1, 0.8];
p3 = [-0.1, 0.6];
objint = [0, 1.5];

%% first run to find initial fold
% zero problems
funcs = {@obv, @obv_dt, @obv_dx, @obv_dy, @obv_dp};
coll_args = [funcs, {[0; 1], [0 0;0 0], [0;0] {'l1', 'l2', 'l3'}, [0.1;0.1;0.1]}];
prob1 = ddaecoll_isol2seg(prob, '', coll_args{:});
[data, uidx] = coco_get_func_data(prob1, 'ddaecoll', 'data', 'uidx');
maps = data;
bc_funcs = {@obv_bc, @obv_bc_du, @obv_bc_dudu};
prob1 = coco_add_func(prob1, 'bc', bc_funcs{:}, [], 'zero', ...
  'uidx', uidx([maps.T_idx; maps.x0_idx; maps.x1_idx; maps.p_idx; maps.T0_idx]));
prob1 = alg_dae_isol2seg(prob1, '', @gfunc, @gfunc_dt, @gfunc_dx, @gfunc_dy, @gfunc_dp);
prob1 = coco_add_func(prob1, 'obj', @objhan, @objhan_du, @objhan_dudu, [], ...
  'inactive', 'obj', 'uidx', uidx(maps.p_idx));

% adjoints
prob1 = adjt_ddaecoll_isol2seg(prob1, '');

[data, axidx] = coco_get_adjt_data(prob1, 'ddaecoll', 'data', 'axidx');
opt = data.ddaecoll_opt;
prob1 = coco_add_adjt(prob1, 'bc', 'aidx', ...
  axidx([opt.T_idx; opt.x0_idx; opt.x1_idx; opt.p_idx; opt.T0_idx]));

prob1 = adjt_alg_dae_isol2seg(prob1, '');

prob1 = coco_add_adjt(prob1, 'obj', 'd.obj', 'aidx', axidx(opt.p_idx));

cont_pars =  {'obj', 'l1', 'd.obj', 'd.l2', 'd.l3', 'd.l1'};
coco(prob1, 'obv1', [], 1, cont_pars, {objint,p1});


%% branch switch from fold to grow nontrivial adjoint
bd1   = coco_bd_read('obv1');
BPlab = coco_bd_labs(bd1, 'BP');
lab   = BPlab(2);

% zero problems
prob2 = ddaecoll_sol2seg(prob, '', 'obv1', lab);
[data, uidx] = coco_get_func_data(prob2, 'ddaecoll', 'data', 'uidx');
maps = data;
prob2 = coco_add_func(prob2, 'bc', bc_funcs{:}, [], 'zero', ...
  'uidx', uidx([maps.T_idx; maps.x0_idx; maps.x1_idx; maps.p_idx; maps.T0_idx]));
prob2 = alg_dae_sol2seg(prob2, '', 'obv1', lab);
prob2 = coco_add_func(prob2, 'obj', @objhan, @objhan_du, @objhan_dudu, [], ...
  'inactive', 'obj', 'uidx', uidx(maps.p_idx));

% branch switch data
chart = coco_read_solution('obv1', lab, 'chart');
cdata = coco_get_chart_data(chart, 'lsol');

% adjoints
prob2 = adjt_ddaecoll_sol2seg(prob2, '', 'obv1', lab, 'BP');
[chart, aidx] = coco_read_adjoint('bc', 'obv1', lab, 'chart', 'lidx');
[data, axidx] = coco_get_adjt_data(prob2, 'ddaecoll', 'data', 'axidx');
opt = data.ddaecoll_opt;
prob2 = coco_add_adjt(prob2, 'bc', 'aidx', ...
  axidx([opt.T_idx; opt.x0_idx; opt.x1_idx; opt.p_idx; opt.T0_idx]), ...
  'l0', chart.x, 'tl0', cdata.v(aidx));
prob2 = adjt_alg_dae_sol2seg(prob2, '', 'obv1', lab, 'BP');
[chart, aidx] = coco_read_adjoint('obj', 'obv1', lab, 'chart', 'lidx');
prob2 = coco_add_adjt(prob2, 'obj','d.obj', 'aidx', axidx(opt.p_idx),...
  'l0', chart.x, 'tl0', cdata.v(aidx));

% computational domain
dobj_int = [chart.x(1) 1.1];
prob2 = coco_add_event(prob2, 'opt', 'BP', 'd.obj', '>', 1); 
cont_pars = {'d.obj', 'obj', 'l1', 'd.l2', 'd.l3', 'd.l1'};
coco(prob2, 'obv2', [], cont_pars, dobj_int);


%% continue to let d.l2=0
bd2 = coco_bd_read('obv2');
lab = coco_bd_labs(bd2, 'opt');

% zero problems
prob3 = ddaecoll_sol2seg(prob, '', 'obv2', lab);
[data, uidx] = coco_get_func_data(prob3, 'ddaecoll', 'data', 'uidx');
maps = data;
prob3 = coco_add_func(prob3, 'bc', bc_funcs{:}, [], 'zero', ...
  'uidx', uidx([maps.T_idx; maps.x0_idx; maps.x1_idx; maps.p_idx; maps.T0_idx]));
prob3 = alg_dae_sol2seg(prob3, '', 'obv2', lab);
prob3 = coco_add_func(prob3, 'obj', @objhan, @objhan_du, @objhan_dudu, [], ...
  'inactive', 'obj', 'uidx', uidx(maps.p_idx));

% adjoints
prob3 = adjt_ddaecoll_sol2seg(prob3, '', 'obv2', lab);
chart = coco_read_adjoint('bc', 'obv2', lab, 'chart');
[data, axidx] = coco_get_adjt_data(prob3, 'ddaecoll', 'data', 'axidx');
opt = data.ddaecoll_opt;
prob3 = coco_add_adjt(prob3, 'bc', 'aidx', ...
  axidx([opt.T_idx; opt.x0_idx; opt.x1_idx; opt.p_idx; opt.T0_idx]), ...
  'l0', chart.x);
prob3 = adjt_alg_dae_sol2seg(prob3, '', 'obv2', lab);
chart = coco_read_adjoint('obj', 'obv2', lab, 'chart');
prob3 = coco_add_adjt(prob3, 'obj','d.obj', 'aidx', axidx(opt.p_idx),...
  'l0', chart.x);

% % computational domain
chart = coco_read_adjoint('ddaecoll.pars', 'obv2', lab, 'chart');
dl20 = chart.x(2);
if dl20>0
  dl2_int = [-0.1 dl20];
  prob3 = coco_add_event(prob3, 'opt', 'BP', 'd.l2', '<', 0);
else
  dl2_int = [dl20 0.1];
  prob3 = coco_add_event(prob3, 'opt', 'BP', 'd.l2', '>', 0);
end

cont_pars = {'d.l2', 'l1', 'l2', 'obj', 'd.l3', 'd.l1'};
coco(prob3, 'obv3', [], cont_pars, {dl2_int});


%% continue to let d.l3=0
bd3 = coco_bd_read('obv3');
lab = coco_bd_labs(bd3, 'opt');

% zero problems
prob4 = ddaecoll_sol2seg(prob, '', 'obv3', lab);
[data, uidx] = coco_get_func_data(prob4, 'ddaecoll', 'data', 'uidx');
maps = data;
prob4 = coco_add_func(prob4, 'bc', bc_funcs{:}, [], 'zero', ...
  'uidx', uidx([maps.T_idx; maps.x0_idx; maps.x1_idx; maps.p_idx; maps.T0_idx]));
prob4 = alg_dae_sol2seg(prob4, '', 'obv3', lab);
prob4 = coco_add_func(prob4, 'obj', @objhan, @objhan_du, @objhan_dudu, [], ...
  'inactive', 'obj', 'uidx', uidx(maps.p_idx));

% adjoints
prob4 = adjt_ddaecoll_sol2seg(prob4, '', 'obv3', lab);
chart = coco_read_adjoint('bc', 'obv3', lab, 'chart');
[data, axidx] = coco_get_adjt_data(prob4, 'ddaecoll', 'data', 'axidx');
opt = data.ddaecoll_opt;
prob4 = coco_add_adjt(prob4, 'bc', 'aidx', ...
  axidx([opt.T_idx; opt.x0_idx; opt.x1_idx; opt.p_idx; opt.T0_idx]), ...
  'l0', chart.x);
prob4 = adjt_alg_dae_sol2seg(prob4, '', 'obv3', lab);
chart = coco_read_adjoint('obj', 'obv3', lab, 'chart');
prob4 = coco_add_adjt(prob4, 'obj','d.obj', 'aidx', axidx(opt.p_idx),...
  'l0', chart.x);

% % computational domain
chart = coco_read_adjoint('ddaecoll.pars', 'obv3', lab, 'chart');
dl30 = chart.x(3);
if dl30>0
  dl3_int = [-0.1 dl30];
  prob4 = coco_add_event(prob4, 'opt', 'BP', 'd.l3', '<', 0);
else
  dl3_int = [dl30 0.1];
  prob4 = coco_add_event(prob4, 'opt', 'BP', 'd.l3', '>', 0);
end

cont_pars = {'d.l3', 'obj', 'l1', 'l2', 'l3', 'd.coll.T0', 'd.l1',};
coco(prob4, 'obv4', [], cont_pars, {dl3_int});



%% plot results
% --------- %
figure(2)
bd = coco_bd_read('obv1');
l1 = coco_bd_col(bd,'l1');
l2 = coco_bd_col(bd,'l2');
l3 = coco_bd_col(bd,'l3');
obj = coco_bd_col(bd,'obj');
figure(2)
plot3(l1,l2,l3,'k-','LineWidth',2); hold on
idx = coco_bd_idxs(bd, 'BP');
plot3(l1(idx), l2(idx), l3(idx), 'ro', 'MarkerFaceColor', 'r')
figure(3)
plot3(l1,l2,obj,'k-','LineWidth',2); hold on
plot3(l1(idx), l2(idx), obj(idx), 'ro', 'MarkerFaceColor', 'r')
idx = coco_bd_idxs(bd, 'EP');
[~,idx] = max(obj(idx));
plot3(l1(idx), l2(idx), obj(idx), 'go', 'MarkerFaceColor', 'g')




    
bd = coco_bd_read('obv3');
l1 = coco_bd_col(bd,'l1');
l2 = coco_bd_col(bd,'l2');
l3 = coco_bd_col(bd,'l3');
obj = coco_bd_col(bd,'obj');
figure(2)
plot3(l1,l2,l3,'k-','LineWidth',2); hold on
idx = coco_bd_idxs(bd, 'opt');
plot3(l1(idx), l2(idx), l3(idx), 'bo', 'MarkerFaceColor', 'b')
figure(3)
plot3(l1,l2,obj,'k-','LineWidth',2); hold on
plot3(l1(idx), l2(idx), obj(idx), 'bo', 'MarkerFaceColor', 'b')
    

bd = coco_bd_read('obv4');
l1 = coco_bd_col(bd,'l1');
l2 = coco_bd_col(bd,'l2');
l3 = coco_bd_col(bd,'l3');
obj = coco_bd_col(bd,'obj');
figure(2)
plot3(l1,l2,l3,'k-','LineWidth',2); hold on
idx = coco_bd_idxs(bd, 'opt');
plot3(l1(idx), l2(idx), l3(idx), 'ko', 'MarkerFaceColor', 'k')
figure(3)
plot3(l1,l2,obj,'k-','LineWidth',2); hold on
plot3(l1(idx), l2(idx), obj(idx), 'ko', 'MarkerFaceColor', 'k')

set(gca,'LineWidth',1.2);
set(gca,'FontSize',14);
xlabel('$p_1$','interpreter','latex','FontSize',14);
ylabel('$p_2$','interpreter','latex','FontSize',14);
% zlabel('$p_3$','interpreter','latex','FontSize',14);
% axis([-0.1 3.5 -0.1 0.8 -0.1 0.6])
zlabel('$\mu_{obj}$','interpreter','latex','FontSize',14);
axis([-0.1 3.5 -0.1 0.8 0 1.5])
grid on 
box on

%%

end

function f = obv(t,x,y,p)

x1 = x(1,:);
p1 = p(1,:);
p2 = p(2,:);
p3 = p(3,:);

f(1,:) = y;
f(2,:) = -p1.*exp(x1+p2.*x1.^2+p3.*x1.^4);

end

function f = obv_dt(t,x,y,p)

f = zeros(2,numel(t));

end

function y = obv_dx(t,x,y,p)

x1 = x(1,:);
p1 = p(1,:);
p2 = p(2,:);
p3 = p(3,:);
pp = x1+p2.*x1.^2+p3.*x1.^4;
pd = 1+2*p2.*x1+4*p3.*x1.^3;

y = zeros(2,2,numel(x1));
y(2,1,:) = -p1.*exp(pp).*pd;

end


function f = obv_dy(t,x,y,p)

f = zeros(2,1,numel(t));
f(1,1,:) = 1;

end

function f = obv_dp(t,x,y,p)

x1 = x(1,:);
p1 = p(1,:);
p2 = p(2,:);
p3 = p(3,:);
pp = x1+p2.*x1.^2+p3.*x1.^4;
ep = exp(pp);

f = zeros(2,3,numel(x1));
f(2,1,:) = -ep;
f(2,2,:) = -p1.*ep.*x1.^2;
f(2,3,:) = -p1.*ep.*x1.^4;

end


function [data, y] = obv_bc(prob, data, u) %#ok<INUSL>

T  = u(1);
x0 = u(2:3);
x1 = u(4:5);
p  = u(6:8); %#ok<NASGU>
T0 = u(9);

y = [T-1; x0(1); x1(1); T0];

end

function [data, J] = obv_bc_du(prob, data, u) %#ok<INUSD,INUSL>

J = zeros(4,9);
J(1,1) = 1;
J(2,2) = 1;
J(3,4) = 1;
J(4,9) = 1;

end

function [data, dJ] = obv_bc_dudu(prob, data, u) %#ok<INUSD,INUSL>

dJ = zeros(4,9,9);

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