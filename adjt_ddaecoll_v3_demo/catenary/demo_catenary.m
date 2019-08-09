coco_use_recipes_toolbox ddaecoll_v3 adjt_ddaecoll_v3

% Here we try to solve the problem presented in section 7.3.1 in Recipes.
% We rewrite the x1'=x2 to be x1'=y and y=x2. The Jacobian of
% vector field is computed via numerical differentiation.

t0 = [0; 1];
x0 = [1 0; 1 0.04];
y0 = [0; 0.04];
prob = coco_prob();
prob = coco_set(prob, 'cont', 'h_max', 100);
% prob = coco_set(prob,'ddaecoll','Dpoints','Flipped_Gauss_Radau');
% prob = coco_set(prob,'ddaecoll','Dnodes','Flipped_Gauss_Radau');
% prob = coco_set(prob,'ddaecoll','Apoints','Flipped_Gauss_Radau');
prob = ddaecoll_isol2seg(prob, '', @catenary, @catenary_dt, @catenary_dx,...
 @catenary_dy, @catenary_dp, t0, x0, y0, []); % Build 'coll' continuation problem
data = coco_get_func_data(prob, 'ddaecoll', 'data'); % Extract toolbox data
prob = coco_add_func(prob, 'gfunc',@coup, data,'zero','uidx',[data.xbp_idx;data.ybp_idx]);
prob = coco_add_pars(prob, 'pars', ...
  [data.x0_idx; data.T0_idx; data.T_idx], ...
  {'x1s' 'x2s' 'T0' 'T'});
prob = coco_add_pars(prob, 'obj', data.x1_idx(1), 'x1_end');

% adjoints
prob = adjt_ddaecoll_isol2seg(prob, '');
[data, axidx] = coco_get_adjt_data(prob, 'ddaecoll', 'data', 'axidx');
opt = data.ddaecoll_opt;
% prob = coco_add_adjt(prob, tbid, @adj, @adj_DU, data, 'l0', sol.l0, ...
%   'tl0', sol.tl0, 'adim', opt.adim);
% prob1 = coco_add_adjt(prob1, 'obj', @adj_objhan, @adj_objhan_du, data, ...
%   'd.obj', 'aidx', axidx([opt.xcn_idx; opt.T0_idx; opt.T_idx; opt.p_idx]), ...
%   'remesh', @adj_obj_remesh);
prob = coco_add_adjt(prob, 'gfunc', @adj_coup, data, 'aidx', ...
  axidx([opt.xcn_idx; opt.ycn_idx]));
prob = coco_add_adjt(prob, 'pars', {'d.x1s' 'd.x2s' 'd.T0' 'd.T'}, 'aidx', axidx([opt.x0_idx; opt.T0_idx; opt.T_idx]));
prob = coco_add_adjt(prob, 'obj', 'd.x1_end', 'aidx', axidx(opt.x1_idx(1)));


coco(prob, 'coll1', [], 1, {'x1_end' 'd.x1_end' 'd.x1s' 'x2s' 'd.T0' 'd.T'}, [0 1.2]);

for lab = 1:5
    [sol data] = ddaecoll_read_solution('', 'coll1', lab);
    plot(sol.t,sol.x(:,1)); pause(0.5);hold on
end
plot(sol.t,cosh(sol.t),'ro')
figure(10)
plot(sol.t,sinh(sol.t),'ro'); hold on
plot(sol.ta,sol.ya,'b.');

prob = coco_prob();
prob = coco_set(prob,'ddaecoll','Dpoints','Flipped_Gauss_Radau');
prob = coco_set(prob,'ddaecoll','Dnodes','Flipped_Gauss_Radau');
prob = coco_set(prob,'ddaecoll','Apoints','Flipped_Gauss_Radau');
prob = ddaecoll_sol2seg(prob, '', 'coll1', 5); % Reconstruct 'coll' continuation problem
data = coco_get_func_data(prob, 'ddaecoll', 'data');
prob = coco_add_func(prob, 'gfunc',@coup, data,'zero','uidx',[data.xbp_idx;data.ybp_idx]);
prob = coco_add_pars(prob, 'pars', ...
  [data.x0_idx; data.x1_idx(1); data.T0_idx; data.T_idx], ...
  {'y1s' 'y2s' 'y1e' 'T0' 'T'});
coco(prob, 'coll2', [], 1, {'y1e' 'y2s'}, [0 3]);
for lab = 1:11
    figure(2)
    [sol data] = ddaecoll_read_solution('', 'coll2', lab);
    plot(sol.t,sol.x(:,1)); pause(0.5);hold on
end

coco_use_recipes_toolbox % remove the coll_v1 toolbox from the search path
