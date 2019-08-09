coco_use_recipes_toolbox ddaecoll_v3 adjt_ddaecoll_v3

% Here we try to solve the problem presented in section 7.3.1 in Recipes.
% We rewrite the x1'=x2 to be x1'=y and y=x2. The Jacobian of
% vector field is computed via numerical differentiation.

t0 = [0; 1];
x0 = [1 0; 1 0];
y0 = [0; 0];
Y0 = 0.8;
prob = coco_prob();
% prob = coco_set(prob,'ddaecoll','Dpoints','Flipped_Gauss_Radau');
% prob = coco_set(prob,'ddaecoll','Dnodes','Flipped_Gauss_Radau');
% prob = coco_set(prob,'ddaecoll','Apoints','Flipped_Gauss_Radau');
% zero problem
prob = ddaecoll_isol2seg(prob, '', @catenary, @catenary_dt, @catenary_dx,...
 @catenary_dy, @catenary_dp_v2, t0, x0, y0, 'Y', Y0); % Build 'coll' continuation problem
[data, uidx] = coco_get_func_data(prob, 'ddaecoll', 'data', 'uidx');
prob = coco_add_func(prob, 'gfunc',@coup, data,'zero','uidx',[data.xbp_idx;data.ybp_idx]);
prob = coco_add_func(prob, 'bc', @catn_bc, @catn_bc_du, data, 'zero', 'uidx', ...
  uidx([data.T_idx; data.x0_idx; data.x1_idx; data.p_idx; data.T0_idx]));
prob = coco_add_pars(prob, 'obj', uidx(data.x1_idx(1)), 'x1_end');

% adjoints
prob = adjt_ddaecoll_isol2seg(prob, '');
[data, axidx] = coco_get_adjt_data(prob, 'ddaecoll', 'data', 'axidx');
opt = data.ddaecoll_opt;
prob = coco_add_adjt(prob, 'gfunc', @adj_coup, data, 'aidx', ...
  axidx([opt.xcn_idx; opt.ycn_idx]));
prob = coco_add_adjt(prob, 'bc', 'aidx', ...
  axidx([opt.T_idx; opt.x0_idx; opt.x1_idx; opt.p_idx; opt.T0_idx]));
prob = coco_add_adjt(prob, 'obj', 'd.x1_end', 'aidx', axidx(opt.x1_idx(1)));


bd1 = coco(prob, 'cant1', [], 1, {'x1_end' 'd.x1_end' 'Y'}, [0 1]);



%% branch switch from fold to grow nontrivial adjoint
BPlab = coco_bd_labs(bd1, 'BP');
% zero problems

prob = coco_prob();
prob = ddaecoll_sol2seg(prob, '', 'cant1', BPlab); % Reconstruct 'coll' continuation problem
data = coco_get_func_data(prob, 'ddaecoll', 'data');
prob = coco_add_func(prob, 'gfunc',@coup, data,'zero','uidx',[data.xbp_idx;data.ybp_idx]);
prob = coco_add_func(prob, 'bc', @catn_bc, @catn_bc_du, data, 'zero', 'uidx', ...
  uidx([data.T_idx; data.x0_idx; data.x1_idx; data.p_idx; data.T0_idx]));
prob = coco_add_pars(prob, 'obj', uidx(data.x1_idx(1)), 'x1_end');

% branch switch data
chart = coco_read_solution('cant1', BPlab, 'chart');
cdata = coco_get_chart_data(chart, 'lsol');
% adjoints
prob = adjt_ddaecoll_sol2seg(prob, '', 'cant1', BPlab, 'BP');
[data, axidx] = coco_get_adjt_data(prob, 'ddaecoll', 'data', 'axidx');
opt = data.ddaecoll_opt;
[chart, lidx] = coco_read_adjoint('gfunc', 'cant1', BPlab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'gfunc', @adj_coup, data, 'aidx', ...
  axidx([opt.xcn_idx; opt.ycn_idx]),'l0', chart.x, 'tl0', cdata.v(lidx));
[chart, lidx] = coco_read_adjoint('bc', 'cant1', BPlab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'bc', 'aidx', ...
  axidx([opt.T_idx; opt.x0_idx; opt.x1_idx; opt.p_idx; opt.T0_idx]),...
  'l0', chart.x, 'tl0', cdata.v(lidx));
[chart, lidx] = coco_read_adjoint('obj', 'cant1', BPlab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'obj', 'd.x1_end', 'aidx', axidx(opt.x1_idx(1)),...
    'l0', chart.x, 'tl0', cdata.v(lidx));


bd2 = coco(prob, 'cant2', [], 1, {'d.x1_end' 'x1_end' 'Y'}, [0 1]);



% % coco(prob, 'coll2', [], 1, {'y1e' 'y2s'}, [0 3]);
for lab = 1:11
    figure(2)
    [sol data] = ddaecoll_read_solution('', 'coll2', lab);
    plot(sol.t,sol.x(:,1)); pause(0.5);hold on
end

coco_use_recipes_toolbox % remove the coll_v1 toolbox from the search path
