coco_use_recipes_toolbox ddaecoll_v3

% here we try to solve the following problem
% x1'=x1(t-1), x2'=x1(t-1)+x2(t-0.2), x3'=x2(t)
% with history x1(t)=x2(t)=x3(t)=1 for t<=0.
% In the framework of our code, we have
% y1(t)=x1(t-1) and y2(t)=x2(t-0.2)



t0 = [0; 2];
x0 = [1 1 1; 0.1 0.2 0.4];
y0 = [1 1; 0.1 0.4];

prob = ddaecoll_isol2seg(coco_prob(), '', @lindde, t0, x0, y0, []); % Build 'coll' continuation problem
data = coco_get_func_data(prob, 'ddaecoll', 'data'); % Extract toolbox data
prob = coco_add_func(prob, 'gfunc',@coup, data,'zero','uidx',[data.xbp_idx;data.ybp_idx;data.T0_idx;data.T_idx]);
prob = coco_add_pars(prob, 'pars', ...
  [data.x0_idx; data.x1_idx(1); data.T0_idx; data.T_idx], ...
  {'y1s' 'y2s' 'Y3S' 'y1e' 'T0' 'T'});
coco(prob, 'coll1', [], 1, {'T' 'y1e'}, [0.1 2]);

for lab = 1:1
    [sol data] = ddaecoll_read_solution('', 'coll1', lab);
    plot(sol.t,sol.x(:,1),'r.');hold on
    plot(sol.t,sol.x(:,2),'b.');
    plot(sol.t,sol.x(:,3),'k.'); 
%     pause(2.0)
end
legend('x1,dde23','x2,dde23','x3,dde23','x1,coco','x2,coco','x3,coco');



% prob = ddecoll_sol2seg(coco_prob(), '', 'coll1', 5); % Reconstruct 'coll' continuation problem
% data = coco_get_func_data(prob, 'ddecoll', 'data');
% prob = coco_add_func(prob, 'gfunc',@coup, data,'zero','uidx',[data.xbp_idx;data.ycn_idx]);
% prob = coco_add_pars(prob, 'pars', ...
%   [data.x0_idx; data.x1_idx(1); data.T0_idx; data.T_idx], ...
%   {'y1s' 'y2s' 'y1e' 'T0' 'T'});
% coco(prob, 'coll2', [], 1, {'y1e' 'y2s'}, [0 3]);
% for lab = 1:11
%     figure(2)
%     [sol data] = ddecoll_read_solution('', 'coll2', lab);
%     plot(sol.t,sol.x(:,1)); pause(0.5);hold on
% end

% coco_use_recipes_toolbox % remove the coll_v1 toolbox from the search path
