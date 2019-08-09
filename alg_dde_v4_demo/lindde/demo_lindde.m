coco_use_recipes_toolbox ddaecoll_v3 alg_dde_v4

% here we try to solve the following problem
% x1'=x1(t-1), x2'=x1(t-1)+x2(t-0.2), x3'=x2(t)
% with history x1(t)=x2(t)=x3(t)=1 for t<=0.
% In the framework of our code, we have T0=0, tau1=1 and tau2=0.2
% y1(t)=x1(t-tau1) if t>T0+tau1 and y1(t)=1 if t<T0+tau1, likewise, we have
% y2(t)=x2(t-tau2) if t>T0+tau2 and y2(t)=1 if t<T0+tau2



t0 = [0; 2];
x0 = [1 1 1; 0.1 0.2 0.4];
y0 = [1 1; 0.1 0.4];

prob = coco_prob();
prob = coco_set(prob, 'cont', 'h_max', 10);
% prob = coco_set(prob, 'ddaecoll', 'NTST', 10);
prob = ddaecoll_isol2seg(prob, '', @lindde, t0, x0, y0, {'tau1', 'tau2'}, [1; 0.2]); % Build 'coll' continuation problem
% g1 = {@gfunc,(@tau^x | '[]'),(@tau^t | '[]'),(@phi |'[]'),dg,dtaux,dtaut}
dg1    = {@g1dt, @g1dy, @g1du, @g1dp};
dg12   = {@g12dt, @g12dy, @g12dx, @g12dp};
dg22   = {@g22dt, @g22dy, @g22dx, @g22dp};
dtaux1 = {@dtaux1dT0, @dtaux1dT, @dtaux1dp};
dtaux2 = {@dtaux2dT0, @dtaux2dT, @dtaux2dp};
dphi   = @Dphi;
g_11 = {@g1, [], @taux1, @phi, dg1, [], dtaux1, dphi};
g_12 = {@g12, @taux1, [], [], dg12, dtaux1};
g_21 = {@g1, [], @taux2, @phi, dg1, [], dtaux2, dphi};
g_22 = {@g22, @taux2, [], [], dg22, dtaux2};
prob = alg_dde_isol2seg(prob, '', 'g1', 1, {@s12}, {g_11, g_12});
prob = alg_dde_isol2seg(prob, '', 'g2', 2, {@s22}, {g_21, g_22});

data = coco_get_func_data(prob, 'ddaecoll', 'data'); % Extract toolbox data
% prob = coco_add_func(prob, 'gfunc',@coup, data,'zero','uidx',[data.xbp_idx;data.ybp_idx;data.T0_idx;data.T_idx]);
prob = coco_add_pars(prob, 'pars', ...
  [data.x0_idx; data.x1_idx(1); data.T0_idx; data.T_idx], ...
  {'y1s' 'y2s' 'y3S' 'y1e' 'T0' 'T'});
% coco(prob, 'coll1', [], 1, {'T' 'y1e'}, [1.9 2]);
coco(prob, 'coll1', [], 1, {'T' 'y1e'}, [1.5 2.5]);


% for lab = 2:4
%     [sol data] = ddaecoll_read_solution('', 'coll1', lab);
%     plot(sol.t,sol.x(:,1),'r.');hold on
%     plot(sol.t,sol.x(:,2),'b.');
%     plot(sol.t,sol.x(:,3),'k.'); 
% %     pause(2.0)
% end
% legend('x1,dde23','x2,dde23','x3,dde23','x1,coco','x2,coco','x3,coco');


prob = coco_prob();
prob = coco_set(prob, 'cont', 'h_max', 10);
% prob = coco_set(prob, 'cont', 'NPR', 1);
% prob = coco_set(prob, 'ddaecoll', 'NTST', 10);
prob = ddaecoll_sol2seg(prob, '', 'coll1', 4); % Reconstruct 'coll' continuation problem
prob = alg_dde_sol2seg(prob, '', 'g1', 'coll1', 4);
prob = alg_dde_sol2seg(prob, '', 'g2', 'coll1', 4);
data = coco_get_func_data(prob, 'ddaecoll', 'data');
% prob = coco_add_func(prob, 'gfunc',@coup, data,'zero','uidx',[data.xbp_idx;data.ycn_idx]);
prob = coco_add_pars(prob, 'pars', ...
  [data.x0_idx; data.x1_idx(2); data.T0_idx; data.T_idx], ...
  {'y1s' 'y2s' 'y3s' 'y2e' 'T0' 'T'});
coco(prob, 'coll2', [], 1, {'tau2' 'y2e'}, [0 0.2]);


% for lab = 7
%     figure(2)
%     [sol data] = ddaecoll_read_solution('', 'coll2', lab);
%     plot(sol.t,sol.x(:,1),'r.');hold on
%     plot(sol.t,sol.x(:,2),'b.');
%     plot(sol.t,sol.x(:,3),'k.'); 
% end


coco_use_recipes_toolbox % remove the coll_v1 toolbox from the search path
