coco_use_recipes_toolbox ddaecoll_v3 alg_ddae_v3 % add the coll_v1 toolbox to the search path

% We solve initial value problem of modified DDAEs in model 1 of paper
% 'Periodic solutions of differential algebraic equations with time delays'
% x1(t)' =
% x1(t)^2-0.5*x1(t-tau1)*x2(t)+alpha*sin(x1(t-tau))-8*x1(t)*x2(t-tau2)+1
% with /*modified*/ algebraic equations
% -2*x1(t)^3+x2(t)=0 with history function x1(t)=1 when t<=0
% In our framework, we have x=x1, y1=x1(t-tau1), y2=x2, y3=y2(t-tau2)

% For the sake of comparison, we plug x2=2*x1^3 into the governing equation
% for x1 such that a standard initial-value problem of dde is obtained.
% This IVP can be solved using dde23


alpha = 2.2;
tau1  = 0.5;
tau2  = 0.2;

sol = ref_dde23(1.0,tau1,tau2,alpha);
t0  = sol.x';
x0  = sol.y';
ts  = [-1;-0.6;-0.3;t0];
xs  = [1;1;1;x0];

y10 = interp1(ts,xs,t0-tau1);
y20 = 2*x0.^3;
y2s = [2;2;2;y20];
y30 = interp1(ts,y2s,t0-tau2);
y0  = [y10 y20 y30];


pname = {'alpha', 'tau1', 'tau2'};
p0    = [alpha tau1 tau2]';
%% continuation in alpha
prob = coco_prob();
prob = coco_set(prob, 'cont', 'h_max', 10);
% prob = coco_set(prob, 'ddaecoll', 'NTST', 10);
prob = ddaecoll_isol2seg(prob, '', @model1, t0, x0, y0, pname, p0); % Build 'coll' continuation problem
% g1 = {@gfunc,(@tau^x | '[]'),(@tau^y | '[]'),(@tau^t | '[]'),(@phi |'[]'),dg,dtaux,dtauy,dtaut}
g_11 = {@g11, [], [], @taux1, @phi1};
g_12 = {@g12, @taux1};
g_2  = {@g2};
g_31 = {@g31, [], [], @taux2, @phi2};
g_32 = {@g32, [], @taux2};
prob = alg_ddae_isol2seg(prob, '', 'g1', 1, {@s1}, {g_11, g_12});
prob = alg_ddae_isol2seg(prob, '', 'g2', 1, {}, {{@g2}});
prob = alg_ddae_isol2seg(prob, '', 'g3', 1, {@s3}, {g_31, g_32});

data = coco_get_func_data(prob, 'ddaecoll', 'data'); % Extract toolbox data
prob = coco_add_pars(prob, 'pars', ...
  [data.x0_idx; data.T0_idx; data.T_idx], ...
  {'y1s' 'T0' 'T'});
coco(prob, 'coll1', [], 1, {'alpha' 'T'}, {[1.0 3.0],[]});

lab = 5;
[sol data] = ddaecoll_read_solution('', 'coll1', lab);
figure(1)
plot(sol.t,sol.x(:,1),'ro'); hold on
sol = ref_dde23(1,tau1,tau2,3);
plot(sol.x,sol.y);

%% continuation in tau1
prob = coco_prob();
prob = coco_set(prob, 'cont', 'h_max', 10);
% prob = coco_set(prob, 'cont', 'NPR', 1);
% prob = coco_set(prob, 'ddaecoll', 'NTST', 10);
prob = ddaecoll_sol2seg(prob, '', 'coll1', 5); % Reconstruct 'coll' continuation problem
prob = alg_ddae_sol2seg(prob, '', 'g1', 'coll1', 5);
prob = alg_ddae_sol2seg(prob, '', 'g2', 'coll1', 5);
prob = alg_ddae_sol2seg(prob, '', 'g3', 'coll1', 5);
data = coco_get_func_data(prob, 'ddaecoll', 'data');
prob = coco_add_pars(prob, 'pars', ...
  [data.x0_idx; data.T0_idx; data.T_idx], ...
  {'y1s' 'T0' 'T'});
coco(prob, 'coll2', [], 1, {'tau1' 'tau2'}, [0.4 0.6]);

lab = 6;
[sol data] = ddaecoll_read_solution('', 'coll2', lab);
figure(1)
plot(sol.t,sol.x(:,1),'bo'); hold on
sol = ref_dde23(1,0.4,tau2,3);
plot(sol.x,sol.y);

%% continuation in tau2
prob = coco_prob();
prob = coco_set(prob, 'cont', 'h_max', 10);
% prob = coco_set(prob, 'cont', 'NPR', 1);
% prob = coco_set(prob, 'ddaecoll', 'NTST', 10);
prob = ddaecoll_sol2seg(prob, '', 'coll2', 6); % Reconstruct 'coll' continuation problem
prob = alg_ddae_sol2seg(prob, '', 'g1', 'coll2', 6);
prob = alg_ddae_sol2seg(prob, '', 'g2', 'coll2', 6);
prob = alg_ddae_sol2seg(prob, '', 'g3', 'coll2', 6);
data = coco_get_func_data(prob, 'ddaecoll', 'data');
prob = coco_add_pars(prob, 'pars', ...
  [data.x0_idx; data.T0_idx; data.T_idx], ...
  {'y1s' 'T0' 'T'});
coco(prob, 'coll3', [], 1, {'tau2' 'tau1'}, [0.0 0.3]);

lab = 6;
[sol data] = ddaecoll_read_solution('', 'coll3', lab);
figure(1)
plot(sol.t,sol.x(:,1),'bo'); hold on
sol = ref_dde23(1,0.4,tau2,3);
plot(sol.x,sol.y);

coco_use_recipes_toolbox % remove the coll_v1 toolbox from the search path

