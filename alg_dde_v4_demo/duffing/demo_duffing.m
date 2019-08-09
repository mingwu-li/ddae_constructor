coco_use_recipes_toolbox ddaecoll_v3 alg_dde_v4

% here we try to solve the periodic orbits of the following problem
% ddot{z}(t)+2*zeta*dot{z}(t)+z(t)+mu*z(t)^3=2*a*z(t-alpha)+2*b*dot{z}(t-alpha)
% +fcos(Omega*t+phi);
% let z1(t)=z(t), z2(t)=dot{x}(t), y=[z(t-alpha);dot{z}(t-\alpha)], it follows that the
% problem is formulated as dot{z1}(t)=z2(t), 
% dot{z2}(t) =
% -2*zeta*z2(t)-z1(t)-mu*z1(t)^3+2*a*y1+2*b*y2+fcos(Omega*t+phi), with bcs
% z(T0)-z(T0+T)=0 and T0=0, T=2pi/Omega, z2(0)=0.
% In our framework, we have y(t)=z(t+T-alpha) if T0<=t<=T0+alpha, and 
% y(t)=z(t-tau) if T0+alpha<t<T0+T.
% As an initial solution, we set alpha=mu=0 to reduce the original
% system back to linear odes, and the analytical solution can be obtained
% z(t)=f/[(1-2*a-Omega^2)^2+(2*(zeta-b)*Omega)^2]^0.5*cos(Omega*t), and the
% initial phi is determined by tan(phi)=[2*(zeta-b)*Omega]/[1-2*a-Omega^2]

zeta  = 0.05;
mu    = 0;
a     = 0.05;
b     = -0.05;
alpha = 0;
f     = 0.5;
Omega = 0.5;
phi   = atan(2*(zeta-b)*Omega/(1-2*a-Omega^2));
amp   = f/((1-2*a-Omega^2)^2+(2*(zeta-b)*Omega)^2)^0.5;
t0    = linspace(0,2*pi/Omega,51)';
z10   = amp*cos(Omega*t0);
z20   = -amp*Omega*sin(Omega*t0);
y10   = z10;
y20   = z20;
x0    = [z10 z20];
y0    = [y10 y20];
pname = {'zeta', 'mu', 'a', 'b', 'alpha', 'f', 'Omega', 'phi'};
p0    = [zeta mu a b alpha f Omega phi]';

prob = coco_prob();
prob = coco_set(prob, 'cont', 'h_max', 10);
% prob = coco_set(prob, 'ddaecoll', 'NTST', 10);
prob = ddaecoll_isol2seg(prob, '', @duffing, @duffing_dt, @duffing_dx, @duffing_dy, @duffing_dp, t0, x0, y0, pname, p0); % Build 'coll' continuation problem
% g1 = {@gfunc,(@tau^x | '[]'),(@tau^t | '[]'),(@phi |'[]'),dg,dtaux,dtaut}
dg1    = {@g1dt, @g1dy, @g1dx, @g1dp};
dtaux1 = {@dtaux1dT0, @dtaux1dT, @dtaux1dp};
dtaux2 = {@dtaux2dT0, @dtaux2dT, @dtaux2dp};
g_11 = {@g1, @taux1, [], [], dg1, dtaux1};
g_12 = {@g1, @taux2, [], [], dg1, dtaux2};
prob = alg_dde_isol2seg(prob, '', 'g1', (1:2), {@s2}, {g_11, g_12});
data = coco_get_func_data(prob, 'ddaecoll', 'data'); % Extract toolbox data
prob = coco_add_func(prob, 'bcs', @bc, [], 'zero', 'uidx',...
    [data.x0_idx; data.x1_idx; data.T0_idx; data.T_idx; data.p_idx(7)]);
prob = coco_add_pars(prob, 'pars', data.x0_idx(1), 'amp');
coco(prob, 'coll1', [], 1, {'amp' 'phi' 'mu'}, {[],[],[0 0.05]});


%% continuation along alpha
prob = coco_prob();
prob = coco_set(prob, 'cont', 'h_max', 10);
% prob = coco_set(prob, 'ddaecoll', 'NTST', 10);
prob = ddaecoll_sol2seg(prob, '', 'coll1', 2); % Reconstruct 'coll' continuation problem
prob = alg_dde_sol2seg(prob, '', 'g1', 'coll1', 2);
data = coco_get_func_data(prob, 'ddaecoll', 'data');
prob = coco_add_func(prob, 'bcs', @bc, [], 'zero', 'uidx',...
    [data.x0_idx; data.x1_idx; data.T0_idx; data.T_idx; data.p_idx(7)]);
prob = coco_add_pars(prob, 'pars', data.x0_idx(1), 'amp');

coco(prob, 'coll2', [], 1, {'amp' 'phi' 'alpha'}, {[],[],[0 1.25]});

figure(2)
coco_plot_bd('coll2', 'alpha', 'amp');

%% continuation along omega
profile on
prob = coco_prob();
prob = coco_set(prob, 'cont', 'h_max', 10);
prob = coco_set(prob, 'cont', 'ItMX', 150);
% prob = coco_set(prob, 'ddaecoll', 'NTST', 10);
prob = ddaecoll_sol2seg(prob, '', 'coll2', 2); % Reconstruct 'coll' continuation problem
prob = alg_dde_sol2seg(prob, '', 'g1', 'coll2', 2);
data = coco_get_func_data(prob, 'ddaecoll', 'data');
prob = coco_add_func(prob, 'bcs', @bc, [], 'zero', 'uidx',...
    [data.x0_idx; data.x1_idx; data.T0_idx; data.T_idx; data.p_idx(7)]);
prob = coco_add_pars(prob, 'pars', data.x0_idx(1), 'amp');

coco(prob, 'coll3', [], 1, {'amp' 'phi' 'Omega'}, {[],[],[0.5 2.0]});
profile viewer
figure(3)
coco_plot_bd('coll3', 'Omega', 'amp');

coco_use_recipes_toolbox % remove the coll_v1 toolbox from the search path
