function sol = ref_dde23(T,tau1,tau2,alpha)
% This is example 3 of D.R. Wille' and C.T.H. Baker,
% DELSOL--a numerical code for the solution of systems
% of delay-differential equations, Appl. Numer. Math., 
% 9 (992) 223-234. 

% Copyright 2002, The MathWorks, Inc.

sol = dde23(@exam1f,[tau1, tau2],ones(1,1),[0, T],[],alpha);
% figure
% plot(sol.x,sol.y);
% title('Example 3 of Wille'' and Baker.')
% xlabel('time t');
% ylabel('y(t)');

%-----------------------------------------------------------------------

function yp = exam1f(t,y,Z,alpha)
%EXAM1F  The derivative function for the Example 1 of the DDE Tutorial.
ylag1 = Z(:,1);
ylag2 = Z(:,2);
x2 = 2*y^3;
y2 = 2*ylag2^3;
yp = y^2-0.5*ylag1*x2+alpha*sin(ylag1)-8*y*y2+1;
% yp = [ ylag1(1)
%        ylag1(1) + ylag2(2)
%        y(2)                ];
