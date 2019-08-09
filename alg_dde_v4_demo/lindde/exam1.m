function sol = exam1
% This is example 3 of D.R. Wille' and C.T.H. Baker,
% DELSOL--a numerical code for the solution of systems
% of delay-differential equations, Appl. Numer. Math., 
% 9 (992) 223-234. 

% Copyright 2002, The MathWorks, Inc.

sol = dde23(@exam1f,[1, 0.001],ones(3,1),[0, 2.5]);
figure
plot(sol.x,sol.y);
title('Example 3 of Wille'' and Baker.')
xlabel('time t');
ylabel('y(t)');

%-----------------------------------------------------------------------

function yp = exam1f(t,y,Z)
%EXAM1F  The derivative function for the Example 1 of the DDE Tutorial.
ylag1 = Z(:,1);
ylag2 = Z(:,2);
yp = [ ylag1(1)
       ylag1(1) + ylag2(2)
       y(2)                ];
