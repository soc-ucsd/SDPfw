syms x1 x2 a b;

vartable = [x1, x2];

prog = sosprogram(vartable);

prog = sosdecvar(prog, [a, b]);

th = pi/4;
c = [cos(th); sin(th)];
%obj = cos(th)*a + sin(th)*b;
 obj = 2*a + 3*b;


%bivariate quartic
%f = x1^4 + x2^4 + a*x1^3*x2 + (2-a-b )*x1^2*x2^2 / 2 + 2*b*x1*x2^3;
f = x1^4 + (1+a)*x2^4  + (2-a-b )*x1^2*x2^2 / 2 + 2*b*x1*x2^3;

prog = sosineq(prog, f);

prog = sossetobj(prog, obj);

%call solver
solver_opt.solver = 'sedumi';
prog = sossolve(prog,solver_opt); 
% =============================================
% Finally, get solution
SOLy = sosgetsol(prog,[a,b])