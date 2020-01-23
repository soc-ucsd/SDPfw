%```matlab
%syms x1 x2;
%f=x1^4+x2^4-x1*x2;
%n=2;
%d=4;

%syms x1 x2 x3;
x1 = sdpvar(1);
x2 = sdpvar(1);
%x3 = sdpvar(1);
%f =  1 + x1^4 + x2^4 + x3^4 + x1*x2*x3 + x2;

n = 3;

%obj = -lower;
%F = (f - lower >= 0)
F = [];
obj = f;

opts_mom = sdpsettings('solver', 'moment');

d=2;
[sol,x,momentdata] = solvemoment(F, obj, [], d);
obj2 = value(obj);
[sol3,x3,momentdata3] = solvemoment(F, obj, [], d+1);
obj3 = value(obj);


% d=2;
% solver = 'sedumi';
% [opt,data,status]=blockpop_uncons_first(f,n,d,solver);

%p = 4;
%[opt_high,data_high,status_high]=blockpop_uncons_higher(n,data,solver);
%By default, a monomial basis computed by the Newton polytope method will be used. If we add 'newton',0 to the input,

