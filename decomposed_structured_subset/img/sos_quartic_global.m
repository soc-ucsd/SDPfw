%mpol x(1) x(2) x(3) a b
a = sdpvar(1,1);
b = sdpvar(1,1);
x = sdpvar(3,1);


%f = x(1)^4 + (1+a)*x(2)^4 + (2+b)*x3^4 + (2-a-b)*x(1)^2*x(2)^2 / 2 + 2*b*x(1)*x(2)^3 ...
%    + (a-b)*x(2)^2*x3^2 + 3*a*x(1)*x3^3 + 2*x(1)^2*x(2)*x3;

f = x(1)^4 + (1+a)*x(2)^4 + (2+b)*x(3)^4 + (2-a-b)*x(1)^2*x(2)^2 / 2 + 2*b*x(1)*x(2)^3  + ...
   (a-b)*x(2)^2*x(3)^2 + 3*a*x(1)*x(3)^3 + 2*x(1)^2*x(2)*x(3);


%f = x(1)^4 + (1+a)*x(2)^4  + (2-a-b )*x(1)^2*x(2)^2 / 2 + 2*b*x(1)*x(2)^3;

c = [2, 3];

obj = c * [a; b];

%moment order
Rx = 1000;
Rab = 100;
K = [f >= 0, norm(x, 2)^2 <= Rx, norm([a;b], 2)^2<= Rab];


m=4;
%opts = sdpsettings('solver', 'moment', 'moment.order', m);
opts = sdpsettings('solver', 'SparsePOP', 'moment.relaxOrder', m);
sol = optimize(K, obj, opts);
[value(a), value(b), value(obj)]

% relax = [3,4,5,6];
% %relax = [3,4];
% 
% %ab_list = cell(length(relax), 1);
% %x_list = cell(length(relax), 1);
% %relax =       
% i=4;
% %i = 1:2;
% %for i = 1:length (relax)
%     m = relax(i);
%     [sol, x, momentdata] = solvemoment(K, obj, [], m);
%     x_list{i} = x;
%     ab_list{i} = [value(a), value(b), value(obj)];
%     time_6 = 931.88;
% %end
