M = sdpvar(4, 4);

a = sdpvar(1);
b = sdpvar(1);
%M(2,2) == 2-b,
F = [
    M(1,1) == 1,    
    M(2,2) == 2,
    M(3,3) == 5,
    M(4,4) == 2,        
    M(1,2) == 0.5+a,
    M(2,3) == -2*a,
    M(2,4) == a + b,    
    M(3,4) == b/2,   
];

Mr = reshape(M, [], 1);
rays = dd_extreme_rays(4);


F_dd_star = [F, rays'*Mr >= 0];

Mr3 = reshape(M(2:4, 2:4), [], 1);
Mr2 = reshape(M(1:2, 1:2), [], 1);
rays3 = dd_extreme_rays(3);
rays2 = dd_extreme_rays(2);
F_cdd_star = [F, rays3'*Mr3 >= 0, rays2'*Mr2 >= 0];


F_psd = [F, M >= 0];


c3 = sdpvar( 9, 1 );
c2 = sdpvar(4, 1);
F_cdd = [F, Mr3 == rays3*c3, Mr2 == rays2*c2, c3 >= 0, c2 >= 0];

c = sdpvar(16, 1);
F_dd = [F, Mr == rays*c, c >= 0];


% F_dd = F;
% for i = 1:4
%     k = (1:4);
%     k(i) = [];
%     F_dd = [F_dd, M(i, i) >= sum(abs(M(i, k)))];
% end



%F = [F, M >= 0];

th = sdpvar(1);
obj = cos(th) * a + sin(th)*b;

obj0 = cos(3*pi/4)*a + sin(3*pi/4)*b;

P_dd_star = optimizer(F_dd_star, obj, sdpsettings('solver', 'mosek'), ...
    th, {a, b});

P_cdd_star = optimizer(F_cdd_star, obj, sdpsettings('solver', 'mosek'), ...
    th, {a, b});


P_psd = optimizer(F_psd, obj, sdpsettings('solver', 'mosek'), ...
    th, {a, b});

P_cdd = optimizer(F_cdd, obj, sdpsettings('solver', 'mosek'), ...
    th, {a, b});

P_dd = optimizer(F_dd, obj, sdpsettings('solver', 'mosek'), ...
    th, {a, b});

% [Fddd, objddd] = dualize(F_dd_star, obj0);
% model_dd_star_dual = export(Fddd, -objddd, sdpsettings('solver', 'sedumi'));
%model_dd_star_dual = export(F_dd_star, obj0, sdpsettings('solver', 'sedumi'));

%N = 40;
N = 200;

th_list = linspace(0, 2*pi, N);
out_dd_star = struct;
out_dd_star.a = zeros(N, 1);
out_dd_star.b = zeros(N, 1);

out_cdd_star = struct;
out_cdd_star.a = zeros(N, 1);
out_cdd_star.b = zeros(N, 1);


out_psd = struct;
out_psd.a = zeros(N, 1);
out_psd.b = zeros(N, 1);


out_cdd = struct;
out_cdd.a = zeros(N, 1);
out_cdd.b = zeros(N, 1);


out_dd = struct;
out_dd.a = zeros(N, 1);
out_dd.b = zeros(N, 1);

for i = 1:N
    th_curr = th_list(i);
    xval = P_dd_star(th_curr);
    out_dd_star.a(i) = xval{1};
    out_dd_star.b(i) = xval{2};
    
    xval = P_cdd_star(th_curr);
    out_cdd_star.a(i) = xval{1};
    out_cdd_star.b(i) = xval{2};
    
    
    xval = P_psd(th_curr);
    out_psd.a(i) = xval{1};
    out_psd.b(i) = xval{2};
    
    xval = P_cdd(th_curr);
    out_cdd.a(i) = xval{1};
    out_cdd.b(i) = xval{2};
    
    xval = P_dd(th_curr);
    out_dd.a(i) = xval{1};
    out_dd.b(i) = xval{2};
end

C =  linspecer(8);

figure(55)
clf
plot_region(out_dd_star, C( 8, :))
plot_region(out_psd, [0, 0,  0])
plot_region(out_dd,C(1, :))
axis square
axis off


figure(56)
clf
plot_region(out_dd_star, C( 6, :))
plot_region(out_cdd_star, C( 5, :))
plot_region(out_psd, [0, 0,  0])
plot_region(out_cdd, C( 2, :))
plot_region(out_dd,C(1, :))
axis square
axis off


% th = 3*pi/4;
% 
% obj = cos(th) * a + sin(th)*b;
% 
% model_dd_star_dual2 = export(F, obj, sdpsettings('solver', 'sedumi', 'dualize', 2));
% 
% 
% 
% 
% sol = optimize(F, obj);
% 
