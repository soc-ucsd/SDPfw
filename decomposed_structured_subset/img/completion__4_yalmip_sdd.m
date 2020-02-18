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


F_sdd_star = F;
F_csdd_star = F;

for i = 1:4
    for j = i+1:4
        ind = [i, j];
        F_sdd_star = [F_sdd_star, M(ind, ind) >= 0];
        
        if i ~= 1 || (i==1 && j == 2)
            F_csdd_star = [F_csdd_star, M(ind, ind) >= 0];
        end
    end
end
%Mi = sdpvar(4, 4, 6);





F_psd = [F, M >= 0];

th = sdpvar(1);
obj = cos(th) * a + sin(th)*b;

obj0 = cos(3*pi/4)*a + sin(3*pi/4)*b;

P_sdd_star = optimizer(F_sdd_star, obj, sdpsettings('solver', 'mosek'), ...
    th, {a, b});

P_csdd_star = optimizer(F_csdd_star, obj, sdpsettings('solver', 'mosek'), ...
    th, {a, b});

P_psd = optimizer(F_psd, obj, sdpsettings('solver', 'mosek'), ...
    th, {a, b});

% [Fddd, objddd] = dualize(F_dd_star, obj0);
% model_dd_star_dual = export(Fddd, -objddd, sdpsettings('solver', 'sedumi'));
%model_dd_star_dual = export(F_dd_star, obj0, sdpsettings('solver', 'sedumi'));

N = 40;
%N = 200;

th_list = linspace(0, 2*pi, N);



out_sdd_star = struct;
out_sdd_star.a = zeros(N, 1);
out_sdd_star.b = zeros(N, 1);

out_csdd_star = struct;
out_csdd_star.a = zeros(N, 1);
out_csdd_star.b = zeros(N, 1);


out_psd = struct;
out_psd.a = zeros(N, 1);
out_psd.b = zeros(N, 1);


for i = 1:N
    th_curr = th_list(i);
    xval = P_sdd_star(th_curr);
    out_sdd_star.a(i) = xval{1};
    out_sdd_star.b(i) = xval{2};
    
    xval = P_csdd_star(th_curr);
    out_csdd_star.a(i) = xval{1};
    out_csdd_star.b(i) = xval{2};
    
    
    xval = P_psd(th_curr);
    out_psd.a(i) = xval{1};
    out_psd.b(i) = xval{2};   
end

C =  linspecer(8);

figure(55)
clf
plot_region(out_sdd_star, C( 7, :))
plot_region(out_csdd_star, C( 8, :))
plot_region(out_psd, [0, 0,  0])

axis square
axis off
% 
% 
% figure(56)
% clf
% plot_region(out_dd_star, C( 6, :))
% plot_region(out_cdd_star, C( 5, :))
% plot_region(out_psd, [0, 0,  0])
% plot_region(out_cdd, C( 2, :))
% plot_region(out_dd,C(1, :))
% axis square
% axis off
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
