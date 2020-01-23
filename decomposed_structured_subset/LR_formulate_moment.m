lower = sdpvar(1,1);
rng( 40, 'twister');
%N = 60;
%N = 72;
%N = 120;
%N = 36;
N = 8;
%b = floor(N/6);
b = floor(N/3);
CONSTRAINED = 0;

x = sdpvar(N, 1);

%f_R  = sum(100*(x(2:end ) - x(1:end-1) .^2).^2 +  (1 -x(1:end-1 )).^2);

xi = x(1:end-3);
xi1 = x(2:end-2);
xi2 = x(3:end-1);
xi3 = x(4:end);
%there was a typo here
f_R = sum(10*((xi2 + 2*xi1 - xi.^2).^2) + (1-xi -xi3).^2);

%mBasis = monolist(x(1:b),2);
%Q = rand(length(mBasis));

%Q = hadamard(
%c = rand(length(mBasis) ,1);
%Q = Q*Q';
%f_Q = (Q'*mBasis).^2  + c'*mBasis;
%f_Q = (Q'*mBasis).^2 ;
mBasis = x(1:b);
Q_P = gallery('lehmer', length(mBasis));
f_Q = mBasis' *Q_P * mBasis;
%f_Q = 0;
f = f_R + f_Q;
 

if CONSTRAINED  
    blo = -0.5;
    bhi = 0.5;
    for i = 1:N
        F = [F; xi >= blo; xi <= bhi];
    end
else
    F = [];
end
obj =  f;

%Moment relaxation
opts_mom = sdpsettings('solver', 'moment');

d=2;
[sol,x2,momentdata] = solvemoment(F, obj, [], d);
obj2 = value(obj);
[sol3,x3,momentdata3] = solvemoment(F, obj, [], d+1);
obj3 = value(obj);

% %CSP
% opts_csp = sdpsettings('solver','SEDUMI');
% %[model, recoverymodel] = export(F , obj, opts);
% 
% %opts_csp = opts;
% opts_csp.sos.csp = 1;
% 
% 
% [F_m, obj_m, monomials] = sosmodel(F, obj, opts_csp, [lower; c]);
% [model_csp, recoverymodel_csp] = export(F_m, obj_m, opts_csp);
% %[x_csp, y_csp, info_csp] = sedumi(model_csp.A', model_csp.b, model_csp.C, model_csp.K);
% %cost_csp = model_csp.C'*x_csp;
% 
% 
% [model_curr.A, model_curr.b, model_curr.C, model_curr.K, ~] = decomposed_subset(model_csp.A,model_csp.b,model_csp.C,model_csp.K, 20);
% prob_curr = convert_sedumi2mosek(model_curr.A, model_curr.b, model_curr.C, model_curr.K);
% 
% [r,res] = mosekopt('minimize',prob_curr);




% opts2 = sdpsettings('solver', 'Mosek');
% sol = optimize(F_m, obj_m, opts2);
% 
% Ks = model_csp.K.s';
% mKs = max(model_csp.K.s);
% 
% figure(1)
% clf
% hold on
% [N_h,edges] = histcounts(Ks, 'BinMethod','integers');
% yl = [0, max(N_h)];
% plot([30,30], yl,'k--')
% plot([100,100],  yl, 'k-.')
% stem(edges([N_h 0] ~= 0), N_h(N_h ~= 0), '.', 'MarkerSize', 30)
% hold off
% 
% title('LR Clique Sizes', 'fontsize', 18, 'Interpreter', 'latex')
% 
% 
% legend({'Size 20', 'Size 70', 'Cliques'},...
% 'location', 'northeast', 'fontsize', 12)
% %legend({'Cliques'}, 'location', 'northeast', 'fontsize', 12)
% %hold off
% xlabel('Size of Clique')
% ylabel('Number of Cliques')
% 
%save('polynomial_big_clique.mat', 'N', 'b', 'model_csp')

[N_h_unc,edges_unc] = histcounts(model_csp.K.s, 'BinMethod','integers');
 yl_unc = [0, max(N_h_unc)];
 plot([11,11], yl_unc,'k--')
 plot([45,45],  yl_unc, 'k-.')
 plot([100,100], yl_unc, 'k-')
 
 stem([1, edges_unc([N_h_unc 0] ~= 0)+0.5], [model_csp.K.l, N_h_unc(N_h_unc ~= 0)], '.', 'MarkerSize', 40)
 title('Unconstrained clique sizes', 'FontSize', 18, 'Interpreter', 'latex')
 hold off
xlabel('Size of Clique')
ylabel('Number of Cliques')