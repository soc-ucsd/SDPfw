load('LR_120.mat', 'model_unc', 'model_cons_trans', 'model_csp')

%[r, res_out2] = mosekopt('maximize', res_unc.prob);

C = linspecer(5);
figure(3)
clf

subplot(2,1,1)
hold on
[N_h_unc,edges_unc] = histcounts(model_unc.K.s, 'BinMethod','integers');
 yl_unc = [0, max([N_h_unc, model_unc.K.l])];
 plot([12,12], yl_unc,'k--')
 plot([45,45],  yl_unc, 'k-.')
 plot([100,100], yl_unc, 'k-')
 
 stem([1, edges_unc([N_h_unc 0] ~= 0)+0.5], [model_unc.K.l, N_h_unc(N_h_unc ~= 0)], '.', 'MarkerSize', 40)
 title('Unconstrained clique sizes', 'FontSize', 18, 'Interpreter', 'latex')
 hold off
xlabel('Size of Clique')
ylabel('Number of Cliques')
    
 
 subplot(2,1,2)   
hold on
 [N_h,edges] = histcounts(model_cons_trans.K.s, 'BinMethod','integers');
 stem([1, edges([N_h 0] ~= 0)+0.5], [model_cons_trans.K.l, N_h(N_h ~= 0)], '.', 'MarkerSize', 40)
 yl = [0, max([N_h, model_cons_trans.K.l])];
 plot([12,12], yl,'k--')
 plot([45,45],  yl, 'k-.')
 plot([100,100], yl, 'k-')
 
 hold off
 %title('Lehmer-Rosenbrock Clique Sizes with 120 Variables','FontSize', 18,'Interpreter', 'latex')
 %legend('Unconstrained', 'Box-constrained')
 title('Box-constrained $[2,4]^N$ clique sizes','FontSize', 18,'Interpreter', 'latex')
xlabel('Size of Clique')
ylabel('Number of Cliques')

%ylabel('Number of Cliques')
            