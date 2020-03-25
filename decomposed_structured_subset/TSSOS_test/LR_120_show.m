load('LR_120.mat', 'model_unc')


load('LR120_box_1_2.mat', 'model_cons_trans')
%[r, res_out2] = mosekopt('maximize', res_unc.prob);

load('LR_120_sos_csp.mat', 'model_csp')

figure(2)
clf
[N_h_csp,edges_csp] = histcounts(model_csp.K.s, 'BinMethod','integers');
stem([1, edges_csp([N_h_csp 0] ~= 0)+0.5], [model_csp.K.l, N_h_csp(N_h_csp ~= 0)], '.', 'MarkerSize', 40)
box off
xticks([1, 15, 100, max(edges_csp)-0.5])
yticks([0, max([model_csp.K.l, N_h_csp])])
set(gca,'XTickLabel',get(gca,'XTickLabel'), 'fontsize',14)
set(gca,'YTickLabel',get(gca,'YTickLabel'), 'fontsize',14)
xlabel('Size')
ylabel('Frequency')
title('Cliques for Unconstrained CSP Problem','FontSize', 18,'Interpreter', 'latex')


figure(2)
clf
[N_h_csp,edges_csp] = histcounts(model_lower.K.s, 'BinMethod','integers');
stem([1, edges_csp([N_h_csp 0] ~= 0)+0.5], [model_lower.K.l, N_h_csp(N_h_csp ~= 0)], '.', 'MarkerSize', 40)
box off
xticks([11, 60, 100, max(edges_csp)-0.5])
yticks([0, max([model_lower.K.l, N_h_csp])])
set(gca,'XTickLabel',get(gca,'XTickLabel'), 'fontsize',14)
set(gca,'YTickLabel',get(gca,'YTickLabel'), 'fontsize',14)
xlabel('Size')
ylabel('Frequency')
title('Sea Star Clique Sizes ($p=1755$)','FontSize', 18,'Interpreter', 'latex')





%C = linspecer(5);
figure(3)
clf

%ha = tight_subplot(1,2,[.01 .03],[.15 .15],[.1 .1]);
%ha = tight_subplot(2,1,[.03 .03],[.15 .15],[.1 .1]);
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

%ha = 

%axes(ha(1))
subplot(1,2,1)
%subplot(2,1,1)
[N_h_unc,edges_unc] = histcounts(model_unc.K.s, 'BinMethod','integers');
stem([1, edges_unc([N_h_unc 0] ~= 0)+0.5], [model_unc.K.l, N_h_unc(N_h_unc ~= 0)], '.', 'MarkerSize', 40)
box off
xticks([1, 6, 12, 45, 100, max(edges_unc)-0.5])
yticks([0, max([model_unc.K.l, N_h_unc])])
set(gca,'XTickLabel',get(gca,'XTickLabel'), 'fontsize',14)
set(gca,'YTickLabel',get(gca,'YTickLabel'), 'fontsize',14)
xlabel('Size')
ylabel('Frequency')
title('Cliques for Unconstrained Problem','FontSize', 18,'Interpreter', 'latex')

%axes(ha(2))
subplot(1,2,2)
%subplot(2,1,2)
[N_h,edges] = histcounts(model_cons_trans.K.s, 'BinMethod','integers');
stem([1, edges([N_h 0] ~= 0)+0.5], [model_cons_trans.K.l, N_h(N_h ~= 0)], '.', 'MarkerSize', 40)
box off
xticks([1, 6, 12, 45, 100, max(edges)-0.5])
yticks([0, max([model_cons_trans.K.l, N_h])])
set(gca,'XTickLabel',get(gca,'XTickLabel'), 'fontsize',16)
set(gca,'YTickLabel',get(gca,'YTickLabel'), 'fontsize',16)
xlabel('Size')
%ylabel('Frequency')
title('Cliques for Constrained Problem','FontSize', 18,'Interpreter', 'latex')

% 
% subplot(2,1,1)
% hold on
% [N_h_unc,edges_unc] = histcounts(model_unc.K.s, 'BinMethod','integers');
%  yl_unc = [0, max([N_h_unc, model_unc.K.l])];
%  plot([12,12], yl_unc,'k--')
%  plot([45,45],  yl_unc, 'k-.')
%  plot([100,100], yl_unc, 'k-')
%  
%  stem([1, edges_unc([N_h_unc 0] ~= 0)+0.5], [model_unc.K.l, N_h_unc(N_h_unc ~= 0)], '.', 'MarkerSize', 40)
%  title('Unconstrained clique sizes', 'FontSize', 18, 'Interpreter', 'latex')
%  hold off
% xlabel('Size of Clique')
% ylabel('Number of Cliques')
%     
%  
%  subplot(2,1,2)   
% hold on
%  [N_h,edges] = histcounts(model_cons_trans.K.s, 'BinMethod','integers');
%  stem([1, edges([N_h 0] ~= 0)+0.5], [model_cons_trans.K.l, N_h(N_h ~= 0)], '.', 'MarkerSize', 40)
%  yl = [0, max([N_h, model_cons_trans.K.l])];
%  plot([12,12], yl,'k--')
%  plot([45,45],  yl, 'k-.')
%  plot([100,100], yl, 'k-')
%  
%  hold off
%  %title('Lehmer-Rosenbrock Clique Sizes with 120 Variables','FontSize', 18,'Interpreter', 'latex')
%  %legend('Unconstrained', 'Box-constrained')
%  title('Box-constrained $[1,2]^{120}$ clique sizes','FontSize', 18,'Interpreter', 'latex')
% xlabel('Size of Clique')
% ylabel('Number of Cliques')

%ylabel('Number of Cliques')
            