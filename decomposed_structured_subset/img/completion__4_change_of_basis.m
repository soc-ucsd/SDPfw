OPT  = 0;
DRAW = 1;


load('sos_quartic_200.mat', 'N', 'model_split', 'out')

if OPT


th = linspace(0,2*pi, N);

ind_star = 80;
theta_star = th(ind_star);

model_split.c =model_split.C;
model_split.At = model_split.A';
model_split = rmfield(model_split, 'C');
model_split = rmfield(model_split, 'A');


%start with initial feasible point
%[x_feas, y_feas, info0] = sedumi(model_split.At, model_split.b, model_split.c, model_split.K, model_split.pars);

%cone = {'dd', 'dd'};

%[model_feas, x_fake] = basis_change(x_feas, model_split);

reg_0 = draw_feasibility(model_split, 'dd', th);

x1 = reg_0.x{ind_star};
[model_split_1, x_fake] = basis_change(x1, model_split);
reg_1 = draw_feasibility(model_split_1, 'dd', th);

x2 = reg_1.x{ind_star};
[model_split_2, x_fake] = basis_change(x2, model_split_1);
reg_2 = draw_feasibility(model_split_2, 'dd', th);

x3 = reg_2.x{ind_star};
[model_split_3, x_fake] = basis_change(x3, model_split_2);
reg_3 = draw_feasibility(model_split_3, 'dd', th);
end

if DRAW
    C = linspecer(10);
    r = 0.8;
figure(1)
clf
hold on
plot_region(out.psd, [0,0,0])
plot_region(reg_0, C(2, :))
plot_region(reg_1, C(7, :))
plot_region(reg_2, C(8, :))
plot_region(reg_3, C(9, :))
%h = quiver(0.5, 2.5, cos(theta_star)*0.8, sin(theta_star)*0.8, 'k', 'MaxHeadSize', 2, 'LineWidth', 2);
%h.Head.LineStyle = 'solid'; 
% headWidth = 8;
% headLength = 8;
%         ah = annotation('arrow',...
%             'headStyle','cback3','HeadLength',headLength,'HeadWidth',headWidth);
%         set(ah,'parent',gca);
%         set(ah,'position',[0, 2.5, cos(theta_star)*r, sin(theta_star)*r]);

A0 = [0, 2.5];
A1 = A0 + [cos(theta_star)*r, sin(theta_star)*r];
arrow(A0, A1, 'width', 2, 'Length', 40);
text(0., 2.5, '$\langle C, X \rangle$', 'Interpreter', 'latex', 'Fontsize', 50)
% 
% hold off
%xlabel('$a$', 'interpreter', 'latex', 'Fontsize', 16)
%ylabel('$b$', 'interpreter', 'latex', 'Fontsize', 16)
axis square
axis off

end