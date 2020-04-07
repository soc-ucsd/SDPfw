%Sea star topology

%One big clique in the middle 'head'
%N different 'arms', each composed of k 'knuckles'
%knuckles communicate with each other across t <=k nodes.
%% Stability Test: Generating data
rng(62, 'twister')

star_size = 9;
Flag  = 3; %3 for Hinf

SYNTHESIZE = 1;
VISUALIZE = 1;
BIG_REAL = 0;
if star_size == 9
    %wide_small
    head = 70;      %size of central 'head'
    knuckle = 10;    %size of each knuckle
    t = 4;          %#links between head and first knuckle
    t_k = 4;        %#links between subsequent knuckles
    N_arm = 12;          %#arms
    k = 2;          %#knuckles per arm
    size_str = 'wide_med';
elseif star_size == 8
    %wide_small
    head = 18;      %size of central 'head'
    knuckle = 5;    %size of each knuckle
    t = 2;          %#links between head and first knuckle
    t_k = 2;        %#links between subsequent knuckles
    N_arm = 7;          %#arms
    k = 1;          %#knuckles per arm
    size_str = 'wide_small';
elseif star_size == 7
    %wide
    head = 90;      %size of central 'head'
    knuckle = 9;    %size of each knuckle
    t = 3;          %#links between head and first knuckle
    t_k = 3;        %#links between subsequent knuckles
    N_arm = 12;         %#arms
    k = 2;          %#knuckles per arm
    size_str = 'wide';
elseif star_size == 6
    %Giant
    head = 60;      %size of central 'head'
    knuckle = 10;    %size of each knuckle
    t = 3;          %#links between head and first knuckle
    t_k = 3;        %#links between subsequent knuckles
    N_arm = 6;          %#arms
    k = 6;          %#knuckles per arm
    size_str = 'giant';
elseif star_size == 5
    %Huge
    head = 60;      %size of central 'head'
    knuckle = 10;    %size of each knuckle
    t = 4;          %#links between head and first knuckle
    t_k = 2;        %#links between subsequent knuckles
    N_arm = 6;          %#arms
    k = 6;          %#knuckles per arm
    size_str = 'huge';
elseif star_size == 4
    %Very Large
    head = 40;      %size of central 'head'
    knuckle = 10;    %size of each knuckle
    t = 2;          %#links between head and first knuckle
    t_k = 2;        %#links between subsequent knuckles
    N_arm = 6;          %#arms
    k = 5;          %#knuckles per arm
    size_str = 'verylarge';
elseif star_size == 3
    %Large
    head = 30;      %size of central 'head'
    knuckle = 8;    %size of each knuckle
    t = 2;          %#links between head and first knuckle
    t_k = 2;        %#links between subsequent knuckles
    N_arm = 5;          %#arms
    k = 4;          %#knuckles per arm
    size_str = 'large';
elseif star_size == 2
    %Medium
    head = 20;      %size of central 'head'
    knuckle = 6;    %size of each knuckle
    t = 2;          %#links between head and first knuckle
    t_k = 2;        %#links between subsequent knuckles
    N_arm = 6;          %#arms
    k = 3;          %#knuckles per arm
    size_str = 'medium';
elseif star_size == 1
    %Small
    head = 10;      %size of central 'head'
    knuckle = 4;    %size of each knuckle
    t = 2;          %#links between head and first knuckle
    t_k = 2;        %#links between subsequent knuckles
    N_arm = 3;          %#arms
    k = 2;          %#knuckles per arm
    size_str = 'small';    
elseif star_size == 0
    %Tiny
    head = 4;      %size of central 'head'
    knuckle = 3;    %size of each knuckle
    t = 1;          %#links between head and first knuckle
    t_k = 1;        %#links between subsequent knuckles
    N_arm = 2;          %#arms
    k = 1;          %#knuckles per arm
    size_str = 'tiny';
else
    %Micro
    head = 3;      %size of central 'head'
    knuckle = 1;    %size of each knuckle
    t = 1;          %#links between head and first knuckle
    t_k = 1;        %#links between subsequent knuckles
    N_arm = 2;          %#arms
    k = 1;          %#knuckles per arm
    size_str = 'micro';
end


N_state = head + knuckle*k*N_arm;
Gw = sparse(N_state, N_state);

%head is dense
weight_head = 0.9;
weight_knuck = 4;
weight_knuck_k = 4;
weight_knuck_h = 2;

Gw(1:head, 1:head) = weight_head;

i_incr = head;
t_incr = 0;
for i = 1:N_arm
    %head to knuckle
    k_ind = i_incr + (1:knuckle);
    t_ind = t_incr + (1:t);
    Gw(k_ind, k_ind) = weight_knuck;
    Gw(k_ind, t_ind) = weight_knuck_h;
    Gw(t_ind, k_ind) = weight_knuck_h;
    i_incr = i_incr + knuckle;
    t_incr = t_incr + t;

    %knuckle to knuckle
    for j = 1:(k-1)
        i_prev = i_incr + ((1-t_k):0);
        i_next = i_incr + (1:t_k);

        i_k = i_incr + (1:knuckle);
        Gw(i_k, i_k) = weight_knuck;
        Gw(i_prev, i_next) = weight_knuck_k;
        Gw(i_next, i_prev) = weight_knuck_k;

        i_incr = i_incr + knuckle;
    end
end



if SYNTHESIZE
G = (Gw > 0);
Mc    = maximalCliques(G);
%generate system
N = N_state;
n     = randi(10,1,N);
m     = randi(5,1,N);
d     = randi(5,1,N);

Sys   = GenerateDynamics(G,[],n,m,d,Flag);

%% Call Yalmip to solve the problem
epsilon = 0.01;
P = [];
for i = 1:N
    P = blkdiag(P,sdpvar(n(i)));
end

num_input  = size(Sys.globalB, 2);
num_output = size(Sys.globalC, 1);
num_state  = size(Sys.globalA, 2);

Sys.globalD = sparse(num_output, num_input);

Constraint = [P - epsilon*eye(sum(n)) >= 0];
if Flag == 1
    %stability
    Constraint = [Constraint, Sys.globalA*P + P*Sys.globalA'+ epsilon*eye(sum(n)) <= 0];
    Cost = 0;
    title_str = '$A^T P + P A^T \leq -\epsilon I$';
    flag_str = 'Stab';
elseif Flag == 2
    %H2 norm
    %Flag == 2
    Constraint = [Constraint, Sys.globalA*P + P*Sys.globalA' + Sys.globalB*Sys.globalB' + epsilon*eye(sum(n)) <= 0];
    Cost = trace(Sys.globalC*P*Sys.globalC');
    title_str = '$A^T P + P A^T + B B^T  \epsilon I \leq 0$';
    flag_str = 'H2';
else
    %Hinf norm

    
    if BIG_REAL
        gamma2 = sdpvar(1);
        Constraint2 = [[P*Sys.globalA+Sys.globalA'*P, P*Sys.globalB, Sys.globalC'; 
                                   Sys.globalB'*P, -gamma*eye(sum(m)), Sys.globalD';
                                   Sys.globalC, Sys.globalD, -gamma*eye(sum(d))] + epsilon*eye(sum(n)+sum(m)+sum(d)) <= 0];
        Constraint = [Constraint, Constraint2, gamma2 >= 0];
        Cost = gamma2;
        flag_str = 'Hinf1';
    else
        gamma2 = sdpvar(1);
        Constraint2 = [[P*Sys.globalA+Sys.globalA'*P + Sys.globalC'*Sys.globalC, P*Sys.globalB + Sys.globalC'*Sys.globalD; 
                            Sys.globalB'*P + Sys.globalD'*Sys.globalC, Sys.globalD'*Sys.globalD-gamma2*eye(sum(m))] + epsilon*eye(sum(n)+sum(m)) <= 0];
        Constraint = [Constraint, Constraint2, gamma2 >= 0];
        Cost = gamma2;
        flag_str = 'Hinf0';
    end
    
    %title_str = '$\begin{pmatrix}A^T P + P A^T & P^T B & C^T \\ B^T P& -\gamma I & 0 \\ C & 0 & -\gamma I \end{pmatrix}\leq -\epsilon I$';
    title_str = 'Bounded Real Lemma';
end

dir_str = strcat('sea_star_',flag_str,'_',size_str);
opts      = sdpsettings('verbose',1,'solver','sedumi');
[model_dense,~,~,~] = export(Constraint,Cost,opts);

if ~exist(dir_str, 'dir')
    mkdir(dir_str);
end

%MAKE PLOTS
if VISUALIZE
    close all;

    %figure(1)
    figure('units','normalized','outerposition',[0 0 0.5 0.5])
    subplot(1,2,2)
    spy(Gw)
    title('Sea Star Interactions', 'fontsize', 18, 'interpreter', 'latex')
    subplot(1,2,1)
    plot(graph(Gw, 'omitselfloops'), 'layout', 'force', ...
        'Iterations', 6000, 'UseGravity', 'on', 'WeightEffect', 'inverse')
    axis square
    axis off
    title('Sea Star Visualization', 'fontsize', 18, 'interpreter', 'latex')
    box off
    export_fig(strcat(dir_str, '\\network_plot'), '-pdf', '-png');
    %export_fig strcat(dir_str, '\\network_plot.pdf')
end

% A = model.A;
% b = model.b;
% c = model.c;
% K.f = model.K.f;K.l = model.K.l;K.q = model.K.q;K.s = model.K.s;

%yalmip outputs SDP
parCoLO0.domain    = 1;  % dConvCliqueTree  ---> equalities 
parCoLO0.range     = 2;   % rConvMatDecomp   ---> equalities 
parCoLO0.EQorLMI   = 1; % CoLOtoEQform     ---> LMI standard form
parCoLO0.SDPsolver = []; % CoLOtoEQform     ---> LMI standard form       
parCoLO0.quiet     = 1; % Some peace and quiet 

%yalmip outputs LMI
parCoLO.domain    = 2;  % dConvCliqueTree  ---> equalities 
parCoLO.range     = 1;   % rConvMatDecomp   ---> equalities 
parCoLO.EQorLMI   = 2; % CoLOtoEQform     ---> LMI standard form
parCoLO.SDPsolver = []; % CoLOtoEQform     ---> LMI standard form       
    
parCoLO.quiet     = 1; % Some peace and quiet       
J.f = length(model_dense.b);

[~,~,~,cliqueDomain,cliqueRange,LOP] = sparseCoLO(model_dense.A',model_dense.b,model_dense.C,model_dense.K,J,parCoLO); 
%[~,~,~,cliqueDomain,cliqueRange,LOP_dual] = sparseCoLO(model.A',model.b,model.c,model.K,J,parCoLO0); 

model = struct;
model.A = -LOP.A';
model.b = -LOP.c;
model.c = -LOP.b;
model.K = LOP.J;

fname = strcat(strcat(dir_str, '\\sea_star.mat'));

save(fname, 'model', 'Sys', 'G', 'n', 'm', 'd', 'Gw', 'model_dense')
end

if VISUALIZE
    figure('units','normalized','outerposition',[0 0 0.5 0.5])
    clf
    SP = spones(spones(model_dense.C) + sparse(sum(spones(model_dense.A),2)));  % vector of 1s and 0s
    mask1 = reshape(SP(1:model_dense.K.s(1)^2), model_dense.K.s(1), model_dense.K.s(1));
    mask2 = reshape(SP(model_dense.K.s(1)^2 + (1:model_dense.K.s(2)^2)), model_dense.K.s(2), model_dense.K.s(2));

    subplot(1,2,1)
    spy(mask1)
    title('$P \geq \epsilon I$', 'interpreter', 'latex', 'Fontsize', 18)
    subplot(1,2,2)
    spy(mask2)
    title(title_str, 'interpreter', 'latex', 'Fontsize', 18)
    export_fig(strcat(dir_str, '\\lmi'), '-pdf', '-png');


    figure('units','normalized','outerposition',[0 0 0.5 0.5])
    subplot(1,3,1)
    plot(graph(Gw, 'omitselfloops'), 'layout', 'force', ...
        'Iterations', 10000, 'UseGravity', 'on', 'WeightEffect', 'inverse')
    axis square
    axis off
    title('Sea Star Visualization', 'fontsize', 18, 'interpreter', 'latex')
    box off
    subplot(1,3,2)
    spy(mask1)
    title('$P \geq \epsilon I$', 'interpreter', 'latex', 'Fontsize', 18)
    subplot(1,3,3)
    spy(mask2)
    title(title_str, 'interpreter', 'latex', 'Fontsize', 18)
    export_fig(strcat(dir_str, '\\vis_lmi'), '-pdf', '-png');
    %plot(sort(LOP.J.s), '.', 'Markersize', 10)
%     plot(xlim, [11,11], 'k--')
% %     plot(xlim, [60,60], 'k-.')
% %     plot(xlim, [100,100], 'k-')
%     %plot(xlim, [40,40], 'k-.')
%     %plot(xlim, [75,75], 'k-')
%     plot(xlim, [41,41], 'k-.')
%     plot(xlim, [63,63], 'k-')
    %[N_h,edges] = histcounts(model.K.s, 'BinMethod','integers');
%     subplot(2,1,1)
%         hold on
%     [N_h,edges] = histcounts(model_lower.K.s, 'BinMethod','integers');
%     edges = edges+0.5;
%      yl = [0, max(N_h)];
%      plot([11,11], yl,'k--')
%      plot([60,60],  yl, 'k-.')
%      plot([100,100], yl, 'k-')
% 
%     stem(edges([N_h 0] ~= 0), N_h(N_h ~= 0), '.', 'MarkerSize', 30)
%     title(strcat('Sea Star Lower Bound Clique Sizes $(p=', num2str(length(model_lower.K.s)),')$'), 'fontsize', 18, 'Interpreter', 'latex')
%      legend({'Size 11', 'Size 60', 'Size 100','Cliques'},...
%         'location', 'northeast', 'fontsize', 12)
%     %legend({'Cliques'}, 'location', 'northeast', 'fontsize', 12)
%     %hold off
%     xlabel('Size of Clique')
%     ylabel('Number of Cliques')
%     
%    subplot(2,1,2)
     figure('units','normalized','outerposition',[0 0 0.5 0.4])
    clf
 
    hold on
    [N_h,edges] = histcounts(model.K.s, 'BinMethod','integers');
    edges = edges+0.5;
    yl = [0, max(N_h)];
     plot([11,11], yl,'k--')
     plot([60,60],  yl, 'k-.')
     plot([100,100], yl, 'k-')
   stem(edges([N_h 0] ~= 0), N_h(N_h ~= 0), '.', 'MarkerSize', 30)

     legend({'Size 11', 'Size 60', 'Size 100','Cliques'},...
        'location', 'northeast', 'fontsize', 12)
    xlabel('Size of Clique')
    ylabel('Number of Cliques')
    export_fig(strcat(dir_str, '\\clique_size'), '-pdf', '-png');
    
    %xlabel('Clique  #', 'fontsize', 12)
    %ylabel('Size of Clique', 'fontsize', 12) 
end