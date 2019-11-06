%Sea star topology

%One big clique in the middle 'head'
%N different 'arms', each composed of k 'knuckles'
%knuckles communicate with each other across t <=k nodes.
%% Stability Test: Generating data
rng(62, 'twister')

head = 60;      %size of central 'head'
knuckle = 10;    %size of each knuckle
t = 4;          %#links between head and first knuckle
t_k = 2;        %#links between subsequent knuckles
N = 6;          %#arms
k = 6;          %#knuckles per arm

% head = 20;      %size of central 'head'
% knuckle = 6;    %size of each knuckle
% t = 2;          %#links between head and first knuckle
% t_k = 2;        %#links between subsequent knuckles
% N = 6;          %#arms
% k = 3;          %#knuckles per arm

N_state = head + knuckle*N*k;

Gw = sparse(N_state, N_state);

%head is dense
weight_head = 0.8;
weight_knuck = 4;
weight_knuck_k = 4;
weight_knuck_h = 2.5;


Gw(1:head, 1:head) = weight_head;

i_incr = head;
t_incr = 0;
for i = 1:N
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

 figure(1)
subplot(1,2,1)
spy(Gw)
title('Sea Star Interactions', 'fontsize', 18, 'interpreter', 'latex')
subplot(1,2,2)
 plot(graph(Gw, 'omitselfloops'), 'layout', 'force', ...
     'Iterations', 6000, 'UseGravity', 'on', 'WeightEffect', 'inverse')
 axis square
 title('Sea Star Visualization', 'fontsize', 18, 'interpreter', 'latex')



G = (Gw > 0);
Mc    = maximalCliques(G);
%generate system
N = N_state;
n     = randi(10,1,N);
m     = randi(5,1,N);
d     = randi(5,1,N);
Flag  = 3;
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
elseif Flag == 2
    %H2 norm
    %Flag == 2
    Constraint = [Constraint, Sys.globalA*P + P*Sys.globalA' + Sys.globalB*Sys.globalB' + epsilon*eye(sum(n)) <= 0];
    Cost = trace(Sys.globalC*P*Sys.globalC');
    title_str = '$A^T P + P A^T + B B^T \leq -\epsilon I$';
else
    %Hinf norm
    gamma = sdpvar(1);
    
    
    Constraint2 = [[P*Sys.globalA+Sys.globalA'*P, P*Sys.globalB, Sys.globalC'; 
                           Sys.globalB'*P, -gamma*eye(sum(m)), Sys.globalD';
                           Sys.globalC, Sys.globalD, -gamma*eye(sum(d))] + epsilon*eye(sum(n)+sum(m)+sum(d)) <= 0];
    Constraint = [Constraint, Constraint2];
    Cost = gamma;
    %title_str = '$\begin{pmatrix}A^T P + P A^T & P^T B & C^T \\ B^T P& -\gamma I & 0 \\ C & 0 & -\gamma I \end{pmatrix}\leq -\epsilon I$';
    title_str = 'Bounded Real Lemma';
end
% by SeDuMi
opts      = sdpsettings('verbose',1,'solver','sedumi');
[model,~,~,~] = export(Constraint,Cost,opts);

% A = model.A;
% b = model.b;
% c = model.C;
% K.f = model.K.f;K.l = model.K.l;K.q = model.K.q;K.s = model.K.s;


parCoLO.domain    = 1;  % dConvCliqueTree  ---> equalities 
parCoLO.range     = 2;   % rConvMatDecomp   ---> equalities 
parCoLO.EQorLMI   = 1; % CoLOtoEQform     ---> equality standard form
parCoLO.SDPsolver = []; % CoLOtoEQform     ---> equality standard form       
parCoLO.quiet     = 1; % Some peace and quiet       
J.f = length(model.b);

[~,~,~,cliqueDomain,cliqueRange,LOP] = sparseCoLO(model.A',model.b,model.C,model.K,J,parCoLO); 

figure(2)
SP = spones(spones(model.C) + sparse(sum(spones(model.A),2)));  % vector of 1s and 0s
mask1 = reshape(SP(1:model.K.s(1)^2), model.K.s(1), model.K.s(1));
mask2 = reshape(SP(model.K.s(1)^2 + (1:model.K.s(2)^2)), model.K.s(2), model.K.s(2));

subplot(1,2,1)
spy(mask1)
title('$P \geq \epsilon I$', 'interpreter', 'latex', 'Fontsize', 18)
subplot(1,2,2)
spy(mask2)
title(title_str, 'interpreter', 'latex', 'Fontsize', 18)

figure(3)
hold on
plot(sort(LOP.K.s), '.', 'Markersize', 10)
plot(xlim, [22,22], 'k:')
plot(xlim, [100,100], 'k--')
title('Sea Star Clique Sizes', 'fontsize', 18, 'Interpreter', 'latex')
legend('Cliques', 'Size 22', 'Size 100', 'location', 'northwest')
hold off
save('sea_star_Hinf.mat', 'model', 'LOP', 'Sys', 'G', 'n', 'm', 'd')
