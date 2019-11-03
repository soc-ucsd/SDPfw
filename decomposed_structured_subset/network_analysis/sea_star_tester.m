%Sea star topology

%One big clique in the middle 'head'
%N different 'tentacles', each composed of k 'knuckles'
%knuckles communicate with each other across t <=k nodes.


head = 30;
knuckle = 10;
t = 3;
N = 3;

N_state = head + knuckle*N;

G = sparse(N_state, N_state);

%head is dense
G(1:head, 1:head) = 1;

i_incr = head;
t_incr = 0;
for i = 1:N
    k_ind = i_incr + (1:knuckle);
    t_ind = t_incr + (1:t);
    G(k_ind, k_ind) = 1;
    G(k_ind, t_ind) = 1;
    G(t_ind, k_ind) = 1;
    i_incr = i_incr + knuckle;
    t_incr = t_incr + t;
end

%spy(G)

%generate system
N = N_state;
n     = randi(10,1,N);
m     = randi(5,1,N);
d     = randi(5,1,N);
Flag  = 2;
Sys   = GenerateDynamics(Gp,[],n,m,d,Flag);