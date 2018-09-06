%%
clc;clear
% load TestEx


% ---------------------------------------------------------------------------- %
%                       Conic program with banded SDP
% ---------------------------------------------------------------------------- %
% Parameters
m   = 3;                      % # constraints
K.s = [1000];                    % PSD cones
bandWidth = [K.s-1];              % bandWidth for SDP cones


% Setup
fprintf('\nSetting up random conic problem with banded SDP cones, m=%i...',m);
tsetup = tic;
[At,b,c,K] = bandedSDP(m,K,bandWidth);
tsetup = toc(tsetup);


% alpha = 10;
% alpha = [2 3 3 2];
% alpha = [5 3 1 1];
opts.NoP = 20 ;  %% number of partition

[Anew, bnew, cnew, Knew] = FactorWidth(At,b,c,K,opts);


%%
% [x,y,info] = sedumi(At,b,c,K);
% [x1,y1,info1] = sedumi(Anew,bnew,cnew,Knew);
% [c'*x,cnew'*x1]
% [info.wallsec,info1.wallsec]

%% Using Mosek
prob0          = convert_sedumi2mosek(At', b, c, K);
Tmosek         = tic;
[rcode0, res0] = mosekopt('minimize info', prob0);
Tmosek0        = toc(Tmosek);

prob1          = convert_sedumi2mosek(Anew', bnew, cnew, Knew);
Tmosek         = tic;
[rcode1, res1] = mosekopt('minimize info', prob1);
Tmosek1        = toc(Tmosek);

[res0.sol.itr.pobjval,res1.sol.itr.pobjval]
[Tmosek0,Tmosek1]
