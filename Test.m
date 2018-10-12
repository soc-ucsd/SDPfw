%%
clc;clear
% load TestEx


% ---------------------------------------------------------------------------- %
%                       Conic program with banded SDP
% ---------------------------------------------------------------------------- %
% Parameters
m   = 2;                      % # constraints
K.s = [100];                    % PSD cones
bandWidth = [K.s-1];              % bandWidth for SDP cones


% Setup
fprintf('\nSetting up random conic problem with banded SDP cones, m=%i...',m);
tsetup = tic;
[At,b,c,K] = bandedSDP(m,K,bandWidth);
tsetup = toc(tsetup);


opts.NoP = 13;%K.s;  %% number of partition
opts.socp = 1;

tic
[Anew, bnew, cnew, Knew] = FactorWidth(At,b,c,K,opts);
toc

% opts.socp = 0;
% tic
% [Anew1, bnew1, cnew1, Knew1] = FactorWidth(At,b,c,K,opts);
% toc

%%
% [x,y,info] = sedumi(At,b,c,K);
% %[x,y,info] = sedumi(Anew1,bnew1,cnew1,Knew1);
% [x1,y1,info1] = sedumi(Anew,bnew,cnew,Knew);
% [c'*x,cnew'*x1]
% [info.wallsec,info1.wallsec]

%% Using Mosek
% tic
% [F,h] = sedumi2yalmip(At,b,c,K);   % use yalmip
% opts = sdpsettings('solver','mosek');
% [Model,RECOVERYMODEL,DIAGNOSTIC,INTERNAL] = export(F,h,opts);   % to mosek format
% toc;

prob0          = convert_sedumi2mosek(At', b, c, K);
Tmosek         = tic;
[rcode0, res0] = mosekopt('minimize info', prob0);
Tmosek0        = toc(Tmosek);


prob1          = convert_sedumi2mosek(Anew', bnew, cnew,Knew); 
Tmosek         = tic;
[rcode1, res1] = mosekopt('minimize info', prob1);
Tmosek1        = toc(Tmosek);

[res0.sol.itr.pobjval,res1.sol.itr.pobjval]
[Tmosek0,Tmosek1]
