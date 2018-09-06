%%
clc;clear
% load TestEx


% ---------------------------------------------------------------------------- %
%                       Conic program with banded SDP
% ---------------------------------------------------------------------------- %
% Parameters
m   = 3;                      % # constraints
K.s = [4];                    % PSD cones
bandWidth = [K.s-1];              % bandWidth for SDP cones


% Setup
fprintf('\nSetting up random conic problem with banded SDP cones, m=%i...',m);
tsetup = tic;
[At,b,c,K] = bandedSDP(m,K,bandWidth);
tsetup = toc(tsetup);


% alpha = 10;
% alpha = [2 3 3 2];
% alpha = [5 3 1 1];
opts.NoP = K.s ;

[Anew, bnew, cnew, Knew] = FactorWidth(At,b,c,K,opts);

[x,y,info] = sedumi(At,b,c,K);

[x1,y1,info1] = sedumi(Anew,bnew,cnew,Knew);

[c'*x,cnew'*x1]

[info.wallsec,info1.wallsec]