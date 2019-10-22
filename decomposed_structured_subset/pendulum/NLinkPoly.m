%if (nargin<1) N=max(ceil(8*rand),2); end


N = 2;
p = PlanarNLink(N);
%v = p.constructVisualizer();
%v.axis=.5*[-1 1 -1 1]+N*1.5*[-1 1 -1 1];

Q = blkdiag(10*eye(N),eye(N)); R = eye(N-1);
x0 = Point(p.getStateFrame,[pi;zeros(2*N-1,1)]);
[c,V0] = tilqr(p,x0,Point(p.getInputFrame,zeros(N-1,1)),Q,R);

sys=feedback(p,c);
p = taylorApprox(p,0,x0,[],3);
f = clean(p.getPolyDynamics());

x = V0.getFrame.getPoly;
u = sys.getInputFrame.getPoly;

feed = c.D*x;

f_cl = f;
for i = 1:length(u)
    f_cl = subs(f_cl, u(i), feed(i));
end


%copying from maxROAFeedback.m
S0 = V0.S;
K = c.D;
%save(strcat('N_link_system_',num2str(N),'.mat'), 'f', 'f_cl', 'S0', 'x', 'u', 'K')


% prog = mssprog;
% 
% %[prog,Phi] = new(prog,length(x),'psd');
% %V = x'*Phi*x;
% %Vdot = diff(V,x)*f;
% %V = clean(V);
% %Vdot = clean(Vdot);
% V = x'*V0.S*x/2;
% Vdot = x'*V0.S*f;


% [prog,rho] = new(prog,1,'pos');
% 
% % Declare multipliers
%  L1m = monomials(x,2:2);
%  [prog,l1] = new(prog,length(L1m),'free');
%  L = l1'*L1m;
%  prog.sos = L;
% %L = 1;
% 
% % Declare SOS conditions
% %prog.eq = trace(Phi) - trace(S0);
% %prog.sos = -Vdot + L1*(V-rho);
% prog.sos = (x'*x)*(V - rho) + L*Vdot;
% 
% pars.fid = 1;
% %[prog,info] = sedumi(prog,-rho,1,pars,1);
% [model.A,model.b,model.c,model.K,h]=mssprog2sedumi(prog,-rho);
% %[xs, ys, info] = sedumi(model.A', model.b, model.c, model.K, pars);
% 
% % 
% parCoLO.domain    = 1;  % dConvCliqueTree  ---> equalities 
% parCoLO.range     = 2;   % rConvMatDecomp   ---> equalities 
% parCoLO.EQorLMI   = 1; % CoLOtoEQform     ---> equality standard form
% parCoLO.SDPsolver = []; % CoLOtoEQform     ---> equality standard form       
% parCoLO.quiet     = 1; % Some peace and quiet       
% J.f = length(model.b);
% [~,~,~,cliqueDomain,cliqueRange,LOP] = sparseCoLO(model.A,model.b,model.c,model.K,J,parCoLO); 
% lks = LOP.K.s;
% [xs, ys, info] = sedumi(LOP.A', LOP.b, LOP.c, LOP.K, pars);

%max rho
%st
% (x'x)(V(x)-rho) + L(x) Vdot(x) >= 0
%V(x) = x' S x >= 0
%
%V, L of degree 2