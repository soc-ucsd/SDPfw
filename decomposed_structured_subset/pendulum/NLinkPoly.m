%if (nargin<1) N=max(ceil(8*rand),2); end


N = 3;
p = PlanarNLink(N);
%v = p.constructVisualizer();
%v.axis=.5*[-1 1 -1 1]+N*1.5*[-1 1 -1 1];

Q = blkdiag(10*eye(N),eye(N)); R = eye(N-1);
x0 = Point(p.getStateFrame,[pi;zeros(2*N-1,1)]);
[c,V] = tilqr(p,x0,Point(p.getInputFrame,zeros(N-1,1)),Q,R);

sys=feedback(p,c);
p = taylorApprox(p,0,x0,[],3);
f = clean(p.getPolyDynamics());

x = V.getFrame.getPoly;
u = sys.getInputFrame.getPoly;

feed = -c.D*x;

f_cl = f;
for i = 1:length(u)
    f_cl = subs(f_cl, u(i), feed(i));
end


%copying from maxROAFeedback.m
S0 = V.S;
save(strcat('N_link_system_',num2str(N),'.mat'), 'f_cl', 'V0')


prog = mssprog;

[prog,Phi] = new(prog,length(x),'psd');
V = x'*Phi*x;
Vdot = diff(V,x)*f;
V = clean(V);
Vdot = clean(Vdot);
[prog,rho] = new(prog,1,'pos');

% Declare multipliers
L1m = monomials(x,2:2);
[prog,l1] = new(prog,length(L1m),'free');
L = l1'*L1m;
prog.sos = L;

% Declare SOS conditions
prog.eq = trace(Phi) - trace(S0);
%prog.sos = -Vdot + L1*(V-rho);
prog.sos = (x'*x)*(V - rho) + L*Vdot;

pars.fid = 1;
[prog,info] = sedumi(prog,-rho,1,pars,1);
%max rho
%st
% (x'x)(V(x)-rho) + L(x) Vdot(x) >= 0
%V(x) = x' S x >= 0
%
%V, L of degree 2