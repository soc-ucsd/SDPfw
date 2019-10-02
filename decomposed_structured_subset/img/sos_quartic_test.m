sos_quartic_prep_3;

load('quartic_process_3.mat');
%load('quartic_process0.mat');
%original
%f = x1^4 + x2^4 + a*x1^3*x2 + (2-a-b )*x1^2*x2^2 / 2 + 2*b*x1*x2^3;

%sparse, missing x1^3*x2 term
%f = x1^4 + (1+a)*x2^4  + (2-a-b )*x1^2*x2^2 / 2 + 2*b*x1*x2^3;
pars.fid = 0;
%c0 = c;

%theta = pi/6;
N =  10;
%N = 4;
%N = 300;  
th = linspace(0, 2*pi, N);

%output dumps 
PSD.a = zeros(N, 1);
PSD.b = zeros(N, 1);
PSD.Q = cell(N, 1);

AGLER.a  = zeros(N, 1);
AGLER.b  = zeros(N, 1);
AGLER.Q  = cell(N, 1);

%iterate through solvers

for i = 1:N
    theta = th(i);
    model.c(1) = cos(theta);
    model.c(2) = sin(theta);

    model2.c(1) = cos(theta);
    model2.c(2) = sin(theta);

    
    %PSD
    [x,y,info] = sedumi(model.At,model.b,model.c,model.K,model.pars);
    PSD.a(i) = x(1);
    PSD.b(i) = x(2);
    K = model.K;
    %PSD.Q{i} = reshape(x(K.f+1:end), K.s, K.s );
    PSD.Q{i} = x(K.f+1:end);
    %SDD(E, ?)
    
    [x,y,info] = sedumi(model2.At,model2.b,model2.c,model2.K,model2.pars);
    K = model2.K;
    AGLER.a(i) = x(1);
    AGLER.b(i) = x(2);
    AGLER.Q{i} = x(K.f+1:end);
    %AGLER.Q1{i} = reshape(x((K.f+1):(K.f+K.s(1)^2)), K.s(1), K.s(1));
    %AGLER.Q2{i} = reshape(x((K.f+K.s(1)^2 +1):end), K.s(2), K.s(2));
    
end

figure(1)
clf

%convex hulls
PSD.conv  = convhull(PSD.a, PSD.b);
AGLER.conv = convhull(AGLER.a, AGLER.b);

hold on 
plot(PSD.a(PSD.conv), PSD.b(PSD.conv))
plot(AGLER.a(AGLER.conv), AGLER.b(AGLER.conv))
hold off
