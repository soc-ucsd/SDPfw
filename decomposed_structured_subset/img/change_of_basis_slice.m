%Plots of (decomposed) structured subsets

N = 150;
%rng(2019,'twister')
rng(200,'twister')
%rng(402, 'twister')

%Diagonal

D = [1 0 0;
     0 1 0;
     0 0 1;
     0 0 0;
     0 0 0];

%Diagonally Dominant 
 
DD_aug = [1  1 0  0;
      1  1 1  1;
      0  0 1  1;
      1 -1 0  0;
      0  0 1 -1]/2;

DD = [D DD_aug];  
  
  
CDD_aug = [1   1  1  1;
           1   1  1  1;
           1   1  1  1;
           1   1 -1 -1;
           1  -1  1 -1]/2;
 
CDD = [DD CDD_aug];   

%PSD
R = 1;
theta = linspace(0, 2*pi, N);
phi = linspace(-pi, pi, N);

Xr = cos(phi)'*cos(theta);
X = reshape(Xr, [], 1);
Yr = cos(phi)'*sin(theta);
Y = reshape(Yr, [], 1);
Z = repmat(sin(phi)', N, 1);

PSD = [X.*X, Y.*Y,  Z.*Z,  X.*Y,  Y.*Z]';

%Find affine subspace through I
I = [1; 1; 1; 0; 0];

A = randn(5, 2);

A0 = abs(randn(5, 1));
b0 = A0'*I;

%constraints and cost
K = null(A');
%c = rand(5, 1);
%c =  ones(5, 1);
%c = [1; 0; 0; 0; 0] ;
%c = [-1; 0; 1; -1; -2] ;
c = A*[ -1.5; 2.5];
%p: weights on extreme rays
c_dd = c'*DD;
At_dd = A0'*DD;
%b_dd = [b0; zeros(size(K, 2),1)];
c_cdd = c'*CDD;
At_cdd = A0'*CDD;

%PSD optimization
K_psd.s = 3;
%prepare for equality constraints 
A0_psd = A0;
A0_psd(4) = A0_psd(4)/2;
A0_psd(5) = A0_psd(5)/2;
A0_psd = reshape(v2m(A0_psd), [], 1);
c_psd = reshape(v2m(c), [], 1);

[x_psd, y_psd, info_psd] = sedumi(A0_psd, b0, c_psd, K_psd);
X_psd = reshape(x_psd, 3, 3);
xv_psd = m2v(X_psd);
cost_psd = c_psd'*x_psd;

%DD optimization
n_dd = size(DD, 2);
K_dd.l = n_dd;

[p_dd, y_dd, ~] = sedumi(At_dd, b0, c_dd, K_dd);
 
x_dd = DD*p_dd;
X_dd = v2m(x_dd)
%L_dd = chol(X_dd);
cost_dd = c'*x_dd;

%CDD optimization
n_cdd = size(CDD, 2);
K_cdd.l = n_cdd;
[p_cdd, y_cdd, ~] = sedumi(At_cdd, b0, c_cdd, K_cdd);
 
x_cdd = CDD*p_cdd;
X_cdd = v2m(x_cdd)
%L_dd = chol(X_dd);
cost_cdd = c'*x_cdd;


%intersection of points in cone with affine subspace
D_proj0    = D.*b0./(A0'*D);
DD_proj0   = DD.*b0./(A0'*DD);
CDD_proj0  = CDD.*b0./(A0'*CDD);
PSD_proj0  = PSD.*b0./(A0'*PSD);


%project onto subspace
DD_proj  = A\DD_proj0;
CDD_proj = A\CDD_proj0;
PSD_proj = A\PSD_proj0;
I_proj = A\I;


opt_DD  = A\x_dd;
opt_CDD = A\x_cdd;
opt_PSD = A\xv_psd;

% D_proj = A\D;
% DD_proj  = A\DD;
% CDD_proj = A\CDD;
% PSD_proj = A\PSD;
% SDD_proj = A\SDD;
% CSDD_proj = A\CSDD;

%convex hulls
kDD =  convhull(DD_proj(1, :),    DD_proj(2, :));
kCDD  = convhull(CDD_proj(1, :),  CDD_proj(2, :));
kPSD  = convhull(PSD_proj(1, :),  PSD_proj(2, :));

%plot
C = linspecer(5);
figure(1)
clf
hold on

%extreme rays
plot(DD_proj(1, kDD), DD_proj(2, kDD),'color', C(1, :), 'linewidth', 2)
plot(CDD_proj(1, kCDD), CDD_proj(2, kCDD),'color', C(2, :),'linewidth', 2)
plot(PSD_proj(1, kPSD), PSD_proj(2, kPSD),'k', 'linewidth', 2)
text(I_proj(1), I_proj(2), 'I', 'interpreter', 'latex', 'Fontsize', 20)

%optimal points
plot(opt_DD(1), opt_DD(2), '*', 'color', C(1, :), 'MarkerSize', 20, 'linewidth', 2)
plot(opt_CDD(1), opt_CDD(2), '*', 'color', C(2, :), 'MarkerSize', 20,'linewidth', 2)
plot(opt_PSD(1), opt_PSD(2), '*k', 'MarkerSize', 20, 'linewidth', 2)


%plot(2*PSD_proj(1, kPSD), 2*PSD_proj(2, kPSD),':k')

%scatter(DD_proj(1, :),DD_proj(2, :), 'b')
%scatter(CDD_proj(1, :),CDD_proj(2, :), 'r')
legend({'DD', 'DD(E, ?)', 'S_+'}, 'location', 'northwest', 'fontsize', 15)
%title('Random affine cut of structured subsets', 'Fontsize', 16)
hold off
axis square
axis off

function m = v2m(v)
    m = zeros(3);
    m(1, 1) = v(1);
    m(2, 2) = v(2);
    m(3, 3) = v(3);
    m(1, 2) = v(4);
    m(2, 1) = v(4);
    m(2, 3) = v(5);
    m(3, 2) = v(5);    
end

function v = m2v(m)
    v = zeros(5, 1);
    v(1) = m(1, 1); 
    v(2) = m(2, 2);
    v(3) = m(3 ,3);
    v(4) = m(1, 2); 
    v(5) = m(2, 3);    
end