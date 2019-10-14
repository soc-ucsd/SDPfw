load('quartic_process_3m.mat')

cone = {'sdd', 'dd'};


model = model_agler;
[model_new.At, model_new.b, model_new.c, model_new.K, model_new.info] = ...
decomposed_subset(model.At',model.b,model.c,model.K,cone);
model_new.pars = model.pars;

c = model_new.c;
    
c(1) = 2;
c(2) = 3;
        
[x,y,info] = sedumi(model_new.At, model_new.b, c, model_new.K, model_new.pars);

x_orig = decomposed_recover(x, model_new.info);

K = model.K;

X1 = reshape(x_orig(K.f + (1:K.s(1)^2)), K.s(1), K.s(1));
d1 = eig(X1);
X2 = reshape(x_orig(K.f + K.s(1)^2 +  (1:K.s(2)^2)), K.s(2), K.s(2));
d2 = eig(X2);


