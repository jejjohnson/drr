

function models = train_KRR(Xtrain,Ytrain,Xvalid,Yvalid)

warning off

% Centrar datos
% if size(Xtrain,1)>1
%     Xtrain = scale(Xtrain);
%     Xvalid = scale(Xvalid);
% end

n = size(Ytrain,1);
n2 = size(Yvalid,1);

%  KRR
Ntrains  = 10;

LAMBDAS1 = [0 logspace(-2,1,Ntrains)];
SIGMAS1  = logspace(-2,3,Ntrains);


i=0;
for sigma1 = SIGMAS1
    Ktrain_x = kernelmatrix('rbf',Xtrain',Xtrain',sigma1);
    Kvalid   = kernelmatrix('rbf',Xvalid',Xtrain',sigma1);
    
    for lambda1 = LAMBDAS1
        i=i+1;
        
        gamma_rbf = (Ktrain_x + lambda1*eye(n))\Ytrain;
        
        Ypred_rbf = Kvalid*gamma_rbf;
        
        res_rbf = sum(sum(abs(Yvalid-Ypred_rbf)));
        
        res(i,:) = [sigma1 lambda1 res_rbf];
    end
end


% Best KRR

[valor idx] = min(res(:,3));
models.krr_sigma1  = res(idx,1);
models.krr_lambda1 = res(idx,2);


models.X = [Xtrain; Xvalid];

Ktrain_x = kernelmatrix('rbf',models.X',models.X',models.krr_sigma1);
models.krr_gamma = (Ktrain_x + models.krr_lambda1*eye(size(Ktrain_x,1))) \ [Ytrain; Yvalid];

