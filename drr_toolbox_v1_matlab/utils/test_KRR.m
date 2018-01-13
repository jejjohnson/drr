

function [Ypred_krr] = test_KRR(models,Xtest)


% Centrar datos
% if size(Xtest,1)>1
%     Xtest = scale(Xtest);
% end

Ktest = kernelmatrix('rbf',Xtest',models.X',models.krr_sigma1);
Ypred_krr = Ktest*models.krr_gamma;
