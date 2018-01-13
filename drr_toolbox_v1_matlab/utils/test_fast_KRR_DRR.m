

function yp2 = test_fast_KRR_DRR(dat_test,models)

 
K_folds = length(models);
aux_yp = zeros(size(dat_test,1),K_folds);

for kk = 1:K_folds
    aux_yp(:,kk) = test_KRR(models(kk).model,dat_test);
end
    
yp2 = mean(aux_yp,2);
    

