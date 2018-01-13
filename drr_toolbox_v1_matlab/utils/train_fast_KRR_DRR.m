

function models = train_fast_KRR_DRR(dat_train1,Y1,dat_train2,Y2)


[NN DD] = size(dat_train1); 
K_folds = floor(NN/2000)+1;
aa = crossvalind('Kfold',size(dat_train1,1),K_folds);

for kk = 1:K_folds
    models(kk).model = train_KRR(dat_train1(find(aa==kk),:),Y1(find(aa==kk)),dat_train2(find(aa==kk),:),Y2(find(aa==kk)));
end
