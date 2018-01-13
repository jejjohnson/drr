

% clear all
% close all
% 
% input_folder = '/media/disk/vista/Papers/aranyas/JMLR2011/databases/DATABASES_MATS/';
% 
% output_folder = '/media/disk/vista/Papers/aranyas/JMLR2011/DR_results/RESULTS/';
% 
% cd([input_folder])
% 
% load([input_folder 'DATA_EUM'])
% 


function datT2 = apply_DRR_method(dat,model)

datm = dat - repmat(model(1).mm,1,size(dat,2));

%dat2 = model(1).DD*model(1).V'*dat;
dat2 = model(1).V'*datm;

datT2 = zeros(size(dat2));

aux_dat2 = dat2;
Nsamples = size(dat2,2);
for n_dim=1:size(dat2,1)-1
    
    ii = model(n_dim).ii;
    %%
    ii2 = setdiff(1:size(aux_dat2,1),ii);
    
    dat_train1 = aux_dat2(ii2,:);
    dat_train2 = aux_dat2(ii2,:);
    
    Y1 = aux_dat2(ii,:);
    Y2 = aux_dat2(ii,:);
    
    [Ypred_krr] = test_fast_KRR_DRR(dat2(ii2,:)',model(n_dim).models);
    
    aux_dat2 = aux_dat2(ii2,:);
    datT2(end-n_dim+1,:) = dat2(ii,:)-Ypred_krr';
    
end

datT2(1,:) = aux_dat2;


