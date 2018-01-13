

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


function inv_dat = inv_DRR_normaliz(datT2,model,DIM)


dim_ori = size(model(1).V,1);

if ~exist('DIM')
   DIM = size(datT2,1);
end

datT2 = [datT2(1:DIM,:) ; zeros(dim_ori-DIM,size(datT2,2))];

aux_dat2 = datT2(1,:);


for n_dim = 1:size(datT2,1)-1
    
    [sigma_pred_krr] = test_fast_KRR_DRR(aux_dat2',model(end-n_dim+1).models_sigma);
    [Ypred_krr] = test_fast_KRR_DRR(aux_dat2',model(end-n_dim+1).models_coef);
    
    aux_dat2 = [aux_dat2; Ypred_krr' + sigma_pred_krr'.*datT2(n_dim+1,:)];
    
end


inv_dat = inv(model(1).V')*aux_dat2;

inv_dat = inv_dat + repmat(model(1).mm,1,size(inv_dat,2));


