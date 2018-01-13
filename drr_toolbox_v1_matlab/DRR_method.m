

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


function [datT2 model] = DRR_method(dat,flag_entropy,flag_cross,flag_trans)

if ~exist('flag_entropy')
    flag_entropy = 0;
end

if ~exist('flag_cross')
    flag_cross = 1;
end

if ~exist('flag_trans')
    flag_trans = 0;
end

mm = mean(dat,2);
model.mm = mm;

datm = dat - repmat(mm,1,size(dat,2));

if flag_trans == 0
    
    C = cov(datm');
    [V D] = eig(C);
    d = diag(D);
    [ss ii] = sort(d,'descend');
    DD = diag(ss.^(-0.5));
    
    V = V(:,ii);
    
elseif flag_trans == 1
    % ICA
    aux_datm = datm;
    for qq=1:size(datm,1)-1
        
        % ICA
        [A W] = fastica(aux_datm,'approach','symm','verbose','off');
        
        Wn = W / (abs(det(W))^(1/size(W,1)));
        %Wn = my_GramS(W);
        
        dat2 = Wn*aux_datm;
        Nsamples = size(dat2,2);
        hx = [];
        for n = 1:size(dat2,1)
            [p R]=hist(dat2(n,:),sqrt(Nsamples));
            delta = R(3)-R(2);
            hx(n)=entropy_mm(p)+log2(delta);
        end
        
        [ss ii] = sort(hx,'ascend');
        
        VV = my_GramS(Wn(ii,:)');
        
        VV = fliplr(VV)';
        
        lala = VV*aux_datm;
        
        aux_datm = lala(1:size(datm,1)-qq,:);
        
        if qq == 1
            Vant = VV;
        else
            VV(end+1,end+1) = 1;
            Vant = VV*Vant;
        end
        
    end
    
    V = Vant';
elseif flag_trans == 2
    
     % MCV
    aux_datm = datm;
    for qq=1:size(datm,1)-1
        
        % MCV
        [auxW W] = minimum_cond_var_direct(aux_datm);
       
        
        VV = (auxW); %flipud
        
        lala = VV*aux_datm;
        
        aux_datm = lala(1:size(datm,1)-qq,:);
        
        if qq == 1
            Vant = VV;
        else
            VV(end+1,end+1) = 1;
            Vant = VV*Vant;
        end
        
    end
    
    V = Vant';
    
end

dat2 = V'*datm;

model.V = V;


if flag_entropy == 1
    %% Conditional entropies
    % MI data
    tic
    [datT,Trans] = MI_auto_RBIG_2(dat2,1000);
    toc
    CC = cumsum(cat(1,Trans.I));
    I_x1_xn = CC(end);
end

datT2 = zeros(size(dat2));

aux_dat2 = dat2;
Nsamples = size(dat2,2);
for n_dim=1:size(dat2,1)-1
    tic
    if flag_entropy == 1
        %% Entropies
        
        % marginal
        h_dat2 = [];
        for n=1:size(aux_dat2,1)
            [p R] = hist(aux_dat2(n,:),round(sqrt(Nsamples)));
            delta = R(3)-R(2);
            h_dat2(n) = entropy_mm(p)+log2(delta);
        end
        
        % Conditionals
        
        h_y_dd_x1_x2 = [];
        for n=1:size(aux_dat2,1)
            
            ii = setdiff(1:size(aux_dat2,1),n);
            
            tic
            [datT,Trans] = MI_auto_RBIG_2(aux_dat2(ii,:),200);
            [n toc]
            CC = cumsum(cat(1,Trans.I));
            
            I_x2_xn = CC(end);
            
            h_y_dd_x1_x2(n) = I_x2_xn - I_x1_xn + h_dat2(n);
        end
        
        [mm ii] = min(h_y_dd_x1_x2);
    else
        ii = size(dat2,1)-n_dim+1;
    end
    %%
    ii2 = setdiff(1:size(aux_dat2,1),ii);
    
    if flag_cross == 1
        NN_t = round(size(aux_dat2,2)/2);
        rand('seed',1)
        ii_2 = randperm(size(aux_dat2,2));
        
        dat_train1 = aux_dat2(ii2,ii_2(1:NN_t));
        dat_train2 = aux_dat2(ii2,ii_2(NN_t+1:end));
        
        Y1 = aux_dat2(ii,ii_2(1:NN_t));
        Y2 = aux_dat2(ii,ii_2(NN_t+1:end));
    else
        dat_train1 = aux_dat2(ii2,:);
        dat_train2 = aux_dat2(ii2,:);
        
        Y1 = aux_dat2(ii,:);
        Y2 = aux_dat2(ii,:);
    end
    
    % models = train_KRR(dat_train1',Y1',dat_train2',Y2');
    models = train_fast_KRR_DRR(dat_train1',Y1',dat_train2',Y2');
    
    Ypred_krr = test_fast_KRR_DRR(dat2(ii2,:)',models);
    
    
    model(n_dim).models = models;
    model(n_dim).ii = ii;
    
    aux_dat2 = aux_dat2(ii2,:);
    datT2(end-n_dim+1,:) = dat2(ii,:)-Ypred_krr';
    [n_dim toc]
end

datT2(1,:) = aux_dat2;



function H = entropy_mm(p)

% mle estimator with miller-maddow correction

c = 0.5 * (sum(p>0)-1)/sum(p);  % miller maddow correction
p = p/sum(p);               % empirical estimate of the distribution
idx = p~=0;
H = -sum(p(idx).*log2(p(idx))) + c;     % plug-in estimator of the entropy with correction

