function [data_x,data_r,raspas_x,raspas_r,struct_ppa]=genera_casquete(angulosE1,polin1,anguloE2,polin2,lim_alpha1,lim_alpha2,disp_3rd_dim,N,std_noise,num_alpha1,num_alpha2);

% GENERA_CASQUETE generates a 2d noisy curved manifold embedded in a 3d space
% by inverting a PPA model
% 
% Parameters of the PPA model 
%
%    E1 = 3*3 rotation matrix defined by the angles:
%              psi   -around x3-
%              theta -around x2-
%              phi   -around x1-
%
%    W1 = 2*3 coeffcient matrix for the first polynomial
%              example W1=[0 0 0;0 0 0.08]; 
%
%    E2 = 2*2 rotation matrix defined by the angle:
%              gamma  -around z1-
%
%    W2 = 1*3 coeffcient matrix for the first polynomial
%              example W1=[0 0 0.05];
%
% Range of the 2d surface in the PPA domain
%
%   lim_alpha1 = [min1 max1]
%   lim_alpha2 = [min2 max2]
%
% Displacement in the 3rd dimension in the PPA domain
%
%   disp_3rd_dim = 
%
% Number of points =
%
% std_noise =
%
% 2d reticle in the PPA domain
%
%    num_alpha1 =
%    num_alpha2 =
%

%% PPA MODEL

    psi=angulosE1(1);    % Eje x3
    tecta=angulosE1(2);  % Eje x2
    fi=angulosE1(3);     % Eje x1

    E1=[cos(tecta)*cos(psi) -cos(fi)*sin(psi)+sin(fi)*sin(tecta)*cos(psi)  sin(fi)*sin(psi)+cos(fi)*sin(tecta)*cos(psi);
        cos(tecta)*sin(psi)  cos(fi)*cos(psi)+sin(fi)*sin(tecta)*sin(psi) -sin(fi)*cos(psi)+cos(fi)*sin(tecta)*sin(psi);
        -sin(tecta)                  sin(fi)*cos(tecta)                    cos(fi)*cos(tecta)];
    
    E2=[cos(anguloE2) -sin(anguloE2);sin(anguloE2) cos(anguloE2)];

    % NON-LINEAR

    struct_ppa(1).Res.m=[0 0 0];
    struct_ppa(1).Res.V=E1';
    struct_ppa(1).Res.W=polin1'; 

    struct_ppa(2).Res.m=[0 0];
    struct_ppa(2).Res.V=E2';
    struct_ppa(2).Res.W=polin2';
    
    struct_ppa(3).Res=0;
   
%% DATA IN THE PPA DOMAIN

% lim_alpha1,lim_alpha2,disp_3rd_dim,N,std_noise,num_alpha1,num_alpha2);

r = [lim_alpha1(1)+rand(1,N)*(lim_alpha1(2)-lim_alpha1(1));lim_alpha2(1)+rand(1,N)*(lim_alpha2(2)-lim_alpha2(1));disp_3rd_dim*ones(1,N)];

data_r = r + std_noise*randn(3,N);

raspas_r = [];
alfas_2=linspace(lim_alpha2(1),lim_alpha2(2),num_alpha2);
for i=1:num_alpha2
    raspas_r = [raspas_r [linspace(lim_alpha1(1),lim_alpha1(2),num_alpha1);alfas_2(i)*ones(1,num_alpha1);disp_3rd_dim*ones(1,num_alpha1)]];
end

raspas_x=inv_PPA(raspas_r,struct_ppa);
data_x=inv_PPA(data_r,struct_ppa);    
    
%     
%     
% 
%     % Raspas 1
% 
%     raspas_1=[];
%     val_alfa2 = linspace(-8,8,10);
% 
%     for i=1:10
% 
%         alfa2 = val_alfa2(i)*ones(1,100);   % Raspa 1 esta contenida en el plano alfa2=0;
%         x2 = zeros(1,100);      % No noise!
%         m2 = W2*[ones(1,100);alfa2;alfa2.^2/2];
% 
%         alfa1=linspace(-10,10,100);
%         x1 = E2'*[alfa2;x2+m2];
%         m1 = W1*[ones(1,100);alfa1;alfa1.^2/2];
% 
%         x0 = E1'*[alfa1;x1+m1];
% 
%         raspa_1 = x0;
%         raspas_1 = [raspas_1 raspa_1];
% 
%     end
% 
%     % Raspas 2
% 
%     val_alfa1=linspace(-8,8,7);
% 
%     raspas_2=[];
% 
%     for i=1:7
% 
%     alfa2 = linspace(-10,10,50);   % Cada Raspa 2 recorre el rango de alfa2; 
%     x2 = zeros(1,50);      % No noise!
%     m2 = W2*[ones(1,50);alfa2;alfa2.^2/2];
% 
%     alfa1=val_alfa1(i)*ones(1,50);
%     x1 = E2'*[alfa2;x2+m2];
%     m1 = W1*[ones(1,50);alfa1;alfa1.^2/2];
% 
%     x0 = E1'*[alfa1;x1+m1];
% 
%     raspa_2 = x0;
%     raspas_2 = [raspas_2 raspa_2];
% 
%     %hold on,plot3(raspa_2(1,:),raspa_2(2,:),raspa_2(3,:),'b.')
% 
%     end
%     %xlabel('x_1'),ylabel('x_2'),zlabel('x_3'),
%     %box on
%     %axis equal
% 
%     % Datos:
%     
%     x_dat = raspas_1 + std_noise*randn(3,length(raspas_1(1,:)));
%     x_dat2 = raspas_2 + std_noise*randn(3,length(raspas_2(1,:)));
%     
%     x_dat = [x_dat x_dat2];