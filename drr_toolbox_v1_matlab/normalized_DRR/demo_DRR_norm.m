% NORMALIZED DRR FOR HETEROSCEDASTIC DATA
% 
%

% DATA
%%%%%%%%%%

% [x_facil,xff,xsf,x_chungo,xfc,xsc] = data_delicado_heteroscedastic;
% 
% figure,plot3(x_facil(1,1:10:end),x_facil(2,1:10:end),x_facil(3,1:10:end),'b.')
% hold on,plot3(xff(1,:),xff(2,:),xff(3,:),'r.')
% hold on,plot3(xsf(1,:),xsf(2,:),xsf(3,:),'r.')
% axis equal 
% 
% figure,plot3(x_chungo(1,1:10:end),x_chungo(2,1:10:end),x_chungo(3,1:10:end),'b.')
% hold on,plot3(xfc(1,:),xfc(2,:),xfc(3,:),'r.')
% hold on,plot3(xsc(1,:),xsc(2,:),xsc(3,:),'r.')
% axis equal

[x_facil,x_first] = helix_heteroscedastic;

figure,plot3(x_facil(1,1:10:end),x_facil(2,1:10:end),x_facil(3,1:10:end),'b.')
hold on,plot3(x_first(1,1:500:end),x_first(2,1:500:end),x_first(3,1:500:end),'r.')
xlabel('x_1'),ylabel('x_2'),zlabel('x_3'),
axis equal
title('Heteroscedastic Manifold')

% Train DRR
%%%%%%%%%%%%%%%%%

[datT Model] = DRR_normaliz(x_facil(:,1:10:end));
[datT_c Model_c] = DRR_method(x_facil(:,1:10:end));


% Apply DRR to new data
%%%%%%%%%%%%%%%%%%%%%%%%%%%

datT2 = apply_DRR_normaliz(x_facil(:,2:10:end),Model);
datT2_c = apply_DRR_method(x_facil(:,2:10:end),Model_c);

%% Transformed data
%%%%%%%%%%%%%%%%%%%%%

figure
plot3(datT(1,:),datT(2,:),datT(3,:),'b.')
axis equal
title(['DRR (heteroscedastic), Transformed data, Train'])

D_n = diag(std(datT'))

figure
plot3(datT_c(1,:),datT_c(2,:),datT_c(3,:),'b.')
axis equal
title(['DRR (symplectic), Transformed data, Train'])

D = diag(std(datT_c'))

%% 

figure
plot3(datT2(1,:),datT2(2,:),datT2(3,:),'b.')
axis equal
title(['DRR norm, Transformed data, Test'])

figure
plot3(datT2_c(1,:),datT2_c(2,:),datT2_c(3,:),'b.')
axis equal
title(['DRR, Transformed data, Test'])

%% Jacobian and Metric features for the helix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xf = x_first(:,1:500:end);

factJ = 0.75;
factM = 0.75;
factn = 0.75;

D(1,1)=0.1*D(1,1);
D_n(1,1)=0.1*D_n(1,1);

figure(100),plot3(x_facil(1,1:10:end),x_facil(2,1:10:end),x_facil(3,1:10:end),'c.')
xlabel('x_1'),ylabel('x_2'),zlabel('x_3'),
axis equal
title('Features from \nabla R^{-1} (symplectic)')

figure(101),plot3(x_facil(1,1:10:end),x_facil(2,1:10:end),x_facil(3,1:10:end),'c.')
xlabel('x_1'),ylabel('x_2'),zlabel('x_3'),
title('Features from \nabla R^{T}.\nabla R (symplectic)')
axis equal

figure(102),plot3(x_facil(1,1:10:end),x_facil(2,1:10:end),x_facil(3,1:10:end),'c.')
xlabel('x_1'),ylabel('x_2'),zlabel('x_3'),
axis equal
title('Features from \nabla R^{-1} (heteroscedastic)')

figure(103),plot3(x_facil(1,1:10:end),x_facil(2,1:10:end),x_facil(3,1:10:end),'c.')
xlabel('x_1'),ylabel('x_2'),zlabel('x_3'),
title('Features from \nabla R^{T}.\nabla R (heteroscedastic)')
axis equal


for i=1:length(xf(1,:))
    
    J = DRR_jacobian(xf(:,i),0.001,Model,1);
    scale = 1; %(det(J))^(-1/3)
    % normas = sqrt(sum(iJ.^2));
    % iJ = [iJ(:,1)/normas(1) iJ(:,2)/normas(2) iJ(:,3)/normas(3)];
    Mn = J'*inv(D_n*D_n)*J;
    [Bn,Ln]=eigs(inv(Mn));
    Bn = Bn*(Ln.^(1/2)); % This would be the right scaling but I am going
                          % to increase the small vectors for the sake of visibility
    iJ = inv(inv(D_n)*J);
    %Bn = Bn*(Ln.^(1/4));
    %iJ = inv(inv(D_n.^(1/2))*J);
    
    Jc = DRR_jacobian(xf(:,i),0.001,Model_c,0);
    M = Jc'*inv(D*D)*Jc;
    [B,L]=eigs(inv(M));
     B = B*(L.^(1/2));
    iJc = inv(inv(D)*Jc);
    %B = B*(L.^(1/4));
    %iJc = inv(inv(D.^(1/2))*Jc);

    figure(100),
    hold on,plot3([xf(1,i) xf(1,i)+factJ*iJc(1,1)],[xf(2,i) xf(2,i)+factJ*iJc(2,1)],[xf(3,i) xf(3,i)+factJ*iJc(3,1)],'r-','linewidth',3)
    hold on,plot3([xf(1,i) xf(1,i)+factJ*iJc(1,2)],[xf(2,i) xf(2,i)+factJ*iJc(2,2)],[xf(3,i) xf(3,i)+factJ*iJc(3,2)],'g-','linewidth',3,'color',[0 0.5 0])
    hold on,plot3([xf(1,i) xf(1,i)+factJ*iJc(1,3)],[xf(2,i) xf(2,i)+factJ*iJc(2,3)],[xf(3,i) xf(3,i)+factJ*iJc(3,3)],'b-','linewidth',3)
    hold on,plot3([xf(1,i) xf(1,i)-factJ*iJc(1,1)],[xf(2,i) xf(2,i)-factJ*iJc(2,1)],[xf(3,i) xf(3,i)-factJ*iJc(3,1)],'r-','linewidth',3)
    hold on,plot3([xf(1,i) xf(1,i)-factJ*iJc(1,2)],[xf(2,i) xf(2,i)-factJ*iJc(2,2)],[xf(3,i) xf(3,i)-factJ*iJc(3,2)],'g-','linewidth',3,'color',[0 0.5 0])
    hold on,plot3([xf(1,i) xf(1,i)-factJ*iJc(1,3)],[xf(2,i) xf(2,i)-factJ*iJc(2,3)],[xf(3,i) xf(3,i)-factJ*iJc(3,3)],'b-','linewidth',3)
    
    figure(101),
    hold on,plot3([xf(1,i) xf(1,i)+factM*B(1,1)],[xf(2,i) xf(2,i)+factM*B(2,1)],[xf(3,i) xf(3,i)+factM*B(3,1)],'r-','linewidth',3)
    hold on,plot3([xf(1,i) xf(1,i)+factM*B(1,2)],[xf(2,i) xf(2,i)+factM*B(2,2)],[xf(3,i) xf(3,i)+factM*B(3,2)],'g-','linewidth',3,'color',[0 0.5 0])
    hold on,plot3([xf(1,i) xf(1,i)+factM*B(1,3)],[xf(2,i) xf(2,i)+factM*B(2,3)],[xf(3,i) xf(3,i)+factM*B(3,3)],'b-','linewidth',3)
    hold on,plot3([xf(1,i) xf(1,i)-factM*B(1,1)],[xf(2,i) xf(2,i)-factM*B(2,1)],[xf(3,i) xf(3,i)-factM*B(3,1)],'r-','linewidth',3)
    hold on,plot3([xf(1,i) xf(1,i)-factM*B(1,2)],[xf(2,i) xf(2,i)-factM*B(2,2)],[xf(3,i) xf(3,i)-factM*B(3,2)],'g-','linewidth',3,'color',[0 0.5 0])
    hold on,plot3([xf(1,i) xf(1,i)-factM*B(1,3)],[xf(2,i) xf(2,i)-factM*B(2,3)],[xf(3,i) xf(3,i)-factM*B(3,3)],'b-','linewidth',3)

    figure(102),
    hold on,plot3([xf(1,i) xf(1,i)+factn*iJ(1,1)],[xf(2,i) xf(2,i)+factn*iJ(2,1)],[xf(3,i) xf(3,i)+factn*iJ(3,1)],'r-','linewidth',3)%,'color',[1 0.5 0])
    hold on,plot3([xf(1,i) xf(1,i)+factn*iJ(1,2)],[xf(2,i) xf(2,i)+factn*iJ(2,2)],[xf(3,i) xf(3,i)+factn*iJ(3,2)],'g-','linewidth',3,'color',[0 0.5 0])
    hold on,plot3([xf(1,i) xf(1,i)+factn*iJ(1,3)],[xf(2,i) xf(2,i)+factn*iJ(2,3)],[xf(3,i) xf(3,i)+factn*iJ(3,3)],'b-','linewidth',3)%,'color',[0 0.3 1])
    hold on,plot3([xf(1,i) xf(1,i)-factn*iJ(1,1)],[xf(2,i) xf(2,i)-factn*iJ(2,1)],[xf(3,i) xf(3,i)-factn*iJ(3,1)],'r-','linewidth',3)%,'color',[1 0.5 0])
    hold on,plot3([xf(1,i) xf(1,i)-factn*iJ(1,2)],[xf(2,i) xf(2,i)-factn*iJ(2,2)],[xf(3,i) xf(3,i)-factn*iJ(3,2)],'g-','linewidth',3,'color',[0 0.5 0])
    hold on,plot3([xf(1,i) xf(1,i)-factn*iJ(1,3)],[xf(2,i) xf(2,i)-factn*iJ(2,3)],[xf(3,i) xf(3,i)-factn*iJ(3,3)],'b-','linewidth',3)%,'color',[0 0.3 1])
    
    figure(103),
    hold on,plot3([xf(1,i) xf(1,i)+factn*scale*Bn(1,1)],[xf(2,i) xf(2,i)+factn*scale*Bn(2,1)],[xf(3,i) xf(3,i)+factn*scale*Bn(3,1)],'r-','linewidth',3)%,'color',[1 0.5 0])
    hold on,plot3([xf(1,i) xf(1,i)+factn*scale*Bn(1,2)],[xf(2,i) xf(2,i)+factn*scale*Bn(2,2)],[xf(3,i) xf(3,i)+factn*scale*Bn(3,2)],'g-','linewidth',3,'color',[0 0.5 0])
    hold on,plot3([xf(1,i) xf(1,i)+factn*scale*Bn(1,3)],[xf(2,i) xf(2,i)+factn*scale*Bn(2,3)],[xf(3,i) xf(3,i)+factn*scale*Bn(3,3)],'b-','linewidth',3)%,'color',[0 0.3 1])
    hold on,plot3([xf(1,i) xf(1,i)-factn*scale*Bn(1,1)],[xf(2,i) xf(2,i)-factn*scale*Bn(2,1)],[xf(3,i) xf(3,i)-factn*scale*Bn(3,1)],'r-','linewidth',3)%,'color',[1 0.5 0])
    hold on,plot3([xf(1,i) xf(1,i)-factn*scale*Bn(1,2)],[xf(2,i) xf(2,i)-factn*scale*Bn(2,2)],[xf(3,i) xf(3,i)-factn*scale*Bn(3,2)],'g-','linewidth',3,'color',[0 0.5 0])
    hold on,plot3([xf(1,i) xf(1,i)-factn*scale*Bn(1,3)],[xf(2,i) xf(2,i)-factn*scale*Bn(2,3)],[xf(3,i) xf(3,i)-factn*scale*Bn(3,3)],'b-','linewidth',3)%,'color',[0 0.3 1])
    
end

figure(1),view([-46 18])
figure(100),view([-46 18])
figure(101),view([-46 18])
figure(102),view([-46 18])
figure(103),view([-46 18])

figure(1),view([-70 10])
figure(100),view([-70 10])
figure(101),view([-70 10])
figure(102),view([-70 10])
figure(103),view([-70 10])
