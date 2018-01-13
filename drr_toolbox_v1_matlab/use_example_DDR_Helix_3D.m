

clear all
close all


% addpath(genpath('./DDR_toolbox/'))


%% DATA: NOISY HELIX

a = 3;
b = 2;
t = linspace(0,4*pi,1000);
x = a*cos(t);
y = a*sin(t);
z = b*t;

helix = [x; y; z];

ruido1 = [];
ruido2 = [];
sigma = linspace(0.1,1.8,length(x));
for i=1:length(x)
    ruido1 = [ruido1 sigma(i)*randn(3,1)]; 
    ruido2 = [ruido2 sigma(i)*randn(3,1)]; 
end

dat = helix + ruido1;
dat2 = helix + ruido2;

% show data
figure
plot3(dat(1,:),dat(2,:),dat(3,:),'b.')
axis equal

%% TRAINING

[datT Model] = DRR_method(dat);
%[datT Model] = DRR_normaliz(dat);

%% APPLY DRR to new data

datT2 = apply_DRR_method(dat2,Model);
%datT2 = apply_DRR_normaliz(dat2,Model);


%% Transformed data

figure
plot3(datT(1,:),datT(2,:),datT(3,:),'b.')
axis equal
title(['DRR norm, Transformed data, Train'])

%% 

figure
plot3(datT2(1,:),datT2(2,:),datT2(3,:),'b.')
axis equal
title(['DRR norm, Transformed data, Test'])


%% DRR to 1D

idat = inv_DRR_method([datT(1,:); zeros(size(datT,1)-1,size(datT,2))],Model);
idat2 = inv_DRR_method([datT2(1,:); zeros(size(datT2,1)-1,size(datT2,2))],Model);

%idat = inv_DRR_normaliz([datT(1,:); zeros(size(datT,1)-1,size(datT,2))],Model);
%idat2 = inv_DRR_normaliz([datT2(1,:); zeros(size(datT2,1)-1,size(datT2,2))],Model);


% show results

figure
plot3(dat(1,:),dat(2,:),dat(3,:),'b.')
axis equal
hold on
plot3(idat(1,:),idat(2,:),idat(3,:),'r.')

title(['DRR, 1D reduction, Train'])

figure
plot3(dat(1,:),dat(2,:),dat(3,:),'b.')
axis equal
hold on
plot3(idat2(1,:),idat2(2,:),idat2(3,:),'r.')

title(['DRR, 1D reduction, Test'])

%% DRR to 2D

idat = inv_DRR_method([datT(1:2,:); zeros(size(datT,1)-2,size(datT,2))],Model);
idat2 = inv_DRR_method([datT2(1:2,:); zeros(size(datT2,1)-2,size(datT2,2))],Model);

%idat = inv_DRR_normaliz([datT(1:2,:); zeros(size(datT,1)-2,size(datT,2))],Model);
%idat2 = inv_DRR_normaliz([datT2(1:2,:); zeros(size(datT2,1)-2,size(datT2,2))],Model);

%%
figure
plot3(dat(1,:),dat(2,:),dat(3,:),'b.')
axis equal
hold on
plot3(idat(1,:),idat(2,:),idat(3,:),'r.')

title(['DRR, 2D reduction, Train'])



%%
figure
plot3(dat(1,:),dat(2,:),dat(3,:),'b.')
axis equal
hold on
plot3(idat2(1,:),idat2(2,:),idat2(3,:),'r.')

title(['DRR, 2D reduction, Test'])

