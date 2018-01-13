%
% This script generates manifolds (inspired in Delicado01) that fulfil these 
% two restrictions to a bigger or lesser degree:
% (1) Conditional Mean Independence assumption
% (2) Stationarity along the Principal Curve
%
% Data:                       x_data
% First Principal Curve:      x_first
% Secondary Principal Curves: x_secondary

% DATA with no torsion

N = 10000;
factor = 1;
theta_max = factor*pi/2;
ro = 1;
sigma_1 = factor*0.25;
sigma_2 = factor*0.01;  % Noise in the 3rd principal direction
A = factor*0.15; 
slope = 1/factor;     % Slope of the torsion angle

theta = linspace(0,theta_max,N);
theta_n = theta + (10*theta_max/N)*randn(1,N);

n = rand(1,N);
xx = (ro-sigma_1*(2*n-1)).*cos(theta_n);
yy = (ro-sigma_1*(2*n-1)).*sin(theta_n);
zz = A*sin((2*n-1)*6*pi/6) + sigma_2*randn(1,N);

x = [xx;yy;zz];

cada = 5;
figure,plot3(x(1,1:cada:end),x(2,1:cada:end),x(3,1:cada:end),'b.'),axis equal

% First Principal Curve

xx = (ro).*cos(theta);
yy = (ro).*sin(theta);
zz = 0*theta;

x_first = [xx;yy;zz];

% Secondary principal curves

x_second = [];
theta_ras = [];
theta_raspa = linspace(0.075*theta_max,0.925*theta_max,10);
x_secondary_guay= [];
theta_ras_guay = [];

rr = linspace(-1,1,71);
for i=1:length(theta_raspa)
    xx = (ro-sigma_1*rr).*cos(theta_raspa(i));
    yy = (ro-sigma_1*rr).*sin(theta_raspa(i));
    zz = A*sin(rr*4*pi/4);
    x_second = [x_second [xx;yy;zz]];
    theta_ras = [theta_ras theta_raspa(i)*ones(1,length(rr))];
 
        x_secondary_guay = [x_secondary_guay [xx(6:66);yy(6:66);zz(6:66)]];
        theta_ras_guay = [theta_ras_guay theta_raspa(i)*ones(1,length(xx(1,6:66)))];

end

hold on,
plot3(x_first(1,1:cada:end),x_first(2,1:cada:end),x_first(3,1:cada:end),'r.')
plot3(x_second(1,:),x_second(2,:),x_second(3,:),'r.')


% Torsion along the principal curve

xtor_dat = [];
for i=1:length(theta)
    
    x0 = x_first(:,i);
    delta_x = x(:,i) - x0;
    R_mtheta = [cos(-theta(i)) -sin(-theta(i)) 0;sin(-theta(i)) cos(-theta(i)) 0;0 0 1];
    R_theta = [cos(theta(i)) -sin(theta(i)) 0;sin(theta(i)) cos(theta(i)) 0;0 0 1];
    R_alpha = [cos(slope*theta(i)) 0 -sin(slope*theta(i));0 1 0;sin(slope*theta(i)) 0 cos(slope*theta(i))];
    
    xtor_dat = [xtor_dat R_theta*R_alpha*R_mtheta*delta_x+x0];
    
end

xtor = [];
for i=1:length(theta_ras)
    
    x0 = [ro.*cos(theta_ras(i));ro.*sin(theta_ras(i));0];
    delta_x = x_second(:,i) - x0;
    R_mtheta = [cos(-theta_ras(i)) -sin(-theta_ras(i)) 0;sin(-theta_ras(i)) cos(-theta_ras(i)) 0;0 0 1];
    R_theta = [cos(theta_ras(i)) -sin(theta_ras(i)) 0;sin(theta_ras(i)) cos(theta_ras(i)) 0;0 0 1];
    R_alpha = [cos(slope*theta_ras(i)) 0 -sin(slope*theta_ras(i));0 1 0;sin(slope*theta_ras(i)) 0 cos(slope*theta_ras(i))];
    
    xtor = [xtor R_theta*R_alpha*R_mtheta*delta_x+x0];
    
end

xtor_g = [];
for i=1:length(theta_ras_guay)
    
    x0 = [ro.*cos(theta_ras_guay(i));ro.*sin(theta_ras_guay(i));0];
    delta_x = x_secondary_guay(:,i) - x0;
    R_mtheta = [cos(-theta_ras_guay(i)) -sin(-theta_ras_guay(i)) 0;sin(-theta_ras_guay(i)) cos(-theta_ras_guay(i)) 0;0 0 1];
    R_theta = [cos(theta_ras_guay(i)) -sin(theta_ras_guay(i)) 0;sin(theta_ras_guay(i)) cos(theta_ras_guay(i)) 0;0 0 1];
    R_alpha = [cos(slope*theta_ras_guay(i)) 0 -sin(slope*theta_ras_guay(i));0 1 0;sin(slope*theta_ras_guay(i)) 0 cos(slope*theta_ras_guay(i))];
    
    xtor_g = [xtor_g R_theta*R_alpha*R_mtheta*delta_x+x0];
    
end

x_data = xtor_dat;
x_secondary = xtor;
x_secondary_guay = xtor_g;

figure,plot3(x_data(1,1:cada:end),x_data(2,1:cada:end),x_data(3,1:cada:end),'b.')
hold on,
plot3(x_first(1,1:cada:end),x_first(2,1:cada:end),x_first(3,1:cada:end),'r.')
plot3(x_secondary(1,:),x_secondary(2,:),x_secondary(3,:),'r.')
plot3(x_secondary_guay(1,:),x_secondary_guay(2,:),x_secondary_guay(3,:),'r.')
axis equal

figure,plot3(x_data(1,1:cada:end),x_data(2,1:cada:end),x_data(3,1:cada:end),'b.')
hold on,
plot3(x_first(1,1:cada:end),x_first(2,1:cada:end),x_first(3,1:cada:end),'r.')
plot3(x_secondary_guay(1,:),x_secondary_guay(2,:),x_secondary_guay(3,:),'r.')
axis equal
