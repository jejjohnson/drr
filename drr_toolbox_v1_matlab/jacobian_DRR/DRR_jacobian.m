function J = DRR_jacobian(x,delta,Model,norm)

% DRR_jacobian computes the Jacobian matrix of DRR at the point x (d*1 vector)
%
% J = DRR_jacobian(x,delta_x,Model,norm)
% 
%     x       = point where the Jacobian will be computed (d*1 vector)
%     delta_x = size of increment in each dimension for the numerical partial derivative
%     Model   = DRR model obtained using DRR_method.m or DRR_normaliz.m
%     norm    = 0 or 1 to choose between regular and normalized DRR
%

d = length(x);

J = zeros(d,d);

if norm ==1
   y_x = apply_DRR_normaliz(x,Model);
else
   y_x = apply_DRR_method(x,Model);    
end

for i = 1:d
    for j = 1:d
          
        delta_x = zeros(d,1);
        delta_x(j) = delta;

        if norm ==1
           y_xmdx = apply_DRR_normaliz(x+delta_x,Model);
        else
           y_xmdx = apply_DRR_method(x+delta_x,Model);            
        end

        J(i,j) = (y_xmdx(i) - y_x(i))/delta;
    end
end



