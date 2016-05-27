function [] = aviprada_hw4(n,gq)
%clear screen
    clc;
    k = 1;
for i=n:5:2000
%assign x-coordinates to the number of nodes    
    x = linspace(0, 10, i);
%calculate the corresponding y-coordinates for the given function    
    y = cos(x);
%create the connectivity matrix    
    conn = [1:i-1;2:i];
%initialize the summation for gauss quadrature    
    sum1 = 0;
    sum2 = 0;
    
%gauss quadrature
if gq == 1 
    quad_pts = [0; 2]; 
elseif gq == 2 
    quad_pts = [[-1, 1]/sqrt(3); 1, 1]; 
elseif gq == 3 
    quad_pts = [[1,0,-1]*sqrt(3/5); [5,8,5]/9]; 
elseif gq == 4 
    x1 = sqrt(3/7 - 2/7*sqrt(6/5)); 
    x2 = sqrt(3/7 + 2/7*sqrt(6/5)); 
    w1 = (18+sqrt(30)) / 36; 
    w2 = (18-sqrt(30)) / 36; 
    quad_pts = [x1,-x1, x2, -x2; w1, w1, w2, w2]; 
end

%perform gauss quadrature on all the elements
    for c = conn
%extract the nodal x and y coordinates and create jacobian        
        xe = x(1,c);
        ye = y(1,c);
        le = xe(1,2) - xe(1,1);
        %le = sqrt((xe(1,2) - xe(1,1))^2 + (ye(1,2) - ye(1,1))^2);
%take derivative of the shape functions        
        dNxdx = [-1/le;1/le];
%take 1 quadrature point at a time        
        for q = quad_pts
%create the shape function matrix in parent coordinates            
            N = 0.5*[1-q(1); 1+q(1)];
%convert the parent coordinates to global coordinates            
            xc = xe*N;
%create the shape function matrix in global coordinates            
            Nx = [(xc-xe(1,2))/(-1*le);(xc-xe(1,1))/le];
%calculate the analytical and numerical values    
            fa1 = cos(xc);
            fa2 = -1*sin(xc);
            fn1 = ye*Nx;
            fn2 = ye*dNxdx;
%calculate errors e1 and e2            
            sum1 = sum1 + (fa1 - fn1)^2*(le/2)*q(2);
            sum2 = sum2 + (fa2 - fn2)^2*(le/2)*q(2);
        end
    end
    e1(k) = sqrt(sum1);
    e2(k) = sqrt(sum2);
    steps(k) = le;
    k = k + 1;
end
%plot and analysis
%rate of convergence is the line slope
p1 = polyfit(log(steps), log(e1), 1);
p2 = polyfit(log(steps), log(e2), 1);

figure();
loglog(steps, e1, '-r*');
s1 = sprintf('Log-log plot of displacement error vs element step size (Rate of convergence: %0.3f)',p1(1));
title(s1);
xlabel('Step size');
ylabel('L2 norm of displacement error');
grid on;

figure();
loglog(steps, e2, '-bo');
s2 = sprintf('Log-log plot of strain error vs element step size (Rate of convergence: %0.3f)',p2(1));
title(s2);
xlabel('Step size');
ylabel('L2 norm of strain error');
grid on;
end