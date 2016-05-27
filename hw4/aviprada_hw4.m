function [e1, e2] = aviprada_hw4(n)
%clear screen
    clc;
%assign x-coordinates to the number of nodes    
    x = linspace(0, 10, n);
%calculate the corresponding y-coordinates for the given function    
    y = cos(x);
%create the connectivity matrix    
    conn = [1:n-1;2:n];
%initialize the summation for gauss quadrature    
    sum1 = 0;
    sum2 = 0;
%assume a 3-point gauss quadrature    
    quad_pts = [[-1,1]/sqrt(3);1,1];
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
            xc = xe*N
%create the shape function matrix in global coordinates            
            Nx = [(xc-xe(1,2))/(-1*le);(xc-xe(1,1))/le];
%calculate the analytical and numerical values    
            fa1 = cos(xc);
            fa2 = -1*sin(xc)
            fn1 = ye*Nx
            fn2 = ye*dNxdx
%calculate errors e1 and e2            
            sum1 = sum1 + (fa1 - fn1)^2*(le/2)*q(2);
            sum2 = sum2 + (fa2 - fn2)^2*(le/2)*q(2);
        end
    end
    e1 = sqrt(sum1);
    e2 = sqrt(sum2);
%display results    
    e1
    e2
end