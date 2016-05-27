function [f] = aviprada_hw5(N, L)
%Clear screen
    clc;
%Define mesh    
    mesh.x = linspace(0, L, N);    
    mesh.conn = [1:N-1; 2:N];
%Define body force    
    b = @(x) 5*sin(2*pi*x/L);
    f = zeros(N, 1);
    quad_pts = [[-1, 1]/sqrt(3); 1, 1]; 
%Evaluate the force vector   
    for c = mesh.conn
        xe = mesh.x(c);
        for q = quad_pts
            S = shape(q);
            dNdp = gradshape(q);
            Je = xe * dNdp;
            dNdx = dNdp / Je;
            x = xe*S;
            f(c) = f(c) + S*b(x) * Je*q(2);
        end
    end
%f = 10N at x = 0    
    f(1) = f(1) - 10;
end

%Evaluate shape function
function [S] = shape(q)
    S = [0.5*(1-q(1)); 0.5*(1+q(1))];
end

%Evaluate the gradient of shape function
function [dNdp] = gradshape(q)
    dNdp = [-0.5; 0.5];
end