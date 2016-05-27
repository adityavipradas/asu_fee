function [f] = aviprada_hw8(xe, p)
    %Clear screen
    clc;
    %Define the force vector
    f = zeros(length(xe)*2,1);
    %Two-point gauss quadrature
    quad_pts = [[-1, 1]/sqrt(3); 1, 1];
    %Iterate
    for q = quad_pts
        %Shape function
        S = shape(q(1));
        %Gradient
        dNdp = gradshape(q(1));
        %Jacobian
        Je = xe([1,2],[1,8,4]) * dNdp;
        %Update force vector
        f([1,15,7],1) = f([1,15,7],1) + S*p*Je(2)*q(2);
        f([2,16,8],1) = f([2,16,8],1) - S*p*Je(1)*q(2);
    end
end

%Evaluate shape function
function [S] = shape(q)
    S = 0.5*[q*(q-1); -2*(q+1)*(q-1); (q+1)*q];
end

%Evaluate the gradient of shape function
function [dNdp] = gradshape(q)
    dNdp = 0.5*[2*q-1; -4*q; 2*q+1];
end
