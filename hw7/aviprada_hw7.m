function [K] = aviprada_hw7(nh)
%clear screen
    clc;
%create mesh    
    mesh = hole_mesh(nh, 2.5, 0.4);
    %disp(mesh.conn);
    %disp(mesh.x);
%parent coordinates
    E = [-1 1 1 -1];
    n = [-1 -1 1 1];
%gauss-quadrature    
    q = [[-1, 1]/sqrt(3); [-1, 1]/sqrt(3); ...
        [1, 1]; [1, 1]];
%initialize    
    K = zeros(length(mesh.x));
%form conductivity matrix    
    for c = mesh.conn
        xe = mesh.x(:,c);
        %disp(xe);
        Ke = zeros(4);
        Kcond = [205 0;0 205];
        for i = 1:1:2
            for j = 1:1:2
                N = shape([q(1,i),q(2,j)], E, n);
                dNdp = gradshape([q(1,i),q(2,j)], E, n);
                Je = xe * dNdp;
                dNdx = dNdp / Je;
                %disp(dNdx);
                Ke = Ke + dNdx*Kcond*dNdx'*det(Je)*q(3,i)*q(4,j);
            end
        end
        K(c,c) = K(c,c) + Ke;
    end
    spy(K);
end

%four-node quadrilateral element
function [N] = shape(q, E, n)
    for i = 1:1:4
        N(i,1) = (1 + E(i)*q(1))*(1 + n(i)*q(2))/4;
    end
end

function [dNdp] = gradshape(q, E, n)
    for i = 1:1:4
        dNdp(i,1) = E(i)*(1 + n(i)*q(2))/4;
        dNdp(i,2) = (1 + E(i)*q(1))*n(i)/4;
    end
end
