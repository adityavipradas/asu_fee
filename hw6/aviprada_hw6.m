%TENSOR PRODUCT METHOD TO EVALUATE THE SHAPE FUNCTIONS OF A 9-NODED
%QUADRILATERAL ELEMENT

function [N, dNdp] = aviprada_hw6(p)
%clear screen
    clc;
    E = shape(p(1));
    n = shape(p(2));
    
    N = [E(1)*n(1);
        E(3)*n(1);
        E(3)*n(3);
        E(1)*n(3);
        E(2)*n(1);
        E(3)*n(2);
        E(2)*n(3);
        E(1)*n(2);
        E(2)*n(2)];
    
    gE = gradshape(p(1));
    gn = gradshape(p(2));
    
    dNdp = [gE(1)*n(1), E(1)*gn(1);
        gE(3)*n(1), E(3)*gn(1);
        gE(3)*n(3), E(3)*gn(3);
        gE(1)*n(3), E(1)*gn(3);
        gE(2)*n(1), E(2)*gn(1);
        gE(3)*n(2), E(3)*gn(2);
        gE(2)*n(3), E(2)*gn(3);
        gE(1)*n(2), E(1)*gn(2);
        gE(2)*n(2), E(2)*gn(2);];
    N
    dNdp
end

%shape functions for a quadratic 1-D element
function [S] = shape(v)
     S = [0.5*v*(v-1); 1-v^2; 0.5*v*(1+v)]; 
end

%shape function gradients for a quadratic 1-D element
function [G] = gradshape(v)
    G = [v-0.5; -2*v; v+0.5]; 
end