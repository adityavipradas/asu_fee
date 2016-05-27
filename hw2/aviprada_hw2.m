function [conn, K] = aviprada_hw2(h, t, D, E1, E2)
    %h, t, D are in meters and E1, E2 are in Pa
    %clear screen
    clc;
    
    %create element cross-sectional area(m2) and length(m) array
    A([1 5]) = h*t;
    L([1 5]) = 4e-2;
    A(2:4) = pi*D^2/4;
    L(2:4) = 3e-2;
    
    %create element Young's modulus array
    E([1 2 4 5]) = E1;
    E(3) = E2;
    
    %create the connectivity matrix
    conn = [1 2 2 2 3;2 3 3 3 4];
    
    %initialize the global stiffness matrix
    K = zeros(4,4);
    
    for i=1:1:5
        %define the element stiffness matrix
        Kele{i} = (E(i)*A(i)/L(i))*[1 -1;-1 1];
        c = conn(:,i);
        %assemble in the global stiffness matrix by direct assembly method
        K(c,c) = K(c,c) + Kele{i};
    end
    conn
    K
end
