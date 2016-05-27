function [K] = aviprada_hw3(L, H, nb)
%clear screen
    clc;
%N ---> number of nodes
%Ne --> number of elements
%generate mesh connectivity and nodal coordinates using make_truss.m
    mesh = make_truss(L, H, nb);
%initialize the global stiffness matrix K with dim 2*N    
    K = zeros(2*length(mesh.x));
%choose one element at a time    
    for c = mesh.conn
%extract the x, y coordinates of each node of the element        
        xe = mesh.x(:,c);
%find x, y lengths of the element        
        dx = xe(:,2) - xe(:,1);
%determine the rotation matrix using direction cosines        
        Re = [dx' 0 0; 0 0 dx']/norm(dx);
%find the element stiffness matrix in local coordinates
        K_elem = [1 -1;-1 1]*(pi*0.01^2)*(70e9)/norm(dx);
%find the element stiffness matrix in global coordinates        
        Ke = Re'*K_elem*Re;
%determine the indices of x displacements        
        pos(1:2:4) = 2*c - 1;
%determine the indices of y displacements        
        pos(2:2:4) = 2*c;
%add the element stiffness matrix in the global stiffness matrix        
        K(pos,pos) = K(pos,pos) + Ke;
    end
%convert K to a sparse matrix
disp(eig(K));
K = sparse(K);
end