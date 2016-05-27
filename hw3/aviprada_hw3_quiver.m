function [] = aviprada_hw3_quiver(L, H, nb)
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
%determine eigenvalues and eigenvectors of K
    [ev d] = eig(K);
    disp(sparse(d));
%find the x,y,u and v vectors from eigenvectors for quiver plots    
    for i = 1:1:2*length(mesh.x)
        k = 1;
        figure();
        for j = mesh.conn
            x(k) = ev(2*j(1)-1,i);
            y(k) = ev(2*j(1),i);
            u(k) = ev(2*j(2)-1,i) - ev(2*j(1)-1,i);
            v(k) = ev(2*j(2),i) - ev(2*j(1),i);
            k = k + 1;
%plot the coordinates            
            plot([ev(2*j(1)-1,i), ev(2*j(2)-1,i)], [ev(2*j(1),i), ...
                ev(2*j(2),i)]);
            hold on;
        end
        figure();
        quiver(x,y,u,v,0);
    end
%plot the original truss    
    figure();
    for j = 1:1:length(mesh.conn)
        plot([mesh.x(1,mesh.conn(1,j)), mesh.x(1,mesh.conn(2,j))], ...
            [mesh.x(2,mesh.conn(1,j)), mesh.x(2,mesh.conn(2,j))]);
        hold on;
    end
%convert K to a sparse matrix
    K = sparse(K);
end