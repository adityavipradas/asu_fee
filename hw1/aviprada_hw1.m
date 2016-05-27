function [nodes, connectivity] = aviprada_hw1(N, R)
%make sure that N is always an integer
%clear screen
    clc;
   
%determine angle increment
    angle_incr = 2*pi/N;
    angle_step = 0;
    
%store the co-ordinates of nodes in the nodes matrix and the node numbers 
%of elements in the connectivity matrix
    for i = 1:1:N
        nodes(1,i) = R*cos(angle_step);
        nodes(2,i) = R*sin(angle_step);
        connectivity(1,i) = i;
        connectivity(2,i) = i+1;
        if i==N
            connectivity(2,i) = 1;
        end
        angle_step = angle_step + angle_incr;
    end
    
%display nodes and connectivity matrix in the command window    
    nodes
    connectivity
    
%plot the approximated circle and number the nodes   
    figure();
    str = sprintf('Discretizing the circle of radius %0.2f into %d elements',R,N);
    title(str);
    hold on;
    for i = 1:1:N-1
        plot([nodes(1,i), nodes(1,i+1)], [nodes(2,i), nodes(2,i+1)], 'b-o');
        str1 = sprintf('%d',i);
        str2 = sprintf('%d',i+1);
        text(nodes(1,i)+0.01, nodes(2,i)+0.01, str1);
        text(nodes(1,i+1)+0.01, nodes(2,i+1)+0.01, str2);
    end
    plot([nodes(1,N), nodes(1,1)], [nodes(2,N), nodes(2,1)], 'b-o');
end