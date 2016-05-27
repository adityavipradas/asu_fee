function [d] = aviprada_hw9(num_nodes, p)
    clc;
    %Flexural modulus
    EI = 2e7;
    %Generate mesh with beam elements
    mesh.x = linspace(0, 8, num_nodes);
    mesh.conn = [1:num_nodes-1; 2:num_nodes];
    %Compute K and f
    quad_pts = [[1,0,-1]*sqrt(3/5); [5,8,5]/9];
    %Initialize
    K = zeros(num_nodes*2);
    f = zeros(num_nodes*2,1);
    d = zeros(num_nodes*2,1);
    %Generate K and f matrices
    for c = mesh.conn
        index = [2*c(1)-1; 2*c(1);2*c(2)-1; 2*c(2)];
        xe = mesh.x(:,c);
        le = xe(2) - xe(1);
        Ke = zeros(4);
        for q = quad_pts
            N = shape(q(1), le);
            B = gradshape(q(1), le);
            Ke = Ke + EI*(B'*B)*(le/2)*q(2);
            x = xe*interpolate(q(1));
            f(index) = f(index) + N'*p(x)*(le/2)*q(2);
        end
        K(index,index) = K(index,index) + Ke;
        %f(index) = f(index) + [1;le/6;1;-le/6]*p(x)*le/2;
    end
    %Boundary conditions
    disp(K);
    d(1) = 0;
    d(2*num_nodes-1) = 0;
    d(num_nodes) = 0;
    d(num_nodes+1) = 0;
    %Solve
    K = K([2:num_nodes-1,num_nodes+2:2*num_nodes-2,2*num_nodes], ...
        [2:num_nodes-1,num_nodes+2:2*num_nodes-2,2*num_nodes]);
    f = f([2:num_nodes-1,num_nodes+2:2*num_nodes-2,2*num_nodes]);
    d([2:num_nodes-1,num_nodes+2:2*num_nodes-2,2*num_nodes]) = K\f;
    %Plot
    plot_solution(mesh, d);
end

function [N] = shape(q, le)
    N(1) = (1/4)*(1-q)^2*(2+q);
    N(2) = (le/8)*(1-q)^2*(1+q);
    N(3) = (1/4)*(1+q)^2*(2-q);
    N(4) = (le/8)*(1+q)^2*(q-1);
end

function [B] = gradshape(q,le)
    B(1) = 6*q/le;
    B(2) = 3*q-1;
    B(3) = -6*q/le;
    B(4) = 3*q+1;
    B = (1/le)*B;
end

function [S] = interpolate(q)
    S = [0.5*(1-q(1)); 0.5*(1+q(1))];
end

function [f] = plot_solution(mesh, d)
    position = [];
    displacement = [];
    for c = mesh.conn
        xe = mesh.x(c);
        h = xe(2) - xe(1);
        sctr(1:2:4) = 2*c-1;
        sctr(2:2:4) = 2*c;
        de = d(sctr);
        for q = linspace(-1,1)
            N2 = 0.5*[1-q(1); 1+q(1)];
            position(end+1) = xe*N2;
            N = shape(q(1), h);
            displacement(end+1) = N*de;
        end
    end
    plot(position, displacement); hold on;
    plot(mesh.x, d(1:2:end), 'o');
end