% Returns a mesh structure containing the following:
%    conn:  2xNe matrix of element nodes.
%    nodes: 2xN  matrix of nodal coordinates.
% The input arguments are:
%        L: the length of the truss.
%        H: the height of the truss at its highest point.
%        n: The number of nodes along a half-span of the truss.
function [mesh] = make_truss(L, H, n)       
    % X-coordinates of base.
    xx = linspace(0, L, 2*n-1);
    % Angular coordinates along top.
    qq = acos(xx/(0.5*L)-1);
    mesh.x = [xx,             xx(2:end-1);
              zeros(1,2*n-1), H*sin(qq(2:end-1))];
    % Connectivity of perimeter elements.
    p = [1:(2*n-1),       1,2*n:(4*n-5);
         2:(2*n-1),4*n-4, (2*n):(4*n-4)];
    % Connectivity of vertical elements.    
    v = [2:2*n-2;   
        (2:2*n-2)+2*n-2];
    % Connectivity of cross elements.
    c = [3:n,          n:2*n-3;
        (3:n)+2*n-3,  (n:2*n-3) + 2*n-1];
    mesh.conn = [p, v, c];
    disp(mesh.conn);
end