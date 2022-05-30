## Difusion on single graphs
graphics_toolkit();
available_graphics_toolkits();
loaded_graphics_toolkits();
graphics_toolkit("gnuplot")

addpath(genpath(pwd))
pkg load statistics

N = 20; #Number of nodes
A = zeros(N, N); 
Adj = zeros(N * N, N * N); % The adjacency matrix
size(Adj);


% Use 8 neighbors, and fill in the adjacency matrix (as an example)
dx = [- 1, 0, 1, - 1, 1, - 1, 0, 1];
dy = [- 1, - 1, - 1, 0, 0, 1, 1, 1];
for x = 1:N
    for y = 1:N
        index = (x - 1) * N + y;
        for ne = 1:length(dx)
            newx = x + dx(ne);
            newy = y + dy(ne);
            if newx > 0 && newx <= N && newy > 0 && newy <= N
                index2 = (newx - 1) * N + newy;
                Adj(index, index2) = 1;
            end
        end
    end
end

## COMPUTES THE SOLUTION TO THE DIFFERENTIAL EQUATION
Deg = diag(sum(Adj, 2)); % Compute the degree matrix
L = Deg - Adj; % Compute the laplacian matrix in terms of the degree and adjacency matrices
[V, D] = eig(L); % Compute the eigenvalues/vectors of the laplacian matrix
D = diag(D);

## Initial condition (place a few large positive values around and make everything else zero)
C0 = zeros(N, N);
C0(2:5, 2:5) = 5;
C0(10:15, 10:15) = 10;
C0(2:5, 8:13) = 7;
C0 = C0(:);

C0V = V'*C0; % Transform the initial condition into the coordinate system of the eigenvectors
count = 0;
red1=unifrnd(0,1); green1=unifrnd(0,1); blue1=unifrnd(0,1);
for t = 0:0.1:1;
    count = count + 1;
    % Loop through times and decay each initial component
    Phi1 = t;
    Phi2 = C0V .* exp(- D * t); % Exponential decay for each component
    hold on
    plot(Phi1,Phi2,'ko','color',[red1 green1 blue1],'Markersize',14)
    #pause
end

