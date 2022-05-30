% Make a random MxN adjacency bipartite matrix
m = 5
n = 4
a = rand(m,n)>.2;

% Expand out to symmetric (M+N)x(M+N) matrix
adj = [zeros(m,m), a;a', zeros(n,n)];


%g = graph(biga);
% Plot
%h = plot(g);
% Make it pretty
%h.XData(1:m) = 1;
%h.XData((m+1):end) = 2;
%h.YData(1:m) = linspace(0,1,m);
%h.YData((m+1):end) = linspace(0,1,n);
