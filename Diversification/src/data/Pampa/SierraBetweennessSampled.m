# This function returns the betweenness of all nodes in an adjacency (transformed from a bipartite) 
# Accounts for site centrality
# Betweenness centrality measure: number of shortest paths running through a vertex.

#INPUTS: Raw Pampa data - spatiotemporal pollination in n sierras : DistanceMatrixPampaSampledDouble.csv
#INPUTS: Adjacency matrix (nxn)

#OUTPUTS: betweeness vector for each node (1xn)

#Routines used: numNodes.m, adj2adjL.m, simpleDijkstra.m, findAllShortestPaths.m
#-------------------------------------------------------------------------------------

# https://en.wikipedia.org/wiki/Betweenness_centrality
# Example with a bipartite graph 

pkg load statistics

#function betw = SierraBetweennessSampled(adj)

#dis = dlmread('../Data/Pampa/DistanceMatrixPampaSampledDouble.csv');
DistanceMatrixPampaSampled;
adj = zeros(length(dis),length(dis));
%Convert adj matrix threshold <=8km
for i = 1:length(dis);
    for j = i+1:length(adj);
        if dis(i,j) <= 10;%arbitrary distance threshold
           adj(i,j) = 1;
           adj(j,i) = 1;
        else
           adj(i,j) = 0;
        end
    end
end

n = numNodes(adj);
w = size(adj);
%coord
xy = zeros(w(1,1),2);
for z = 1:length(xy);
    xy(z,1) = unidrnd(10);
    xy(z,2) = unidrnd(10);
end
%plot
clf;
gplot (adj, xy, "o-");
set (get (gca, ("children")), "markersize", 12);
title ("gplot() of bit2uni Adjacency matrix");

spaths=inf(n,n);
betw = zeros(1,n);  
adjL = adj2adjL(adj);

for i=1:n
    uB = simpleDijkstra(adj,i);
    for j=1:n
    %i
    %j
        if i==j; continue; end

        [allPaths, ~] = findAllShortestPaths(adjL,i,j, uB(j), allPaths={},path=[]);
        spaths(i,j) = length(allPaths);
        
        % for all paths in allPaths, parse out the path:
        for p=1:length(allPaths)
            path = strsplit(allPaths{p},'-');
            pathvec = [];
            for x=2:length(path)  % skip the first one
                pathvec = [pathvec str2num(path{x})];
            end
            betw(pathvec(2:length(pathvec)-1)) = betw(pathvec(2:length(pathvec)-1)) + 1/spaths(i,j);
        end
    end  % end of j=1:n
end      % end of i=1:n

#Some outputs
#Total distance vs centrality
D=sum(dis,2);b=betw';
plot(D,b,'o')

   
