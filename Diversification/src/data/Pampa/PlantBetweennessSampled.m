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

n = numNodes(adj);
w = size(adj);
%coord
xy = zeros(w(1,1),2);
for z = 1:length(xy);
    xy(z,1) = unidrnd(100);
    xy(z,2) = unidrnd(100);
end
%plot
clf;
gplot (adj, xy, "o-");
set (get (gca, ("children")), "markersize", 12);
title ("gplot() of bit2uni Adjacency matrix");

spaths=inf(n,n);
betw = zeros(1,n);  
adjL = adj2adjL(adj);

for i = iix;#Only species from mbetween
    uB = simpleDijkstra(adj,i);
    for l=1:20;%n
        if i==l; continue; end

        [allPaths, ~] = findAllShortestPaths(adjL,i,l, uB(l), allPaths={},path=[]);
        spaths(i,l) = length(allPaths);
        
        % for all paths in allPaths, parse out the path:
        for p=1:length(allPaths)
            path = strsplit(allPaths{p},'-');
            pathvec = [];
            for x=2:length(path)  % skip the first one
                pathvec = [pathvec str2num(path{x})];
            end
            betw(pathvec(2:length(pathvec)-1)) = betw(pathvec(2:length(pathvec)-1)) + 1/spaths(i,l);
        end
    end  % end of j=1:n
end      % end of i=1:n

#Some outputs
#Total distance vs centrality
#D=sum(dis,2);b=betw';
#plot(D,b,'o')

   
