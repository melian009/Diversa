#----------------------------------------------------------------------------------
#Multilayer rarefaction random
#Betweenness centrality from spatial-interaction-biodiversity data
# Melian@CH, Multilayer networks OCT NOV 2021
#----------------------------------------------------------------------

#------------------------------------------------------ 
#Betweenness centrality based on 
#https://en.wikipedia.org/wiki/Centrality
#------------------------------------------------------

addpath(genpath(pwd))
pkg load statistics

# Make a random spatial matrix
s1 = 10;s2 = 10;adjS = rand(s1,s2)>.75;

# betwenneess sites
S = numNodes(adjS);spaths=inf(s1,s2);
betwS = zeros(1,s1);adjL = adj2adjL(adjS);

for i=1:s1;
    uB = simpleDijkstra(adjS,i);
    for j=1:s2;
        if i==j; continue; end
        [allPaths, ~] = findAllShortestPaths(adjL,i,j, uB(j), allPaths={},path=[]);
        spaths(i,j) = length(allPaths);
        # for all paths in allPaths, parse out the path:
        for pat=1:length(allPaths)
            path = strsplit(allPaths{pat},'-');
            pathvec = [];
            for x=2:length(path) # skip the first one
                pathvec = [pathvec str2num(path{x})];
            end
            betwS(pathvec(2:length(pathvec)-1)) = betwS(pathvec(2:length(pathvec)-1)) + 1/spaths(i,j);
        end
    end  # end of j=1:s1
end      # end of i=1:s2

#P/A plant 
pl = 10;PA = rand(pl,s1)>.75
BC = zeros(pl,s1);

for ii = 1:pl;#loop plant i
ii
    red=unifrnd(0,1); green=unifrnd(0,1); blue=unifrnd(0,1);#def color each species i
    for jj = 1:s1;
        if PA(ii,jj) == 1;        
           a = unidrnd(15);A = rand(pl,a)>.75;# Make a random MxN adjacency bipartite matrix
           adj = [zeros(pl,pl), A;A', zeros(a,a)];# Expand out to symmetric (M+N)x(M+N) matrix
           PlantBetweennessRandom;
           Bc(ii,jj) =  betwS(1,jj) * betw(1,ii);
        else
           PA(ii,jj) = 0;
        end
    end#site
Bci = sort(Bc(ii,:),2,'descend');
rank = 1:length(Bci);
hold on
plot(rank,Bci,'k','color',[red green blue],'Markersize',12)
set(gca,'FontSize',12)
xlabel('Rank species centrality','fontsize',12)
ylabel('Betweenness','fontsize',12) 
end




