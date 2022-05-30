#----------------------------------------------------------------------------------
#Multilayer rarefaction 
#Betweenness centrality from spatial-interaction-biodiversity data
# Melian@CH, Multilayer networks OCT NOV 2021
#----------------------------------------------------------------------

#------------------------------------------------------ 
#Betweenness centrality based on 
#https://en.wikipedia.org/wiki/Centrality
#------------------------------------------------------

addpath(genpath(pwd))
pkg load statistics
#pkg load dataframe
#pkg load symbolic
#pkg load matgeom

# 1. DATA --------------------------------------
fid = fopen('./data/Pampa/Pampa_2020.csv');
in = fscanf(fid,'%c');
fclose(fid);

linesR = regexp(in,'(^,|)([^\n]+)', 'match');
All = char(linesR);
Species = cell(length(All),1);
Sierras = cell(length(All),1);
wR = unique_no_sort(linesR);
D = char(wR);
#----------------------------

# 2. Forest type (species) ---
for ksp = 1:length(All);
    wHF = regexp(All(ksp,:),',');
    Species(ksp,1) = All(ksp,wHF(1,2)+1:wHF(1,3)-1);#Plants --- How about insects?
end
wSp = unique_no_sort(Species);#pause
#----------------------------
 
# 3. Forest type (sierras) ---
for ksi = 1:length(All);
    wHF = regexp(All(ksi,:),',');
    Sierras(ksi,1) = All(ksi,1:wHF(1,1)-1);
end
wSi = unique_no_sort(Sierras);
#---------------------------

# 4. Betweenness centrality each SAMPLED sierra ------------------------------
#Get betweenn... each sierra from SierraBetweennessSampled.m
#Build BetweenSierra.csv
fid = fopen('./data/Pampa/BetweenSierra.csv');#20, 30 Kms thresholds
in = fscanf(fid,'%c');
fclose(fid);
linesS = regexp(in,'(^,|)([^\n]+)', 'match');
AllB = char(linesS);
AllBSize = size(AllB);
SierrasB = cell(AllBSize(1,1),1);
SierrasC = cell(AllBSize(1,1),1);
a=size(AllB);


# 5. Betweenness centrality each sierra
for ksii = 1:a(1,1);
    wHS = regexp(AllB(ksii,:),',');
    SierrasC(ksii,2) = AllB(ksii,wHS(1,1)+1:wHS(1,2)-1);#<=10
    SierrasC(ksii,1) = AllB(ksii,1:wHS(1,1)-1);#<=20
end
%--------------------------------------------

# 5. TODO: Bc each sierra ALL sampled-non sampled sierras

# 6. Betwenn... species i for all sierras j ---------------------------------
Bc=zeros(length(wSp)-1,length(wSi)-1);#Multilayer centrality
for i = 2:length(wSp);#...for species

wSp(1,i);
i;


cR = cellstr(wSp(1,i));
xR = regexpcell(linesR,cR);#presence species i all

AllSierras = char(xR);
AllxR = size(AllSierras);
LocateSierras = cell(AllxR(1,1),1);

for ksi = 1:length(xR);
    wHF = regexp(All(xR(1,ksi),:),',');
    LocateSierras(ksi,1) = All(xR(1,ksi),1:wHF(1,1)-1);
end
SSpi = unique_no_sort(LocateSierras);
     for j = 1:length(SSpi);#presence sierras species i
         #Bipartite sierra j
         Sierra = zeros(1,2);xR1 = zeros(1,1);wHR1 = zeros(1,1);FT1 = zeros(1,1);
         cR1 = cellstr(SSpi(1,j));xR1 = regexpcell(linesR,cR1);
         Sierra = cell(length(xR1),1:2);
                
        for k = 1:length(xR1);
            wHR1 = regexp(All(xR1(1,k),:),',');
            Sierra(k,1) = All(xR1(1,k),wHR1(1,2)+1:wHR1(1,3)-1);%Plant species
            Sierra(k,2) = All(xR1(1,k),wHR1(1,4)+1:wHR1(1,5)-1);%Insect species
        end
        
        #Remove No Visit
        removeIndex = strcmp(Sierra(:,2),'No visit');
        Sierra(removeIndex,:) = [];

        
        listSPi = regexpcell(Sierra(:,1),wSp(1,i))
        if ~isempty(listSPi);
        
        #Unique   
        FT = strcat(Sierra(:,1),{' '},Sierra(:,2));
        A = unique(FT);
        plants=unique(Sierra(:,1));
        insects=unique(Sierra(:,2));
        Adj = zeros(length(plants),length(insects));
        
        #presence species i bipartite sierra j
        xRPre = regexpcell(A,cR);#presence species i all
        for k1 = 1:length(plants);
            listPla = regexpcell(A,plants(k1,1));
            for k2 = 1:length(listPla);#locate insects int add 1 to Adj k1,k2    
                listPoll = regexpcell(A(listPla,1),insects(k2,1));
                if ~isempty(listPoll);
                   #Check
                   Adj(k1,k2) = 1;
                   #A(listPla,1);insects(k2,1);pause
                else
                   Adj(k1,k2) = 0;
                end
            end
        end
        
        # Convert bipartite to unipartite bip2uni
        #length(plants)#length(insects)
        adj = [zeros(length(plants),length(plants)), Adj;Adj', zeros(length(insects),length(insects))];
        q = size(adj)
        
        #if q(1,1) <= 50;
        #Call between from SierraBetweenessSampled for species j
        iix = i;
        PlantBetweennessSampled;
         
         CEN = zeros(1,1);
         SierraSpj = cellstr(SSpi(1,j));
         CEN = regexpcell(SierrasC,SierraSpj); 
         
        #Test        
        #SSpi(1,j) LocateSierras SierrasC
        wSp(1,i) 
        
        #Find SSpi in SierrasC 
        Bc(i-1,j) = str2double(SierrasC(CEN,2)) * betw(1,i-1);
        #end#if size


        else
           Bc(i-1,j) = 0;
        end

     end#sierras j
end#species i


