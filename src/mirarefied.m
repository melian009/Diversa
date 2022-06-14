#-----------------------------------------------------------------------
# Multilayer networks - Melian@CH - Multilayer rarefaction OCT NOV 2021 
# ----------------------------------------------------------------------

pkg load statistics
addpath(genpath(pwd))

# 0. DATA ------------------------
fid = fopen('./data/Pampa/Pampa_2020.csv');#Open raw data
in = fscanf(fid,'%c');fclose(fid);
linesR = regexp(in,'(^,|)([^\n]+)', 'match');
All = char(linesR);Forest = cell(length(All),1);
wR = unique_no_sort(linesR);R = char(wR);

#Distance matrix
DistanceMatrixPampaSampled;
#---------------------------------

# 1. List sampled sites
for k = 1:length(All);wHF = regexp(All(k,:),',');
    Forest(k,1) = All(k,1:wHF(1,1)-1);
end
wF = unique_no_sort(Forest);#pause


count = 0;#for intersection

# 2. Loop shared interactions sites(n,m)
for n = 2:length(wF)-1;
red1=unifrnd(0,1); green1=unifrnd(0,1); blue1=unifrnd(0,1);
wF(1,n)
n
        # Site n
        ForestType1 = zeros(1,2);xR1 = zeros(1,1);wHR1 = zeros(1,1);FT1 = zeros(1,1);
        cR1 = cellstr(wF(1,n));xR1 = regexpcell(linesR,cR1);
        ForestType1 = cell(length(xR1),1:2);
                
        for j = 1:length(xR1);
            wHR1 = regexp(All(xR1(1,j),:),',');
            ForestType1(j,1) = All(xR1(1,j),wHR1(1,2)+1:wHR1(1,3)-1);%Plant species
            ForestType1(j,2) = All(xR1(1,j),wHR1(1,4)+1:wHR1(1,5)-1);%Insect species
        end
        
    for m = n+1:length(wF)-1;
    m
        red2=unifrnd(0,1); green2=unifrnd(0,1); blue2=unifrnd(0,1);
        red3=unifrnd(0,1); green3=unifrnd(0,1); blue3=unifrnd(0,1);
        wF(1,m) 
       
        
        # Site m 
        ForestType2 = zeros(1,2);xR2 = zeros(1,1);wHR2 = zeros(1,1);FT2 = zeros(1,1);
        cR2 = cellstr(wF(1,m));xR2 = regexpcell(linesR,cR2);
        ForestType2 = cell(length(xR2),1:2);
                
        for j = 1:length(xR2);
            wHR2 = regexp(All(xR2(1,j),:),',');
            ForestType2(j,1) = All(xR2(1,j),wHR2(1,2)+1:wHR2(1,3)-1);#Plant species
            ForestType2(j,2) = All(xR2(1,j),wHR2(1,4)+1:wHR2(1,5)-1);#Insect species
        end
        
   #Remove No Visit
   removeIndex1 = strcmp(ForestType1(:,2),'No visit');
   ForestType1(removeIndex1,:) = [];
   
   removeIndex2 = strcmp(ForestType2(:,2),'No visit');
   ForestType2(removeIndex2,:) = [];
   #pause
    
    #Concatenate interactions to unique
    FT1 = strcat(ForestType1(:,1),{' '},ForestType1(:,2));
    FT2 = strcat(ForestType2(:,1),{' '},ForestType2(:,2));
    #W = intersect(FT1,FT2);#shared unique
        
    # 3. Pairs alpha-beta interaction curves
    M=zeros(1,1);W12=zeros(1,2);U1=zeros(1,2);U2=zeros(1,2);UN1 = zeros(1,1);UN2 = zeros(1,1);w=zeros(1,1);
    S1 = zeros(1,1);S2 = zeros(1,1);L1 = cell(1,1);L2 = cell(1,1);In = zeros(1,2);

    M = min(length(FT1),length(FT2));#account equal sampling size
    for r = 1:M;
    
        S1 = unidrnd(length(FT1),1,1)';L1(r,1) = FT1(S1,1);UN1 = unique(L1(:,1));
        U1(r,1) = r;U1(r,2) = length(UN1); 
   
        S2 = unidrnd(length(FT2),1,1)';L2(r,1) = FT2(S2,1);UN2 = unique(L2(:,1));
        U2(r,1) = r;U2(r,2) = length(UN2); 
        
        w = intersect(UN1,UN2);#shared unique
        W12(r,1) = length(w);
        W12(r,2) = dis(n,m);#Distance ij
        
        
    end
    
    # 4. Intersections interlayer
    count = count + 1;q = zeros(1,1);
    In(count,1) = dis(n,m); 
    intersection=find(U1(:,2)==U2(:,2));
    q = find(intersection) >= 10;
    In(count,2) = nnz(q);
    #pause
    
    #. 5. PLOT
    hold on
    subplot(2,2,1)
    plot(U1(:,1),U1(:,2),'k','color',[red1 green1 blue1],'Markersize',14)
    hold on
    plot(U2(:,1),U2(:,2),'k','color',[red2 green2 blue2],'Markersize',14)
    set(gca,'FontSize',12)
    xlabel('Number of flower visitors','fontsize',12)
    ylabel('Number of unique interactions','fontsize',12)
    
    hold on
    subplot(2,2,2)
    plot(U2(:,1),W12(:,1),'k','color',[red3 green3 blue3],'Markersize',14)
    set(gca,'FontSize',12)
    xlabel('Number of flower visitors','fontsize',12)
    ylabel('Number of shared interactions','fontsize',12)
    
    hold on   
    subplot(2,2,3)
    plot(W12(:,2),W12(:,1),'ko','color',[red3 green3 blue3],'Markersize',14)
    set(gca,'FontSize',12)
    xlabel('Pairwise distance','fontsize',12)
    ylabel('Number of shared interactions','fontsize',12)   
    
    hold on
    subplot(2,2,4)
    plot(In(:,1),In(:,2),'ko','color',[red3 green3 blue3],'Markersize',14)
    set(gca,'FontSize',12)
    xlabel('Pairwise distance','fontsize',12)
    ylabel('Number of interlayer intersections','fontsize',12) 
 
    #hold off
   end#m Sites       
end#n Sites

#print -color muinrarefaction_v2.eps

