%------------------------------------------------------ 
%Baldo & Melian, JULY 2020, Horw, CH


%Hill Numbers based on 

%Gotelli, N. J., Hsieh, T. C., Sander, E. L., & Colwell, R. K. (2014). Rarefaction and Extrapolation with Hill Numbers: A Framework for Sampling and Estimation in Species Diversity Studies. Ecological Monographs, 84(1), 45–67. https://doi.org/10.1890/13-0133.1
%------------------------------------------------------

addpath(genpath(pwd))
pkg load statistics

%DATA ---------------------
fid = fopen('./src/data/Pampa/Pampa_2020.csv');
#fid = fopen('Matrix complete_2020.csv');
in = fscanf(fid,'%c');
fclose(fid);
linesR = regexp(in,'(^,|)([^\n]+)', 'match');
All = char(linesR);
Forest = cell(length(All),1);
wR = unique_no_sort(linesR);
R = char(wR);

#categorize forest type
for k = 1:length(All);
    wHF = regexp(All(k,:),',');
    Forest(k,1) = All(k,1:wHF(1,1)-1);
end
wF = unique_no_sort(Forest);

for n = 2:length(wF); 
TypeOutput = zeros(length(R),4);
Type = zeros(1,1);
ForestType = zeros(1,2);
InsectType = zeros(1,1);
xR = zeros(1,1);
wHR1 = zeros(1,1);

wF(1,n)   
    cR = cellstr(wF(1,n));
    xR = regexpcell(linesR,cR);
    ForestType = cell(length(xR),1:2);
    %InsectType = cell(length(xR),1);
    
    for j = 1:length(xR);
        wHR1 = regexp(All(xR(1,j),:),',');
        ForestType(j,1) = All(xR(1,j),wHR1(1,2)+1:wHR1(1,3)-1);%Plant species
        ForestType(j,2) = All(xR(1,j),wHR1(1,4)+1:wHR1(1,5)-1);%Insect species
    end

    %RarefactionI 
    N = zeros(1,1);
    Type = unique(ForestType);
    N = length(ForestType);
    UNI = length(Type);
    
    %All --- 
    A = ForestType;
    [B,I,J] = unique(ForestType);
    
    for u = 1:length(A);
        New(u,1) = J(u,1);
        New(u,2) = J(u+length(A),1);
    end
    Total = unique(New,"rows");
    red=unifrnd(0,1); green=unifrnd(0,1); blue=unifrnd(0,1);
    
    R=10;
    for rep = 1:R; 
        U = zeros(1,2);
        for r = 1:N;
            S = unidrnd(length(New),r,1);
            Sampled = New(S,1:2);
            L = unique(Sampled,"rows");
            U(r,1) = r;
            U(r,2) = length(L);  
        end
        hold on
        plot(U(:,1),U(:,2),'k','color',[red green blue],'Markersize',36)
    end%rep  
    #pause   
end
set(gca,'FontSize',20)
xlabel('Number of flower visitors','fontsize',30)
ylabel('Number of unique interactions','fontsize',30)
