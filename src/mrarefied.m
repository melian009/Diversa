#----------------------------------------------------------------------------------
#Plot multilayer rarefied curves from biodiversity data
#CJ Melian, Baldo F., Multilayer networks MAR2020-OCT 2021
#----------------------------------------------------------------------

%------------------------------------------------------ 
%Hill Numbers based on 

%Gotelli, N. J., Hsieh, T. C., Sander, E. L., & Colwell, R. K. (2014). Rarefaction and Extrapolation with Hill Numbers: A Framework for Sampling and Estimation in Species Diversity. Ecological Monographs, 84(1), 45–67. https://doi.org/10.1890/13-0133.1
%------------------------------------------------------

pkg load statistics
pkg load dataframe
pkg load symbolic

%DATA ---------------------
fid = fopen('../Data/Pampa/Pampa_2020.csv');
in = fscanf(fid,'%c');
fclose(fid);
linesR = regexp(in,'(^,|)([^\n]+)', 'match');
All = char(linesR);
Forest = cell(length(All),1);
wR = unique_no_sort(linesR);

%1. Print ALL raw data ----------------------------
D = char(wR);
%----------------------------

#categorize forest type (sites)
for k = 1:length(All);
    wHF = regexp(All(k,:),',');
    Forest(k,1) = All(k,1:wHF(1,1)-1);
end


%2. List sampled sites ----------------------------------------------
wF = unique_no_sort(Forest);%pause
%----------------------------------------------

for n = 2:4;%length(wF); 

    for m = n+1:5;%:length(wF);%sampling beta

%1st sierra
TypeOutput1 = zeros(length(D),4);
Type1 = zeros(1,1);
ForestType1 = zeros(1,1);
xR1 = zeros(1,1);
wHR1 = zeros(1,1);

%2nd sierra
TypeOutput2 = zeros(length(D),4);
Type2 = zeros(1,1);
ForestType2 = zeros(1,1);
xR2 = zeros(1,1);
wHR2 = zeros(1,1);


wF(1,n);%first sierra
    cR1 = cellstr(wF(1,n));
    xR1 = regexpcell(linesR,cR1);
    ForestType1 = cell(length(xR1),1);
    
    for j = 1:length(xR1);
        wHR1 = regexp(All(xR1(1,j),:),',');
        ForestType1(j,1) = All(xR1(1,j),wHR1(1,2)+1:wHR1(1,3)-1);
    end
    
    
wF(1,m);%second sierra
    cR2 = cellstr(wF(1,m));
    xR2 = regexpcell(linesR,cR2);
    ForestType2 = cell(length(xR2),1);
    
    for j = 1:length(xR2);
        wHR2 = regexp(All(xR2(1,j),:),',');
        ForestType2(j,1) = All(xR2(1,j),wHR2(1,2)+1:wHR2(1,3)-1);
    end
      
   
    
%.3 Size each site -------------------------------------------
wF(1,n);
N1 = length(ForestType1);
N2 = length(ForestType2);
%wFcom(1,n) = N1;
%----------------------------------------------    
    
%m-Rarefaction: 
%A. Pick up more-sampled site
%B. Obtain n-subsamplings of size of less-sampled site
%C. Compare alpha and beta curves along subsamplings of different sizes
%wFcom(1,n) = N;
%MIN = min(wFcom(1,2:12));
MIN = 348;%Check

%4. Compare curves
Type1 = unique(ForestType1);
UNI1 = length(Type1);

Type2 = unique(ForestType2);
UNI2 = length(Type2);

%clf;
%labels = {};
%colororder = get (gca, "colororder");  
  
  red1=unifrnd(0,1); green1=unifrnd(0,1); blue1=unifrnd(0,1);
  red2=unifrnd(0,1); green2=unifrnd(0,1); blue2=unifrnd(0,1);
  red3=unifrnd(0,1); green3=unifrnd(0,1); blue3=unifrnd(0,1);
R=10;
for rep = 1:R; 
      U1 = zeros(1,2);
      U2 = zeros(1,2);
      
      for r = 1:MIN;
          S1 = unidrnd(length(xR1),r,1)';%R replicates
          L1 = unique(ForestType1(S1(1,1:r)));
          U1(r,1) = r;
          U1(r,2) = length(L1); 

%test accumulating sampled species
%S1
%length(L1)
%pause
          S2 = unidrnd(length(xR2),r,1)';%R replicates
          L2 = unique(ForestType2(S2(1,1:r)));
          U2(r,1) = r;
          U2(r,2) = length(L2); 
          
          %pause
          
          %L1
          %L2
          %I
          %pause
          I = intersect(L1,L2,'stable');
          S1 = length(L1); 
          S2 = length(L2);
          %beta : number of shared U(A,B) / (A + B - U(A,B))
          U3(r,1) = r;
          U3(r,2) = length(I) / ((S1 + S2) - length(I));
         
      end
      
    
     
      hold on
      %plot alpha two sierras
      subplot(2,2,1)
      plot(U1(:,1),U1(:,2),'k','color',[red1 green1 blue1],'Markersize',36)
      hold on
      %red=unifrnd(0,1); green=unifrnd(0,1); blue=unifrnd(0,1);
      plot(U2(:,1),U2(:,2),'k','color',[red2 green2 blue2],'Markersize',36)
      
      set(gca,'FontSize',12)
      xlabel('Number Individuals Sampled ','fontsize',12)
      ylabel('\alpha','fontsize',12)
      
      %plot beta two sierras 
      subplot(2,2,2)
      hold on
      %red=unifrnd(0,1); green=unifrnd(0,1); blue=unifrnd(0,1);
      plot(U3(:,1),U3(:,2),'k','color',[red3 green3 blue3],'Markersize',36)
      
      set(gca,'FontSize',12)
      xlabel('Number Individuals Sampled ','fontsize',12)
      ylabel('\beta','fontsize',12)
      %pause
      
end%rep 
 %hold on;
 %set (h, "color", colororder(n,:));
 %labels = {labels{:}, ["Signal ", num2str(n)]};

end%n
end%m

%set(gca,'FontSize',20)
%xlabel('Number Individuals Sampled ','fontsize',22)
%ylabel('Species Number','fontsize',22)
%legend (labels, "location", "southoutside");

%h = legend ({"Difundito"}, "La Paja", "Cinco cerros", "Barrosa", "Amarante", "Difuntos", "La Chata", "Piedra Alta", "Vigilancia", "Volcan", "La Brava", "El Morro");
%legend (h, "location", "northeastoutside");
%set (h, "fontsize", 20);





