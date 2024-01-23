% This script is used to analyse the full suite of networks which are
% capable of perfect resilience. This script uses
% 'AllPossibleStructures.mat' which is combined list of all networks as
% either: the full list of terms, or binary values which represent the
% presence/absence of each parameter/term. This script provides two sets of
% plots: a plot of all network parameters (Supplementary figure 2), and set
% of figures grouped by the interspecific and stimulus interactions
% present. The latter figures include the possible arrangements of growth
% and self-regulation terms for each set of interspecific and stimulus 
% interactions which can result in perfect resilience. These figures are
% also plotted alongside calculations of the RPA quality for all three
% species (RPAqualityAll3SpeciesIOTarget.mat, 
% RPAqualityAll3SpeciesITarget.mat, RPAqualityAll3SpeciesOTarget.mat). 
% These latter figures are also loaded alongside initial parameter sets
% for these networks so that specific parameter values can be examined
% alongside the figures.

%% load in data
% load structure data
load Data\AllPossibleStructures % Networks (equations), Binary (1's are 0's)
% load RPA quality for all species
load Data\RPAqualityAll3SpeciesIOTarget; RPAqualityAllIO = RPAqualityAll;
load Data\RPAqualityAll3SpeciesITarget; RPAqualityAllI = RPAqualityAll;
load Data\RPAqualityAll3SpeciesOTarget; RPAqualityAllO = RPAqualityAll;

% load initial parameter sets for all networks
load Data\InitialStableSystemsAuto3SpeciesIOTargetBigger
load Data\RPAquality3SpeciesIOTargetBigger; parasIO = SavedParameterSets(NotEmptySets(RPAquality>=0.1));
load Data\InitialStableSystemsAuto3SpeciesITargetBigger
load Data\RPAquality3SpeciesITargetBigger; parasI = SavedParameterSets(NotEmptySets(RPAquality>=0.1));
load Data\InitialStableSystemsAuto3SpeciesOTargetBigger
load Data\RPAquality3SpeciesOTargetBigger; parasO = SavedParameterSets(NotEmptySets(RPAquality>=0.1));

%% Data Manipulation
% network data
InterInd = [5,14,3,4,7,9,11,12];
[OrderedMat,OrderInd] = sortrows(double(Binary),InterInd);
[UniqueSets,~,UniqueInd] = unique(OrderedMat(:,InterInd),'rows');

% RPA quality data
RPAQualities = [RPAqualityAllO'; RPAqualityAllI'; RPAqualityAllIO'];
RPAQualities = RPAQualities(OrderInd,:);
% parameter data
AllParas = [parasO; parasI; parasIO];
AllParas = AllParas(OrderInd);

%% Plot all networks
labelstring = {'r_I','a_{II}','a_{IM}','a_{IO}','d_I','r_M','a_{MI}','a_{MM}','a_{MO}','r_O','a_{OI}','a_{OM}','a_{OO}','d_O'};
figure,imagesc([~OrderedMat(1:113,:), ~OrderedMat(114:226,:), ~OrderedMat(227:339,:)])
hold on
for i = 4:4:113
    plot([0.5 42.5],[i+0.5, i+0.5],'color',[0.5 0.5 0.5])
end
for i = 0:41
    plot([i+0.5 i+0.5],[0.5 113.5],'color',[0.8 0.8 0.8])
end
plot([14.5 14.5],[0.5 113.5],'color',[0.4 0.4 0.4],'LineWidth',3)
plot([28.5 28.5],[0.5 113.5],'color',[0.4 0.4 0.4],'LineWidth',3)
set(gca,'XTick',1:42,'XTickLabel',labelstring)
set(gca,'YTick',10:10:113,'YTickLabel',10:10:113)
colormap('bone')

%% Plot split network figures
for i = 1:length(UniqueSets)
figure(i), subplot(1,4,1:3)
PlotMat = OrderedMat(UniqueInd==i,:);
imagesc(PlotMat)
set(gca,'XTick',1:14,'XTickLabel',{'r1','a11','a12','a13','d1','r2','a21','a22','a23','r3','a31','a32','a33','d3'})
set(gca,'YTick',1:sum(UniqueInd==i),'YTickLabel',1:sum(UniqueInd==i))
set(gcf,'Position',[794,217,555,420])
colormap('bone')
subplot(1,4,4)
PlotQuality = RPAQualities(UniqueInd==i,:);
imagesc(PlotQuality<0.1)
set(gca,'XTick',1:3,'XTickLabel',{'I','M','O'})
end

