% this script plots a generic version of all RPA capable networks along
% with projection monomials. Ideally would like to be ableto also perform a
% Ma style analysis to check for transient dynamics
SaveNetworkName = '3SpeciesITarget';
% read in data
load(['InitialStableSystemsAuto' SaveNetworkName]);
load(['Monomials' SaveNetworkName]);
load(['AcceptedEqnsAuto' SaveNetworkName]);
load(['GroebnerBases' SaveNetworkName],'EquationStarts','FullEquation','InputTargets','NumSpp','Species');

SortNetworks = 0; % sort networks based on monomial equation

% Calculate the quality of RPA in Networks
QualityCalc = 0;
if QualityCalc == 1
    RPAquality = PerturbedTimeSeriesV2(NotEmptySets,SavedParameterSets,FullEquation,EquationStarts,InputTargets,NumSpp,Species);
    save(['RPAquality' SaveNetworkName],'RPAquality');
else
    load(['RPAquality' SaveNetworkName],'RPAquality');
end

% NotEmptySets = NotEmptySets(RPAquality>0.3);

% sort based groebner monomials
GBeqns = sum(KnownMonomials.*Monomials(index(NotEmptySets),:),2);
if SortNetworks == 1
    [SortedGBeqns,SortedIndex] = sort(GBeqns);
else
    SortedGBeqns = GBeqns; SortedIndex = 1:length(GBeqns);
end

% start figure
figure, tiledlayout('flow','TileSpacing','compact')

% loop over parameter sets
for i = 1:length(NotEmptySets)
nexttile, hold on

% plot a point with color that determines the quality of RPA (white is very bad)
scatter(1,1,100,([1,1,1]-[0,0.5,1].*RPAquality(SortedIndex(i))),'filled');

% plot letters
text(0,0,'M','HorizontalAlignment',"center",'VerticalAlignment',"middle")
text(0,1,'I','HorizontalAlignment',"center",'VerticalAlignment',"middle")
text(1,1,'S','HorizontalAlignment',"center",'VerticalAlignment',"middle")
text(1,0,'O','HorizontalAlignment',"center",'VerticalAlignment',"middle")

% plot the interactions
if EquationParas(NotEmptySets(SortedIndex(i)),1) ~= 0
    plot([-0.05 -0.05],[1.3 1.1],'g','linewidth',2)
end
if EquationParas(NotEmptySets(SortedIndex(i)),2) ~= 0
    plot([0.05 0.05],[1.3 1.1],'r','linewidth',2)
end
if EquationParas(NotEmptySets(SortedIndex(i)),3) ~= 0
    plot([-0.05 -0.05],[0.9 0.1],'b','linewidth',2)
end
if EquationParas(NotEmptySets(SortedIndex(i)),4) ~= 0
    plot([0.1 0.9],[0.9 0.1],'b','linewidth',2)
end
if EquationParas(NotEmptySets(SortedIndex(i)),EquationStarts(2)) ~= 0
    plot([-0.05 -0.05],[-0.3 -0.1],'g','linewidth',2)
end
if EquationParas(NotEmptySets(SortedIndex(i)),EquationStarts(2)+1) ~= 0
    plot([0.05 0.05],[0.1 0.9],'linewidth',2,'color',"#7E2F8E")
end
if EquationParas(NotEmptySets(SortedIndex(i)),EquationStarts(2)+2) ~= 0
    plot([0.05 0.05],[-0.3 -0.1],'r','linewidth',2)
end
if EquationParas(NotEmptySets(SortedIndex(i)),EquationStarts(2)+3) ~= 0
    plot([0.1 0.9],[0.05 0.05],'linewidth',2,'color',"#7E2F8E")
end
if EquationParas(NotEmptySets(SortedIndex(i)),EquationStarts(3)) ~= 0
    plot([0.95 0.95],[-0.3 -0.1],'g','linewidth',2)
end
if EquationParas(NotEmptySets(SortedIndex(i)),EquationStarts(3)+1) ~= 0
    plot([0.1 1],[1 0.1],'c','linewidth',2)
end
if EquationParas(NotEmptySets(SortedIndex(i)),EquationStarts(3)+2) ~= 0
    plot([0.1 0.9],[-0.05 -0.05],'c','linewidth',2)
end
if EquationParas(NotEmptySets(SortedIndex(i)),EquationStarts(3)+3) ~= 0
    plot([1.05 1.05],[-0.3 -0.1],'r','linewidth',2)
end

% limit plot
axis equal
ylim([-0.3,1.3])
xlim([-0.1,1.1])

% remove the axes
set(gca,'XTick',[], 'YTick', [])

% add x label with monomials
xlabel(string(SortedGBeqns(i)))
end
