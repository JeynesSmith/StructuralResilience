
% define symbols
syms r1 r2 r3 r4 a11 a12 a13 a14 a21 a22 a23 a24 a31 a32 a33 a34 a41 a42 a43 a44 d1 d2 d3 d4 M N I S O

% load in datasets
SaveNetworkName = '3SpeciesOTarget'; % IO, I, O
eval(['load AcceptedEqnsAuto' SaveNetworkName])
eval(['load InitialStableSystemsAuto' SaveNetworkName ' NotEmptySets'])
eval(['load ExtensiveStableSystems' SaveNetworkName])
eval(['load RPAquality' SaveNetworkName])

refs = find(RPAquality > 0.3);

for j = 1:length(refs)
    
    % identify the present parameters
    PresentParas = find(EquationParas(NotEmptySets(refs(j)),1:end-1)~=0);
    
    % find a parameter set
    ParameterSets = SavedParameterSets{j};
    DifferentSignParas = find(~all(ParameterSets(:,PresentParas)>0,1) & ~all(ParameterSets(:,PresentParas)<0,1));
    
    
    figure, clf, tiledlayout('flow','TileSpacing','compact')
    %%% plot histograms
    for i = 1:length(PresentParas)
        nexttile
        histogram(ParameterSets(:,PresentParas((i))))
        xlabel(string(EquationParas(NotEmptySets(refs(j)),PresentParas((i)))))
    end
    %%% plot correlation matrix
    if ~isempty(DifferentSignParas)
        %         imagesc(corr(ParameterSets(:,PresentParas(DifferentSignParas))))
        nexttile
        imagesc(corr(sign(ParameterSets(:,PresentParas(DifferentSignParas)))))
        axis square
        colorbar, colormap(gray)
        ax = gca;
        ax.CLim = [-1,1];
        ax.YTick = 1:length(PresentParas(DifferentSignParas)); ax.XTick = [];
        ax.YTickLabel = string(string(EquationParas(NotEmptySets(refs(j)),PresentParas(DifferentSignParas))));
        ax.XTickLabelRotation = 90;
        xlabel('Correlation Matrix')
        %     eval(['saveas(gcf,''Figures\Histograms\Histogram' SaveNetworkName 'Network' num2str(NotEmptySets(j)) '.fig'');']);
    end
    
end