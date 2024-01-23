
% define symbols
syms r1 r2 r3 r4 a11 a12 a13 a14 a21 a22 a23 a24 a31 a32 a33 a34 a41 a42 a43 a44 d1 d2 d3 d4 M N I S O

% load in datasets
SaveNetworkName = '3SpeciesITarget'; % IO, I, O
eval(['load HPCResults\AcceptedEqnsAuto' SaveNetworkName 'base'])
eval(['load HPCResults\InitialStableSystemsAuto' SaveNetworkName  'Bigger NotEmptySets'])
eval(['load HPCResults\ExtensiveStableSystems' SaveNetworkName 'Bigger'])
eval(['load HPCResults\RPAquality' SaveNetworkName  'Bigger'])
eval(['load HPCResults\GroebnerBases' SaveNetworkName  'base EquationStarts'])

refs = find(RPAquality > 0.1); % double check
tol = 0.25;
TotalNumSets = 0; % Number of networks based on parameter correlation
NumberPairs = 0;
NumberPredPrey = 0;
test1=0;test2=0;

ParaPairs = [EquationStarts(1)+2, EquationStarts(2)+1
    EquationStarts(1)+3, EquationStarts(3)+1
    EquationStarts(2)+3, EquationStarts(3)+2];

for j = 1:length(refs)
    
    % define the number of correlated and uncorrelated paras
    NumUncorr = 0;
    NumCorrGroups = 0;
    
    % identify the present parameters
    PresentParas = find(EquationParas(NotEmptySets(refs(j)),1:end-1)~=0);
    
    % find a parameter set
    ParameterSets = SavedParameterSets{j};
    DifferentSignParas = find(~all(ParameterSets(:,PresentParas)>0,1) & ~all(ParameterSets(:,PresentParas)<0,1));
    
    
    figure(1), clf, tiledlayout('flow','TileSpacing','compact')
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
        CorrMat = corr(sign(ParameterSets(:,PresentParas(DifferentSignParas)))); % calculate the correlation matrix based on sign
        imagesc(CorrMat)
        axis square
        colorbar, colormap(gray)
        ax = gca;
        ax.CLim = [-1,1];
        ax.YTick = 1:length(DifferentSignParas); ax.XTick = [];
        ax.YTickLabel = string(string(EquationParas(NotEmptySets(refs(j)),PresentParas(DifferentSignParas))));
        ax.XTickLabelRotation = 90;
        xlabel('Correlation Matrix')
        %     eval(['saveas(gcf,''Figures\Histograms\Histogram' SaveNetworkName 'Network' num2str(NotEmptySets(j)) '.fig'');']);
        
        % calculating the number of networks
        CorrParas = [];
        CorrCorrMat = CorrMat>=(1-tol) | CorrMat<=(-1+tol);
        for i = 1:length(DifferentSignParas)
            if sum(CorrCorrMat(i,:))==1 % if only the diagonal is correlated
                NumUncorr = NumUncorr + 1; % add to the number of uncorrelated parameters
            else
                CorrParas = [CorrParas, i]; % save index of correlated parameters
            end
        end
        while ~isempty(CorrParas) % repeated until all correlated parameters and seperated into groups
            CorrParas = CorrParas(~ismember(CorrParas,find(CorrCorrMat(CorrParas(1),:)))); % check for parameters which are collectively correlated (grouped together), and remove from index
            NumCorrGroups = NumCorrGroups + 1; % add 1 to the number of grouped correlated parameters
        end
        if NumCorrGroups>1
            warning('Wow! More than 1 group!') % REMOVE %%%%%%%%%%%%%%%%%
        end
    end
    
    % examine predator-prey interactions
    for i = 1:3
        if ismember(ParaPairs(i,1),PresentParas) && ismember(ParaPairs(i,2),PresentParas)
            NumberPairs = NumberPairs + 1;
            if sum(ismember(ParaPairs(i,:),PresentParas(DifferentSignParas)))==2
                tempcorr = corr(sign(ParameterSets(:,ParaPairs(i,:))));
                if tempcorr(1,2)<=(-1+tol)
                    NumberPredPrey = NumberPredPrey + 1;
                else
                    test1=test1+1;
                end
            else
                if sign(mean(sign(ParameterSets(:,ParaPairs(i,1))))) == -1*sign(mean(sign(ParameterSets(:,ParaPairs(i,2)))))
                    NumberPredPrey = NumberPredPrey + 1;
                    else
                    test2=test2+1;
                end
            end
        end
    end
    
    % calculate the number of possible networks based on assigning
    % parameter sign
    NumSets = 2^NumUncorr*2^NumCorrGroups;
    disp(NumSets)
    TotalNumSets = TotalNumSets + NumSets; % update the total number of networks for this stimulus type
    drawnow % plot now
end

disp([NumberPairs, NumberPredPrey NumberPredPrey/NumberPairs])
disp(TotalNumSets)