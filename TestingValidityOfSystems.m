function [] = TestingValidityOfSystems(SaveNetworkName,EquationParas,Species,FullEquation,EquationStarts,InputTargets,NumSpp,Extent,QualityCalc)
% define symbols
syms r1 r2 r3 r4 a11 a12 a13 a14 a21 a22 a23 a24 a31 a32 a33 a34 a41 a42 a43 a44 d d1 d2 d3 d4 M N I S O

% 1. load in the working equations
% SaveNetworkName = '3SpeciesOTarget';
% eval(['load Data\AcceptedEqnsAuto' SaveNetworkName])
% eval(['load Data\GroebnerBases' SaveNetworkName ' GroebnerBasisEquations Species FullEquation EquationStarts InputTargets NumSpp'])

% 2. create 
% - paras are the parameters and match the order of EquationParas. Also
% create a vector which identifies the self-regulation terms for parameter
% search later on.
paras = []; SameSignParas = EquationStarts;
for i = 1:NumSpp
    paras = [paras str2sym(['r' num2str(i)])];
    for j = 1:NumSpp
        paras = [paras str2sym(['a' num2str(i) num2str(j)])];
        if i==j
            SameSignParas = [SameSignParas length(paras)];
        end
    end
    if ~isempty(find(InputTargets==i,1))
        paras = [paras str2sym(['d' num2str(i)])];
    end
end

% - MonomialEquation which is just the prime product identifier
MonomialEquation = EquationParas(:,end);

% - EquationParas has the input terms appended to the end.
EquationParas(:,end) = [];

% - NumberEquations and ParaLength (number of parameters)
[NumberEquations,ParaLength] = size(EquationParas); 

% - simpeqnvector which is simplifiedequation but summed up
% - fulleqnvector which is fullequation but summed up
simpeqnvector = sym(zeros(NumSpp,1));
fulleqnvector = sym(zeros(NumSpp,1));
for j = 1:(NumSpp-1)
    simpeqnvector(j) = simplify(sum(FullEquation(EquationStarts(j):(EquationStarts(j+1)-1))./Species(j)));
    fulleqnvector(j) = sum(FullEquation(EquationStarts(j):(EquationStarts(j+1)-1)));
end
simpeqnvector(end) = simplify(sum(FullEquation(EquationStarts(end):end)./Species(NumSpp)));
fulleqnvector(end) = sum(FullEquation(EquationStarts(end):end));

% - Jacobian
Jacobian = sym(zeros(NumSpp));
for i = 1:NumSpp
    Jacobian(:,i) = diff(fulleqnvector,Species(i));
end
% - solution based on simplified equation and split into numden
CompleteSolution = solve(simpeqnvector, Species);
SolutionNumer = [];
SolutionDenom = [];
for i = 1:NumSpp
    eval(['[TempN,TempD] = numden(CompleteSolution.' char(Species(i)) ');'])
    SolutionNumer = [SolutionNumer TempN];SolutionDenom = [SolutionDenom TempD];
end
if Extent == 0
% 3. find all possible equations that have stable, feasible solutions
NumberofTests = 4e3; % used to be 1e2 in 1e4
NumberofPsets = 10;

[SavedParameterSets] = StabilityFeasilityTest(NumberofTests,NumberofPsets,NumberEquations,EquationParas,ParaLength,SolutionNumer,SolutionDenom,paras,SameSignParas,Jacobian,Species);

NotEmptySets = [];
for i = 1:NumberEquations
    if ~all(SavedParameterSets{i}==0,'all') && ~isempty(SavedParameterSets{i})
        NotEmptySets = [NotEmptySets i];
    end
end
disp(NotEmptySets)

save(['InitialStableSystemsAuto' SaveNetworkName],'NotEmptySets','SavedParameterSets'); 

elseif Extent == 1
% 4. find a large number of regimes for equations that have stable,
% feasible solutions
eval(['load InitialStableSystemsAuto' SaveNetworkName])
if QualityCalc == 1
    RPAquality = PerturbedTimeSeriesV2(NotEmptySets,SavedParameterSets,FullEquation,EquationStarts,InputTargets,NumSpp,Species);
    save(['RPAquality' SaveNetworkName],'RPAquality');
else
    load(['RPAquality' SaveNetworkName],'RPAquality');
end

NumberofTests = 1e10; % boosting
NumberofPsets = 2000;

[SavedParameterSets] = StabilityFeasilityTest(NumberofTests,NumberofPsets,length(NotEmptySets(RPAquality>0.3)),EquationParas(NotEmptySets(RPAquality>0.3),:),ParaLength,SolutionNumer,SolutionDenom,paras,SameSignParas,Jacobian,Species);

save(['ExtensiveStableSystems' SaveNetworkName],'SavedParameterSets'); 
end
end