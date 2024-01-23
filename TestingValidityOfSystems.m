function [] = TestingValidityOfSystems(SaveNetworkName,EquationParas,Species,FullEquation,EquationStarts,InputTargets,NumSpp,Extent,QualityCalc,Factors)
% this function is used to identify parameter sets for networks which
% enable feasible and stable steady states. This either does a preliminary
% check (Extent==0) or an extensive parameter search (Extent==1).
% define symbols. To find parameter sets, this function calls its 
% dependant, StabilityFeasibilityTest.m. When QualityCalc==1 and Extent==1, 
% this function will call the dependant, ModifiedPerturbedTimeSeries.m. 
% This function simulates a network using the initial parameter sets to
% determine whether the output of a network will react to a stimulus change
% and then return within some accepted region of its setpoint abudnance.
% The intial search identifies valid parameter sets, the associated steady
% states for all species, and a vector identifying which networks had valid
% parameter sets.

% define symbols
syms r1 r2 r3 r4 a11 a12 a13 a14 a21 a22 a23 a24 a31 a32 a33 a34 a41 a42 a43 a44 d d1 d2 d3 d4 M N I S O

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
% - solution based on simplified equation and split into numerator and
% denominator
CompleteSolution = solve(simpeqnvector, Species);
SolutionNumer = [];
SolutionDenom = [];
for i = 1:NumSpp
    eval(['[TempN,TempD] = numden(CompleteSolution.' char(Species(i)) ');'])
    SolutionNumer = [SolutionNumer TempN];SolutionDenom = [SolutionDenom TempD];
end

% check which type of parameter search we are performing
if Extent == 0
% Find all possible equations that have stable, feasible solutions
NumberofTests = 2e4; % maximum number of searches
NumberofPsets = 50; % number required before stopping

% run the initial parameter search
[SavedParameterSets,RPASteadyStates] = StabilityFeasilityTest(NumberofTests,NumberofPsets,NumberEquations,EquationParas,ParaLength,SolutionNumer,SolutionDenom,paras,SameSignParas,Jacobian,Species,Factors);

% identify which networks had valid parameter sets found
NotEmptySets = [];
for i = 1:NumberEquations
    if ~all(SavedParameterSets{i}==0,'all') && ~isempty(SavedParameterSets{i})
        NotEmptySets = [NotEmptySets i];
    end
end
disp(NotEmptySets)

save(['Data\InitialStableSystems' SaveNetworkName],'NotEmptySets','SavedParameterSets','RPASteadyStates'); 

elseif Extent == 1
% 4. find a large number of regimes for equations that have stable,
% feasible solutions
eval(['load Data\InitialStableSystems' SaveNetworkName])

% check the quality of RPA
if QualityCalc == 1
    RPAquality = PerturbedTimeSeries(NotEmptySets,SavedParameterSets,FullEquation,EquationStarts,InputTargets,NumSpp,Species,RPASteadyStates);
    save(['Data\RPAquality' SaveNetworkName],'RPAquality');
else
    load(['Data\RPAquality' SaveNetworkName],'RPAquality');
end

% Perform extensive parameter search
NumberofTests = 1e10; 
NumberofPsets = 2000; 
NewEqns = EquationParas(NotEmptySets(RPAquality>=0.1),:);
NewFact = Factors(NotEmptySets(RPAquality>=0.1));
[SavedParameterSets,~] = StabilityFeasilityTest(NumberofTests,NumberofPsets,length(NotEmptySets(RPAquality>=0.1)),NewEqns,ParaLength,SolutionNumer,SolutionDenom,paras,SameSignParas,Jacobian,Species,NewFact);

save(['Data\ExtensiveStableSystems' SaveNetworkName],'SavedParameterSets'); 
end
end