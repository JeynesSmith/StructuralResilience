%%% Main Function
% this function can be called to analyse networks step by step or
% altogether

%% Specify network
syms I M N O
% species (where last species is always output (O)
Species = [I M O];
% input species (refers to index of species in Species which are affected by stimulus)
InputTargets = [3];
% defineSaveFileName
SaveNetworkName = '3SpeciesOTarget';
%% Specify steps to run
Step1 = 0;
Step2 = 1;
Step3P1 = 0;
Step3P2 = 0;
QualityCalc = 0;


%% Step 1 - Calculating Groebner Bases
if Step1 == 1
    [GroebnerBasisEquations,Species,InputTargets,FullEquation,EquationStarts,lengthcombos,NumSpp,KnownMonomials,Monomials] = ExaminingEcosystemGBLoop(Species,InputTargets,SaveNetworkName);
else
    load(['GroebnerBases' SaveNetworkName],'GroebnerBasisEquations','Species','InputTargets','FullEquation','EquationStarts','lengthcombos','NumSpp');
    load (['Monomials' SaveNetworkName]);
end

%% Step 2 - Identify RPA capability
if Step2 == 1
    [EquationParas,index,Reason] = AutomaticRPAIdentification(SaveNetworkName,Monomials,KnownMonomials,FullEquation,GroebnerBasisEquations,Species);
else
    load(['AcceptedEqnsAuto' SaveNetworkName]);
end

%% Step 3 - Check Stability and Feasibility
% run a rough search, then extensive search (parameter at end)
if Step3P1 == 1
    TestingValidityOfSystems(SaveNetworkName,EquationParas,Species,FullEquation,EquationStarts,InputTargets,NumSpp,0,QualityCalc);
end
if Step3P2 == 1
    TestingValidityOfSystems(SaveNetworkName,EquationParas,Species,FullEquation,EquationStarts,InputTargets,NumSpp,1,QualityCalc);
end





