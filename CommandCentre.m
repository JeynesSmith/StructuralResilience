%%% Main Function
% This function performs an extensive search for perfect resilience in
% networks of any size for a specified method in which the stimulus affects
% the network. The function can be called to analyse networks step by step 
% or altogether by specifying the steps to run, where a 1 indicates that
% the step should be run. To specify the species in the network, they must
% be specified as a symbol and then in a vector for the species, where the
% output must be the last species listed. The species affected by the input
% should then referenced in order in InputTargets, and finally the unique
% naming convention for that input type is updated. Dependent functions
% assume that the governing equations take the form of the generalised
% Lotka-Volterra equations.

%% Specify network
syms I M N O
% species (where last species is always output (O)
Species = [I M O];
% input species (refers to index of species in Species which are affected by stimulus)
InputTargets = [1,2,3];
% defineSaveFileName
SaveNetworkName = '3SpeciesIMOTarget';
%% Specify steps to run
Step1 = 0;
Step2 = 1;
Step3P1 = 0;
Step3P2 = 0;
QualityCalc = 0;


%% Step 1 - Calculating Groebner Bases and Identify the Monomials in bases
if Step1 == 1
    [GroebnerBasisEquations,Species,InputTargets,FullEquation,EquationStarts,lengthcombos,NumSpp,KnownMonomials,Monomials,SecondaryGBEquations] = ExaminingEcosystemGBLoop(Species,InputTargets,SaveNetworkName);
else
    load(['Data\GroebnerBases' SaveNetworkName],'GroebnerBasisEquations','Species','InputTargets','FullEquation','EquationStarts','lengthcombos','NumSpp','SecondaryGBEquations');
    load (['Data\Monomials' SaveNetworkName]);
end

%% Step 2 - Identify RPA capability
if Step2 == 1
    [EquationParas,index,Reason,Factors] = AutomaticRPAIdentification(SaveNetworkName,Monomials,KnownMonomials,FullEquation,GroebnerBasisEquations,Species,SecondaryGBEquations);
else
    load(['Data\AcceptedEqnsAuto' SaveNetworkName]);
end
Factors = Factors(index); % reduce factors based on index

%% Step 3 - Check Stability and Feasibility
% run a rough search, check RPA-quality, then perform extensive search
if Step3P1 == 1
    TestingValidityOfSystems(SaveNetworkName,EquationParas,Species,FullEquation,EquationStarts,InputTargets,NumSpp,0,QualityCalc,Factors);
end
if Step3P2 == 1
    TestingValidityOfSystems(SaveNetworkName,EquationParas,Species,FullEquation,EquationStarts,InputTargets,NumSpp,1,QualityCalc,Factors);
end





