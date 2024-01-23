%% Define the function

%% Define equations
% define a save name
SaveNetworkName = 'SingleNetwork';

% define the symbols for variables and stimulus
syms I M O S real
% define the symbols for parameters
syms rI rM rO aII aIM aIO aMI aMM aMO aOI aOM aOO dI dM dO real

% define the variables, seperated by commas. Include the target last
species = [I M O];

% define the paras, seperated by commas.
paras = [rM aII aIM aMO aOI aOM dO];

% define the equations, each line is a different equation. Make sure that
% the order of your equations matches the same order as how you defined
% your species. Here we use the generalised Lotka-Volterra equations.
equations = [-aII*I^2 + aIM*I*M
    rM*M - aMO*M*O
    -aOI*O*I + aOM*O*M - dO*S*O];

%%%%%%%%% No input from this point on %%%%%%%%%%%
%% calculate groebner basis

% create a vector of all variables
AllVariables = [species(1:(end-1)) S species(end)];
% calculate a Groebner basis
GB = gbasis(equations',AllVariables,'MonomialOrder','lexicographic');
% split into Monomials and Coefficients for the last Groebner basis
% equations
[Coefficients,Monomials] = coeffs(GB(end),AllVariables);
% Find the second GB equation
SecondaryGBEquation = cell(1,1);
if length(GB)>1
    [~,SecondaryGBEquation{1}] = coeffs(GB(end-1),AllVariables);
else
    SecondaryGBEquation{1} = 0;
end

%% automatic detection using function

[EquationParas,~,Reason,Factors] = AutomaticRPAIdentification(SaveNetworkName,ones(1,length(Monomials)),Monomials,equations',{equations', Monomials, Coefficients},species,SecondaryGBEquation);

%% initial stable states
if Reason == 0
    [Pass, Parameters, SteadyStates] = SingleStabilityFeasibility(equations,species,Factors,paras,S);
else
    error('This network does not have perfect resilience - failed to match RPA-polynomial')
end

%% check RPA quality
% WARNING - this section can take considerable time in order to keep it
% completely general for the type of equations used.
if Pass == 1
    NumberofParas = 2; % quick search
    % NumberofParas = sum(Parameters(:,1)~=0); % extensive search
    RPAquality = CheckQuality(Parameters, NumberofParas,equations,paras,species,SteadyStates);
else
    error('This network does not have perfect resilience - failed to find parameter sets')
end

%% Final Output
if RPAquality>=0.1
    disp('This network has perfect resilience')
else
    error('This network does not have perfect resilience - failed to meet standards of perfect resilience')
end




%% Functions
% stability and  feasibility
function [Pass, Parameters, SteadyStates] = SingleStabilityFeasibility(equations,species,Factors,paras,S)
% define the number of successful parameter sets to find
NumberofParameterSets = 50;
% Find the solution
CompleteSolution = solve(equations, species);
SolutionMat = [];
for i = 1:length(species)
    eval(['SolutionMat = [SolutionMat, CompleteSolution.' char(species(i)) '];'])
end
% find solutions with factors
index1 = [];
for i = 1:length(Factors{1})
    [~,TempInd] = ismember(SolutionMat(:,end),-1.*Factors{1}(i));
    index1 = [index1 find(TempInd)];
end
% find solutions with nonzero abundances
if ~isempty(index1)
    index2 = [];
    for i = 1:length(index1)
        if isempty(find(SolutionMat(index1(i),:)==0,1))
            index2 = [index2 index1(i)];
        end
    end
else
    % failed to find solution with factor
    disp('failed to find solution with factor')
    Pass = 0;
end
% find solutions with positive abundances
if ~isempty(index2)
    % Find Jacobian
    Jacobian = jacobian(equations',species);
    % start parameterisation
    CompletedTests = 0;
    SuccessfulTests = 0;
    Parameters = zeros(NumberofParameterSets,length(paras));
    SteadyStates = zeros(NumberofParameterSets,length(species));
    while CompletedTests<4e3 && SuccessfulTests<NumberofParameterSets
        % generate a parameter set that is uniformly distributed
        k = rand(1,length(paras));
        % substitute parameters into solution
        SubbedEqns = subs(SolutionMat(index2,:),[paras S],[k 1]);
        % check if solutions are stable for parameter sets with
        % positive solutions
        index3 = find(~any(SubbedEqns<=0,2));
        for i = 1:length(index3)
            % substitute the solutions and parameters into the Jacobian
            SubbedJacobian = subs(Jacobian,[paras,S,species],[k,1,SubbedEqns(index3(i),:)]);
            % calculate the eigenvalues
            Eigenvalues = eig(SubbedJacobian);
            % check if solution is stable
            if all(real(double(Eigenvalues))<=0)
                % update the number of successful and save the solution
                % and parameters
                SuccessfulTests = SuccessfulTests+1;
                Parameters(SuccessfulTests,:) = k;
                SteadyStates(SuccessfulTests,:) = SubbedEqns(index3(i),:);
                disp(SuccessfulTests)
            end
        end
    end
    % finalise if network found successful parameters
    if all(Parameters==0,'all')
        disp('No parameters found')
        Pass = 0;
    else
        Pass = 1;
    end
else
    % failed to find nonzero solutions
    disp('Failed to find nonzero solutions')
    Pass = 0;
end
end

% Check Quality
function RPAquality = CheckQuality(Parameters, NumberofParas,equations,paras,species,SteadyStates)
RPAquality = 0;
SInputs = [1,2,4];
for i = 1:NumberofParas
    [~,unew] = ode45(@(t,u) odesys(t,u,equations,paras,species,Parameters(i,:),SInputs(1)), [0 1e3], SteadyStates(i,:));
    uinit = unew(end,:);
    for j = 2:length(SInputs)
        % run simulation
        [~,unew] = ode45(@(t,u) odesys(t,u,equations,paras,species,Parameters(i,:),SInputs(j)), [0 1e3], uinit);
        % determine if output is following RPA conditions
        RPAquality = RPAquality + ((unew(end,end)-uinit(end))/uinit(end)<=0.01 && (max(unew(:,end)-uinit(end)))/uinit(end)>=0.01);
        % update init
        uinit = unew(end,:);
    end
    disp([num2str(i) ' parameter set tested'])
end
RPAquality = RPAquality./(NumberofParas*(length(SInputs)-1));
end

% ODE system
function out = odesys(t,u,eqns,paras,vars,paravals,input)
syms S
out = double(subs(eqns,[paras vars S],[paravals u' input]));
end
