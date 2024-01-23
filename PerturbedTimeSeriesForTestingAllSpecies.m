% This script simulates all networks detected as having perfect resilience
% under a particular stimulus regime, in order to test the quality of RPA 
% for all variables. This uses the same requirements as those used for the
% target species, namely a species must react to a change in stimulus, and
% then return within some accepted region of its original abundance.

% define symbols
syms r1 r2 r3 r4 a11 a12 a13 a14 a21 a22 a23 a24 a31 a32 a33 a34 a41 a42 a43 a44 d1 d2 d3 d4 M N I S O

% load in datasets
SaveNetworkName = '3SpeciesITarget'; % IO, I, O
eval(['load Data\AcceptedEqnsAuto' SaveNetworkName])
eval(['load Data\InitialStableSystemsAuto' SaveNetworkName]) 
eval(['load Data\RPAquality' SaveNetworkName]) 
eval(['load Data\GroebnerBases' SaveNetworkName  ' EquationStarts FullEquation InputTargets NumSpp '])

% initialise storage
SearchIndex = NotEmptySets(RPAquality>=0.1);
RPAqualityAll = zeros(NumSpp,size(SearchIndex,2));

% loop over networks
for MainIndex = 1:length(SearchIndex)
% pick a set of equations
index = SearchIndex(MainIndex);

SInputs = [1 2 4]; % input difference must be 1 i.e. (S(i+1)-S(i))/S(i)=1
NumInputs = length(SInputs);

% define matrix equations and parameters
rInd = EquationStarts';
AInd = zeros(NumSpp); dInd = zeros(size(InputTargets));
for i = 1:NumSpp
    AInd(i,:) = (EquationStarts(i)+1):(EquationStarts(i)+NumSpp);
    if i<=length(InputTargets)
        if InputTargets(i)==NumSpp
            dInd(i) = length(FullEquation);
        else
            dInd(i) = EquationStarts(i+1)-1;
        end
    end
end

paras = [];
for i = 1:NumSpp
    paras = [paras str2sym(['r' num2str(i)])];
    for j = 1:NumSpp
        paras = [paras str2sym(['a' num2str(i) num2str(j)])];
    end
    if ~isempty(find(InputTargets==i,1))
        paras = [paras str2sym(['d' num2str(i)])];
    end
end

% loop over parameter sets
CurrentParaSets = SavedParameterSets{index,1};
NonZeroParas = any(CurrentParaSets~=0);
time = [0 1e3];
jloop = tic;
for j = 1:size(CurrentParaSets,1)
    % convert to matrix problem
    r = CurrentParaSets(j,rInd)';
    A = reshape(CurrentParaSets(j,AInd),NumSpp,NumSpp); A(1:NumSpp+1:end) = A(1:NumSpp+1:end)*-1;
    d = zeros(NumSpp,1); d(InputTargets)=CurrentParaSets(j,dInd);
    
    % run simulation (need init to be a feasible stable steady state)
    uinit = RPASteadyStates{index}(j,:); %disp(uinit), disp([r,d,A])
    [t,unew] = ode45(@(t,u) odesys(t,u,r,d,SInputs(1),A), time, uinit);
    uinit = unew(end,:);
    SimulationData = [t unew];
    
    % step over inputs
    for k = 2:NumInputs
        % run simulation
        [t,unew] = ode23s(@(t,u) odesys(t,u,r,d,SInputs(k),A), time, uinit);
        % calculate measurements of RPA
        SimulationData = [SimulationData; t+SimulationData(end,1) unew];
        % determine if output is following RPA conditions
        RPAqualityAll(:,MainIndex) = RPAqualityAll(:,MainIndex) + [((unew(end,:)-uinit)./uinit<=0.01 & (max(unew-uinit,[],1))./uinit>=0.01)]';
        % update init
        uinit = unew(end,:);
    end

end
% calculate the quality and store
RPAqualityAll(:,MainIndex) = RPAqualityAll(:,MainIndex)./(size(CurrentParaSets,1).*(NumInputs-1));
disp(['Quality of RPA determined for ' num2str(round(MainIndex/length(SearchIndex)*100)) '% of networks'])
end

% plot a rough figure of RPA qualities
figure
imagesc(RPAqualityAll)
axis equal

% save data
eval(['save data\RPAqualityAll' SaveNetworkName ' RPAqualityAll'])

% function for ODE
function eqn = odesys(t,u,r,d,S,A)
% input: t,u,r,d,S,A
% output: eqn
% use matrix multiplication to advance solution of ode system. sets
% negative abundances to zero
u(u<0) = 0; u(u>100) = 100;
if any(isnan(u)), keyboard, end
eqn = (r+d.*S).*u + (A*u).*u;
end
