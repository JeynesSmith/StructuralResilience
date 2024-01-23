function [RPAquality] = PerturbedTimeSeries(NotEmptySets,SavedParameterSets,FullEquation,EquationStarts,InputTargets,NumSpp,Species,RPASteadyStates)
% this function is used to ensure the networks with perfect resilience also
% have a valid form of RPA that is not trivial. This requires that the
% target reacts to a change in stimuli and then returns to within some
% accepted region of its original abundance. The stimulus is increase two
% times, and the resposne of the output is measured relative to its
% original abundance. If a parameter set results in a response which passes
% the tests, then it is added to a sum, which is finally output as a
% ratio of the NumberOfParameterSetsWithPass/TotalNumberofParameterSets for
% each network. Later a check of RPAQuality>=0.1 is required for a network
% to be accepted as not trivial RPA.

% initialise
RPAquality = zeros(size(NotEmptySets));

% loop over networks
for MainIndex = 1:length(NotEmptySets)
% pick a set of equations
index = NotEmptySets(MainIndex);

SInputs = [1 2 4]; % input different must be 1 i.e. (S(i+1)-S(i))/S(i)=1
NumInputs = length(SInputs);

% define equations and parameters
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
    
    % run simulation 
    uinit = RPASteadyStates{index}(j,:); disp(uinit), disp([r,d,A])
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
        RPAquality(MainIndex) = RPAquality(MainIndex) + ((unew(end,end)-uinit(end))/uinit(end)<=0.01 && (max(unew(:,end)-uinit(end)))/uinit(end)>=0.01);
        % update init
        uinit = unew(end,:);
    end

end
RPAquality(MainIndex) = RPAquality(MainIndex)./(size(CurrentParaSets,1)*(NumInputs-1));
disp(['Quality of RPA determined for ' num2str(round(MainIndex/length(NotEmptySets)*100)) '% of networks'])
end
end

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
