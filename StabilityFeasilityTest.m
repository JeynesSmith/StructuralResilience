function [SavedParameterSets,RPASteadyStates] = StabilityFeasilityTest(NumberofTests,NumberofPsets,NumberEquations,EquationParas,ParaLength,SolutionNumer,SolutionDenom,paras,SameSignParas,Jacobian,Species,Factors)
% This function identifies parameter sets which enable feasible and stable
% steady states for a set of equations. It requires the factors found in
% the automatic search. This function outputs the parameter sets and can
% output the steady states assoicated with those parameters.

% define variables
syms r1 r2 r3 r4 a11 a12 a13 a14 a21 a22 a23 a24 a31 a32 a33 a34 a41 a42 a43 a44 d d1 d2 d3 d4 M N I S O
% Initialise SavedParameterSets
SavedParameterSets = cell(NumberEquations,1);
RPASteadyStates = cell(NumberEquations,1);

tic
% loop over equations
for i = 1:NumberEquations

    % find non-zero and zero terms in function 
    nonzero = find(EquationParas(i,:));
    zero = 1:ParaLength; zero(nonzero) = [];

    % substitute zeros into solutions
    SubbedSolutionNumer = subs(SolutionNumer,[paras(zero) S],[zeros(size(zero)) 1]);
    SubbedSolutionDenom = subs(SolutionDenom,[paras(zero) S],[zeros(size(zero)) 1]);
    [PossibleSetpoint,PossibleIndex] = ismember(-1.*Factors{i},simplify(SubbedSolutionNumer(:,end)./SubbedSolutionDenom(:,end)));
    
    % find the solutions which match the groebner basis factors
    index1 = [];
    PossibleIndex = PossibleIndex(PossibleIndex~=0);
    if isempty(PossibleIndex)
        disp('No possible factor solutions')
    end
    % check that solutions for all species are non-zero
    for j = 1:length(PossibleIndex)
        if isempty(find([SubbedSolutionDenom(PossibleIndex(j),:),SubbedSolutionNumer(PossibleIndex(j),:)]==0,1))
            index1 = [index1 PossibleIndex(j)];
        end
    end
    
    % if there are non-zero parameter sets, continue to check parameters
    if ~isempty(index1)

        % pre-allocate space for valid parameter sets and initialise checks
        Psets = zeros(NumberofPsets,ParaLength);
        SStates = zeros(NumberofPsets,length(SubbedSolutionDenom(1,:)));
        CompletedTests = 0;
        SuccessfulTests = 0;
        CheckFactorSign = 0;
        CheckSolutionFeasible = 0;
        CheckSolutionStable = 0;

        % Test random parameter sets
        while SuccessfulTests < NumberofPsets && CompletedTests < NumberofTests
            % randomly generate parameters
            k = (rand(1,ParaLength)-0.5)*2;
            k(SameSignParas) = rand(1,length(SameSignParas));
            k(zero) = zeros(size(zero));
            
            % check that there are positive factors for target species
            if any(double(subs(-1.*Factors{i}(PossibleSetpoint),paras,k))>0) 

                % calculated substituted solutions
                SubbedEqns = subs(SubbedSolutionNumer./SubbedSolutionDenom,paras,k);
                
                % check that all solutions are positive
                index2 = [];
                for j = 1:length(index1)
                    if all((real(SubbedEqns(index1(j),:))>0) & (imag(SubbedEqns(index1(j),:)) == 0)) %%%% Maybe introduce a tolerance
                        index2 = [index2 index1(j)];
                    end
                end
                
                % if positive solutions, check stability
                if ~isempty(index2)
                    for j = 1:length(index2)
                        SubbedJac = subs(Jacobian,[Species,paras(zero), S, paras(nonzero)],[SubbedEqns(index2(j),:), zeros(size(zero)),1,k(nonzero)]);
                        Evals = eig(SubbedJac);
                        if all(real(double(Evals))<0) % if stable, add parameters
                            SuccessfulTests = SuccessfulTests + 1;
                            Psets(SuccessfulTests,:) = k;
                            SStates(SuccessfulTests,:) = SubbedEqns(index2(j),:);
                            disp(['Got one: ' num2str(SuccessfulTests/NumberofPsets*100) '%'])
                        else
                            CheckSolutionStable=CheckSolutionStable+1; % add to tally of unstable solutions
                        end
                    end
                else
                    CheckSolutionFeasible=CheckSolutionFeasible+1; % add to tally of infeasible solutions
                end
            else
                CheckFactorSign = CheckFactorSign+1; % add to tally of solutions with infeasible factors
            end

            CompletedTests = CompletedTests + 1;

        end

        SavedParameterSets{i} = Psets;
        RPASteadyStates{i} = SStates;
        disp([CheckFactorSign, CheckSolutionFeasible, CheckSolutionStable])
        disp(CompletedTests)

    else
        disp('no non-zero solutions')
    end

    disp(i), toc

end

end