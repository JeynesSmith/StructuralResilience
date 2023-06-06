function [SavedParameterSets] = StabilityFeasilityTest(NumberofTests,NumberofPsets,NumberEquations,EquationParas,ParaLength,SolutionNumer,SolutionDenom,paras,SameSignParas,Jacobian,Species)
syms r1 r2 r3 r4 a11 a12 a13 a14 a21 a22 a23 a24 a31 a32 a33 a34 a41 a42 a43 a44 d d1 d2 d3 d4 M N I S O
% Initialise SavedParameterSets
SavedParameterSets = cell(NumberEquations,1);

tic
% loop over equations
for i = 1:NumberEquations
    
    % find non-zero and zero terms in function
    nonzero = find(EquationParas(i,:));
    zero = 1:ParaLength; zero(nonzero) = [];
    
    % substitute zeros into solutions
    SubbedSolutionNumer = subs(SolutionNumer,[paras(zero) S],[zeros(size(zero)) 1]);
    SubbedSolutionDenom = subs(SolutionDenom,[paras(zero) S],[zeros(size(zero)) 1]);
    
    index1 = [];
    for j = 1:size(SubbedSolutionDenom,1)
        if isempty(find([SubbedSolutionDenom(j,:),SubbedSolutionNumer(j,:)]==0,1))
            index1 = [index1 j];
        end
    end
    
    if ~isempty(index1)
        
        % pre-allocate space for valid parameter sets
        Psets = zeros(NumberofPsets,ParaLength);
        CompletedTests = 0;
        SuccessfulTests = 0;
        
        % Test random parameter sets
        while SuccessfulTests < NumberofPsets && CompletedTests < NumberofTests
            % randomly generate parameters
            k = (rand(1,ParaLength)-0.5)*2;
            k(SameSignParas) = rand(1,length(SameSignParas)); 
            k(zero) = zeros(size(zero));
            
%             SubbedEqns = subs(SubbedSolutionNumer./SubbedSolutionDenom,paras(nonzero),k(nonzero));
            SubbedEqns = subs(SubbedSolutionNumer./SubbedSolutionDenom,paras,k);
            
            index2 = [];
            for j = 1:length(index1)
                if all((real(SubbedEqns(index1(j),:))>0) & (imag(SubbedEqns(index1(j),:)) == 0)) %%%% Maybe introduce a tolerance
                    index2 = [index2 index1(j)];
                end
            end
            
            if ~isempty(index2)
                for j = 1:length(index2)
                    SubbedJac = subs(Jacobian,[Species,paras(zero), S, paras(nonzero)],[SubbedEqns(index2(j),:), zeros(size(zero)),1,k(nonzero)]);
                    Evals = eig(SubbedJac);
                    if all(real(double(Evals))<0)
                        SuccessfulTests = SuccessfulTests + 1;
                        Psets(SuccessfulTests,:) = k;
                        disp(['Got one: ' num2str(SuccessfulTests/NumberofPsets*100) '%'])
                    end
                end
            end
            
            CompletedTests = CompletedTests + 1;
            
        end
        
        SavedParameterSets{i} = Psets;
        disp(CompletedTests)
        
    end
    
    disp(i), toc
    
end

end