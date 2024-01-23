function [GroebnerBasisEquations,Species,InputTargets,FullEquation,EquationStarts,lengthcombos,NumSpp,KnownMonomials,Monomials,SecondaryGBEquations] = ExaminingEcosystemGBLoop(Species,InputTargets,SaveNetworkName)
% This function identifies the full 

% initialise variables
syms r1 r2 r3 r4 a11 a12 a13 a14 a21 a22 a23 a24 a31 a32 a33 a34 a41 a42 a43 a44 d1 d2 d3 d4 M N I S O real

% 1. define species and input targets
% automatically define gLV equations and parameters
NumSpp = length(Species);
FullEquation = []; paras = []; 
for i = 1:NumSpp
    FullEquation = [FullEquation str2sym(['r' num2str(i)])*Species(i)];
    paras = [paras str2sym(['r' num2str(i)])];
    for j = 1:NumSpp
        paras = [paras str2sym(['a' num2str(i) num2str(j)])];
        if i==j
            FullEquation = [FullEquation str2sym(['a' num2str(i) num2str(j)])*Species(i)*Species(j)*(-1)];
        else
            FullEquation = [FullEquation str2sym(['a' num2str(i) num2str(j)])*Species(i)*Species(j)];
        end
    end
    if ~isempty(find(InputTargets==i,1))
        FullEquation = [FullEquation str2sym(['d' num2str(i)])*Species(i)*S];
        paras = [paras str2sym(['d' num2str(i)])];
    end
end

% define the starting point of every equation
EquationStarts = ones(1,NumSpp);
for i = 1:(NumSpp-1)
    EquationStarts(i+1) = EquationStarts(i)+NumSpp+1+~isempty(find(InputTargets==i,1));
end
% 2. find all possible combinations of terms in equations based on full
% factorial design
Combos = ff2n(NumSpp*(NumSpp+1));
lengthcombos = size(Combos,1);
AllCombos = [];
for i = 1:NumSpp
    AllCombos = [AllCombos Combos(:,(1:NumSpp+1)+(i-1)*(NumSpp+1)) ones(lengthcombos,~isempty(find(InputTargets==i,1)))];
end

% 3. calculate the groebner bases by mutliplying a combination of 
% parameters by the full equation and then defining a vector of equations

% initialise data storage variables
GroebnerBasisEquations = cell(lengthcombos,3);
SecondaryGBEquations = cell(lengthcombos,1);
vars = [Species(1:(end-1)) S Species(end)];

% loop over combinations
for i = 1:lengthcombos
    
    % define variables which we keep and define equations
    AlteredEquation = FullEquation.*AllCombos(i,:);
    p = sym(zeros(1,NumSpp));
    for j = 1:(NumSpp-1)
        p(j) = sum(AlteredEquation(EquationStarts(j):(EquationStarts(j+1)-1)));
    end
    p(end) = sum(AlteredEquation(EquationStarts(end):end));
    
    % calculate Grobner Basis
    GB = gbasis(p,vars,'MonomialOrder','lexicographic');
    
    % split into coefficients and monomials
    [tc,tm] = coeffs(GB(end),vars);
    
    % add data to storage
    GroebnerBasisEquations{i,1} = AlteredEquation; 
    GroebnerBasisEquations{i,2} = tm;
    GroebnerBasisEquations{i,3} = tc;
    
    % if there is more than 1 equation in the gorebner basis, save the
    % second equation monomials, or set it as 0
    if length(GB)>1
        [~,tempMon] = coeffs(GB(end-1),vars);
    else 
        tempMon = 0;
    end
    SecondaryGBEquations{i} = tempMon;
    
    % display an update of the progress and partially save progress
    disp([num2str(i) ' out of ' num2str(lengthcombos) ' completed'])
    if mod(i,100) == 0
        save(['Data\GroebnerBases' SaveNetworkName],'GroebnerBasisEquations','Species','InputTargets','FullEquation','EquationStarts','lengthcombos','NumSpp', '-v7.3')
    end
    
end

% save Groebner Bases
save(['GroebnerBases' SaveNetworkName],'GroebnerBasisEquations','Species','InputTargets','FullEquation','EquationStarts','lengthcombos','NumSpp','SecondaryGBEquations', '-v7.3')

% Find monomials in all equations
[KnownMonomials,Monomials] = MonomialAnalysis(SaveNetworkName,GroebnerBasisEquations,lengthcombos);

end