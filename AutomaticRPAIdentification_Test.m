function [EquationParas,index] = AutomaticRPAIdentification_Test(SaveNetworkName,Monomials,KnownMonomials,FullEquation,GroebnerBasisEquations,Species)
% analysis of the outputted monomial matrix to all equations which have
% equations capable of RPA. To do this, I use manual detection to identify
% what must be removed and then remove it. require that RPA is of form
% p(S)*(O-k)
syms S

%% Load in dataset
% SaveNetworkName = '3SpeciesOTarget';
% eval(['load Data\Monomials' SaveNetworkName])
% eval(['load Data\GroebnerBases' SaveNetworkName ' GroebnerBasisEquations FullEquation Species'])

%% Define bad monomials, find, and remove
index = [];
for i = 1:size(Monomials,1)
    TempMonomial = sum(KnownMonomials.*Monomials(i,:));
    disp(i)
    if has(TempMonomial,Species(end)) && ~has(TempMonomial,Species(1:end-1))
        TempPolynomial = coeffs(TempMonomial,S,'All');
        
        pass = 1; ref = 2;
        [r,q] = polynomialReduce(TempPolynomial(1),1+Species(end));
        while r == 0 && pass == 1 && ref<=length(TempPolynomial)
            [r,q2] = polynomialReduce(TempPolynomial(ref),1+Species(end));
            if r ~= 0 || ~(isequal(q,q2)||q2==0)
                pass = 0;
            end
            ref = ref+1;
        end
        if pass == 1 && r == 0 %&& TempPolynomial(end)==0 && any(factor(sum(TempPolynomial(1:end-1))) == 1+Species(end)) 
            index = [index i];
        end
        
    end
    
end

%% find actual sets
syms r1 r2 r3 r4 a11 a12 a13 a14 a21 a22 a23 a24 a31 a32 a33 a34 a41 a42 a43 a44 d d1 d2 d3 d4 M N I S O real

if ~isempty(index)
    EquationParas = sym(zeros(length(index),length(FullEquation)+1));
    for i = 1:length(index)
        EquationParas(i,:) = [GroebnerBasisEquations{index(i),1}, 0]; % the 0 is a filler so that I don't have to change other functions, but might need to do that later anyway
    end
    eval(['save AcceptedEqnsAuto' SaveNetworkName ' EquationParas index'])
end
end