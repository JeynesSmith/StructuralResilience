function [EquationParas,index,Reason,Factors] = AutomaticRPAIdentification(SaveNetworkName,Monomials,KnownMonomials,FullEquation,GroebnerBasisEquations,Species,SecondaryGBEquations)
% This function automatically detects whether the grobner bases can be
% factorised into the RPA-polynomial form, p(S)*(O-k). 
syms S

% allocate storage
index = []; 
Factors = cell(1,size(Monomials,1));
Reason = zeros(size(Monomials,1),1);

% Loop over equations
for i = 1:size(Monomials,1)
    % define monomial equation
    TempMonomial = sum(KnownMonomials.*Monomials(i,:));
    disp(i)
    % remove anything that isn't just a function of output species and
    % stimulus
    if has(TempMonomial,Species(end)) && ~has(TempMonomial,Species(1:end-1))
        % check if 2nd equation is a function of output and stimulus,
        % remove if so
        if all(SecondaryGBEquations{i}==0) || any(has(SecondaryGBEquations{i},Species(1:end-1)))
            % factorise the groebner basis
            TempFactors = factor(sum(GroebnerBasisEquations{i,2}.*GroebnerBasisEquations{i,3}));
            disp(TempFactors)
            % loop over factors and check for a linear term with output and
            % constant solution
            pass = 0;
            for j = 1:length(TempFactors)
                TempCoef = coeffs(TempFactors(j),Species(end),'All');
                if length(TempCoef)==2
                    % if there is a constant solution, save the factor and
                    % pass the groebner basis
                    if all(~has(TempCoef,S)) && TempCoef(2)~=0
                        Factors{i} = [Factors{i} TempCoef(2)/TempCoef(1)];
                        disp(TempFactors(j))
                        pass = 1;
                    end
                else
                    % fail 1
                    Reason(i) = 1;
                end
            end
            if pass == 1
                index = [index i];
                Reason(i) = 0;
            else
                % there are no factors in equation
                Reason(i) = 2;
            end
        else
            % Second euation is also just a function of O and S
            Reason(i) = 3;
        end
        
    else
        % contained other variables
        Reason(i) = 4;
    end
end

%% find actual equations
syms r1 r2 r3 r4 a11 a12 a13 a14 a21 a22 a23 a24 a31 a32 a33 a34 a41 a42 a43 a44 d d1 d2 d3 d4 M N I S O real

% check that there were valid groebner bases
if ~isempty(index)
    % find equations and save
    EquationParas = sym(zeros(length(index),length(FullEquation)+1));
    for i = 1:length(index)
        EquationParas(i,:) = [GroebnerBasisEquations{index(i),1}, 0]; % the 0 is a filler so that I don't have to change other functions, but might need to do that later anyway
    end
    eval(['save Data\AcceptedEqnsAuto' SaveNetworkName ' EquationParas index Reason Factors'])
end
end