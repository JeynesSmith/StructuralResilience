function [KnownMonomials,Monomials] = MonomialAnalysis(SaveNetworkName,GroebnerBasisEquations,lengthcombos)

% preallocate space
KnownMonomials = GroebnerBasisEquations{1,2};
Monomials = zeros(lengthcombos,length(KnownMonomials));


% loop over combinations, find the monomials and add them to the matrix
for i = 1:lengthcombos
    
    % first find the number of monomials in a given equation
    tempequation = GroebnerBasisEquations{i,2};
    NumMonomials = length(tempequation);
    
    % loop over equation and find monomials
    for j = 1:NumMonomials
        
        % find the location of the monomial in the list of the monomials
        location = find(KnownMonomials==tempequation(j));
        
        if isempty(location) 
            % if the monomials is not known, add it to the list and add a
            % column to the matrix
            KnownMonomials = [KnownMonomials tempequation(j)];
            Monomials = [Monomials zeros(lengthcombos,1)];
            Monomials(i,end) = 1;
        else
            % if the monomial is known, place a 1 in monomials
            Monomials(i,location) = 1;
        end
        
    end
    
    disp([num2str(i/lengthcombos) '% equations analysed for monomials'])
    
end

eval(['save Monomials' SaveNetworkName ' KnownMonomials Monomials'])

end