function [omegaList, gainList, y_r] = refineAll(y_r, omegaList,...
    gainList, S, ant_idx, R_s, R_c)
% SUMMARY:
%   uses refineOne algorithm to refine frequencies and gains of
%   of all sinusoids
% INPUT:
%    y_r - residual measurement after all detected sinusoids have been
%          removed
%    omegaList - list of frequencies of known(detected) sinusoids
%    gainList  - list of gains of known(detected) sinusoids
%    S - measurement matrix
%        if [], measurements are direct i.e S = eye(N), where N is
%       length of sinusoid	
%    ant_idx - translating indexes to phases in definition of sinusoid
%    R_s - number of times each sinusoid in the list is refined
%    R_c - number of cycles of refinement fot all of frequencies that have
%       been estimated till now
%
% OUTPUT:
%       refined versions of inputs omegaList, gainList, y_r

K = length(omegaList); % number of sinusoids


% Total rounds of cyclic refinement is "R_c"
for i = 1:R_c
    
    % chose an ordering for refinement
    
    % % *random* ordering
    % order = randperm(K);
    
    % *sequential* ordering
    % order = 1:K;
    order = K : - 1 : 1;
    
    for j = 1:K
        l = order(j);
        
        
        
        % parameters of the l-th sinusoid
        omega = omegaList(l);
        gain = gainList(l,:);
            
        % refinement repeated "R_s" times per sinusoid
        for kk = 1:R_s

            % refine our current estimates of (gain, omega) of the
            % l-th sinusoid
            [omega, gain, y_r] = refineOne(y_r,...
                       omega, gain, S, ant_idx, false);
                   
        end
        
        omegaList(l) = omega;
        gainList(l,:) = gain;
        % refineOne ensures that (gain, omega) pair are altered iff
        % the energy in the residual measurements y_r does not 
        % increase
        
    end
    
end

end