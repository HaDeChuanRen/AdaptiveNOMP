function [omegaList, gainList, y_r, A] = solveLeastSquares(y , omegaList, ...
    S, ant_idx)
% SUMMARY:
%    Reestimates all gains wholesale
N = length(ant_idx);
if isempty(S) % is this an identity matrix 
    A = exp(1j*ant_idx*omegaList.')/sqrt(N);
else
    A = S * exp(1j*ant_idx*omegaList.')/sqrt(N);
end   

% update gains
gainList = (A'*A)\(A'*y);
% update residues
y_r = y - A*gainList;

% energy in the residual measurement y_r is guaranteed to not increase 
% as a result of this operation. Therefore, we refrain from checking 
% whether the energy in y_r has increased
end