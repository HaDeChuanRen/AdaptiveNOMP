function [omegaList, gainList, residueList] = MNOMP(y, S,...
			      	   tau, overSamplingRate, R_s , R_c)
% code is built based on NOMP 
% SUMMARY:
% 
%   given measurements: y = S * (mixture of sinusoids) + white noise
%          and a parameter tau which relates reduction in residue to
%          sparsity of the explanation (number of sinusoids), 
%          this module returns two **ordered** lists: a list of 
%          estimates of frequencies of the sinusoids in the mixture 
%          and and another list of corresponding gains 
% INPUT:
%    y - measurements
%    S - measurement matrix: can be compressive, or
%        just the identity matrix
%    tau - algorithm parameter which determines what the minimum payoff
%          we expect per sinusoid
%        - should be slightly more than the noise level sigma2
%        We solve for gains and continuous frequencies which minimize:
%        minimize norm(y - sum of sinusoids)^2 + tau * ell0_norm(gain)
%        ** LARGE VALUES OF TAU PROMOTE SPARSITY **% 
%    overSamplingRate (optional) -  used to determine how finely 
%              we sample frequency grid on which we precompute
%              responses when measurements are compressive (using 
%              IFFTs) how fine (in terms of multiples of the FFT 
%              grid) do we want the coarse grid to be?
%              number of grid points = overSamplingRate * N
%              Default value is 4
% 
%    R_s (optional) : number of Newton steps for a Single frequency update
%                           Default value is 1
%    R_c (optional) : number of rounds of Cyclic Refinement for all of the
%                           already detected sinusoids
%                           Default value is 3
% 
% OUTPUT:
%   omegaAns    - frequencies
%   gainAns     - gains of estimated frequencies
%   residueList - trajectory of the energy in the residual
%                 measurements as we add new sinsoids

if ~exist('overSamplingRate','var'), overSamplingRate = 4;
elseif isempty(overSamplingRate), overSamplingRate = 4; end

if ~exist('R_s','var'), R_s = 1;
elseif isempty(R_s),    R_s = 1; end

if ~exist('R_c','var'), R_c = 3;
elseif isempty(R_c),    R_c = 3; end

% Algorithm preprocessing
sampledManifold = preProcessMeasMat(S, overSamplingRate);

if sampledManifold.is_eye
	S = []; % erase identity matrix
end

omegaList = [];
gainList  = [];
y_r = y;

residueList = [ y_r(:)' * y_r(:) ];

while true
    
    % keep detecting new sinusoids until power in residue 
    % becomes small; *** how small *** determined by *** tau ***
    
    % detect gain and frequency of an additional sinusoid
    [omega_new, gain_new, y_r, res_inf_normSq_rot] = ...
        detectNew(y_r, sampledManifold);
    % detecttNew removes the contribution of the newly detected
    % from the old residue  y_r(input) and reports the new residual
    % measurement y_r (output)
    % stopping criterion:
    if res_inf_normSq_rot < tau
        break;
    end
    
    % newly detected sinusoid is coarse - so we refine it to 
    % imitate detection on the continuum
    for i = 1:R_s
        [omega_new, gain_new, y_r] = refineOne(y_r, omega_new, ...
            gain_new, S, sampledManifold.ant_idx, true);
    end
    % refineOne checks whether the refinement step decreases
    % the l-2 norm of the residue y_r
    
    % Add newly detected sinusoid to the ordered lists
    omegaList = [omegaList; omega_new];
    gainList  = [gainList; gain_new];
    
    % refine all frequencies detected so far
    % can be interpreted as a search for better frequency supports
    [omegaList, gainList, y_r] = refineAll(y_r, omegaList,...
        gainList, S, sampledManifold.ant_idx, R_s, R_c);
    % refineAll only uses refineOne to tweak parameters and the energy 
    % in the residual measurements y_r can only decrease as a result

%     % Solve least squares for the dictionary set [Ax(omega)] omega in 
%     % omegaList
    [omegaList, gainList, y_r] = solveLeastSquares(y , omegaList, ...
        S, sampledManifold.ant_idx);    
    % ensures that for the support we have settled on our choice of 
    % gains is optimal (in terms giving us the lowest residual energy)
    % though the practical value of the least-squares step is debatable
    % when the frequency support is well-conditioned (since least-squares
    % is a by-product of refineAll once the frequencies have converged)
    % we need this step for theoretical guarantees on convergence rates
 
    residue_new = y_r(:)'*y_r(:);
    residueList = [residueList; residue_new];
    
end

% revert to standard notion of sinusoid: 
%           exp(1j*(0:(N-1))'*omega)/sqrt(N)
gainList = bsxfun(@times,gainList,exp(1j*sampledManifold.ant_idx(1)*omegaList));

% gainList = gainList .* exp(1j*sampledManifold.ant_idx(1)*omegaList);
omegaList = wrap_2pi(omegaList);

end