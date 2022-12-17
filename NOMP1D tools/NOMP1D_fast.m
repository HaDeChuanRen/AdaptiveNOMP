function [omegaList, gainList, residueList, Threshold_collect] =...
NOMP1D_fast(y, S, alpha_set, training_cells, K_max,...
overSamplingRate, R_s, R_c)
% last update date: 2022.7.30
% code is built based on NOMP-CFAR
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
%    alpha_set - the multiplu factor of CFAR method, which is calculated
%        by $\P_{\rm OE}$ and decided by CFAR method
%    training_cells: number of training cells (N_r)
%    K_max - upper limit of the target number in the scene
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
%   omegaList    - frequencies
%   gainList     - gains of estimated frequencies
%   residueList - trajectory of the energy in the residual
%                 measurements as we add new sinsoids
%   Threshold_collect - threshold collected with CFAR method

% initialize the parameters
if ~exist('overSamplingRate','var'), overSamplingRate = 4;
elseif isempty(overSamplingRate), overSamplingRate = 4; end

if ~exist('R_s','var'), R_s = 1;
elseif isempty(R_s),    R_s = 1; end

if ~exist('R_c','var'), R_c = 3;
elseif isempty(R_c),    R_c = 3; end

% [M, T] = size(y);

% Algorithm preprocessing
sampledManifold = preProcessMeasMat(S, overSamplingRate);

if sampledManifold.is_eye
	S = []; % erase identity matrix
end

omegaList = [];
gainList  = [];
Threshold_collect = [];
y_r = y;

residueList = y_r(:)' * y_r(:);
Khat = 0;
cfar_flag = 0;

% [T_judgement, Threshold_CUT] = CFAR_1D_alpha(y_r,...
% sampledManifold, training_cells, alpha_set, omegaList);

% initialization step: assumed that the number of targets $K$ is
% upper bounded by a known constant $K_{\rm max}$
while (Khat <= K_max) && (cfar_flag < 3)
    [T_judgement, Threshold_CUT] = CFAR_1D_alpha(y_r,...
    sampledManifold, training_cells, alpha_set, omegaList);

    [omega_new, gain_new, y_r, ~] = ...
        detectNew(y_r, sampledManifold);

    [omega_new, gain_new, y_r] = refineOne(y_r, omega_new, ...
        gain_new, S, sampledManifold.ant_idx, true);

    if T_judgement > 0
        omegaList = [omegaList; omega_new];
        gainList  = [gainList; gain_new];
        Threshold_collect = [Threshold_collect; Threshold_CUT];
        Khat = Khat + 1;
    else
        cfar_flag = cfar_flag + 1;
        continue;
    end

    [omegaList, ~, ~] = refineAll(y_r, omegaList, gainList, S,...
        sampledManifold.ant_idx, R_s, R_c);

    [omegaList, gainList, y_r, ~] = solveLeastSquares(y, omegaList, ...
        S, sampledManifold.ant_idx);

    residue_new = y_r(:)' * y_r(:);
    residueList = [residueList; residue_new];

end

if ~isempty(omegaList)
    gainList = bsxfun(@times, gainList,...
    exp(1j*sampledManifold.ant_idx(1) * omegaList));

    % gainList = gainList .* exp(1j*sampledManifold.ant_idx(1)*omegaList);
    omegaList = wrap_2pi(omegaList);
end

end