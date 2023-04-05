function [omegaList, gainList, residueList, Threshold_collect] =...
NOMP1D_low(y, S, alpha_set, training_cells, K_max, flag_max...
OSR, R_s, R_c)
% last update date: 2023.2.1
% code is built based on low-complexity NOMP-CFAR
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
%              just the identity matrix
%    alpha_set - the multiplu factor of CFAR method, which is calculated
%              by $\P_{\rm OE}$ and decided by CFAR method
%    training_cells: number of training cells (N_r)
%    K_max - upper limit of the target number in the scene
%    flag_max - upper limit of number of continuous non-effictive targets
%    OSR (optional) -  used to determine how finely
%              we sample frequency grid on which we precompute
%              responses when measurements are compressive (using
%              IFFTs) how fine (in terms of multiples of the FFT
%              grid) do we want the coarse grid to be?
%              number of grid points = OSR * N
%              Default value is 4
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
if ~exist('OSR','var'), OSR = 4;
elseif isempty(OSR), OSR = 4; end

if ~exist('R_s','var'), R_s = 1;
elseif isempty(R_s),    R_s = 1; end

if ~exist('R_c','var'), R_c = 3;
elseif isempty(R_c),    R_c = 3; end

% [M, T] = size(y);

% Algorithm preprocessing
sampledManifold = preProcessMeasMat(S, OSR);

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

R = OSR * Nx;
yrvec_fft = fft(y, R) / sqrt(N);   % R * T matrix
yrvec_fft = bsxfun(@times, yrvec_fft, exp(-1j * ...
    sampledManifold.coarseOmega(:) * sampledManifold.ant_idx(1)));
% yrvec_prob = sum(abs(yrvec_fft) .^ 2, 2);

% [T_judgement, Threshold_CUT] = CFAR_1D_alpha(y_r,...
% sampledManifold, training_cells, alpha_set, omegaList);

% initialization step: assumed that the number of targets $K$ is
% upper bounded by a known constant $K_{\rm max}$
while (Khat <= K_max) && (cfar_flag < flag_max)
    yrvec_prob = sum(abs(yrvec_fft) .^ 2, 2);
    yrvec_dec = yrvec_prob(1 : OSR : end);
    [T_judgement, Threshold_CUT] = CFAR1D_spec(yrvec_dec, ...
    sampledManifold, training_cells, alpha_set, omegaList);

    % [omega_new, gain_new, y_r, ~] = ...
    %     detectNew(y_r, sampledManifold);

    [~, IDX] = max(yrvec_prob);
    omega_new = sampledManifold.coarseOmega(IDX);
    gain_new = yrvec_fft(IDX, :);

    avec_new = exp(1j * sampledManifold.ant_idx * omega_new) / Nx;
    y_r = y_r - gain_new * avec_new;

    [omega_new, gain_new, y_r] = refineOne(y_r, omega_new, ...
        gain_new, S, sampledManifold.ant_idx, true);

    if T_judgement > 0
        omegaList = [omegaList; omega_new];
        gainList  = [gainList; gain_new];
        Threshold_collect = [Threshold_collect; Threshold_CUT];
        Khat = Khat + 1;
    else
        cfar_flag = cfar_flag + 1;
        yrvec_fft = yrvec_fft - 
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