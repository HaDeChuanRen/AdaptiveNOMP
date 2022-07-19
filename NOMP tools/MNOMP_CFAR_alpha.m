function [omegaList, gainList, residueList, Threshold_collect] = MNOMP_CFAR_alpha(y, S, alpha_set, training_cells, K_max, CFAR_method, overSamplingRate, R_s, R_c)
% last refinement date: 2022.6.30
% 2022.6.30: 1. separate the function file into different parts
% 2022.6.26: 1. change the parameter training_gurad_cell as
%           training_cells (which is equal to N_r), set gurad_cells as 5
%           2. test the range of Target_set
%           3. test the setdiff function.
%           4. change some name of variables
% 2022.6.25: 1. add notes to the functions
% 2022.6.24: 1. make CAFR_method optianal
%           2. add notes to the functions
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
%    CFAR_method (optional): method of CFAR, which can be input as
%        'CA' or 'OS'
%              Default value is 'CA'
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
if ~exist('CFAR_method','var'), CFAR_method = 'CA';
elseif isempty(CFAR_method), CFAR_method = 'CA'; end

if ~exist('overSamplingRate','var'), overSamplingRate = 4;
elseif isempty(overSamplingRate), overSamplingRate = 4; end

if ~exist('R_s','var'), R_s = 1;
elseif isempty(R_s),    R_s = 1; end

if ~exist('R_c','var'), R_c = 3;
elseif isempty(R_c),    R_c = 3; end

[M,T] = size(y);

% Algorithm preprocessing
sampledManifold = preProcessMeasMat(S, overSamplingRate);

if sampledManifold.is_eye
	S = []; % erase identity matrix
end

omegaList = [];
gainList  = [];
y_r = y;

residueList = [y_r(:)' * y_r(:)];
Khat = K_max;

% initialization step: assumed that the number of targets $K$ is
% upper bounded by a known constant $K_{\rm max}$
for cnt_k = 1 : K_max

    % keep detecting new sinusoids until power in residue
    % becomes small; *** how small *** determined by *** tau ***

    % detect gain and frequency of an additional sinusoid
    [omega_new, gain_new, y_r, ~] = ...
        detectNew(y_r, sampledManifold);
    % detecttNew removes the contribution of the newly detected
    % from the old residue y_r(input) and reports the new residual
    % measurement y_r (output)
    % stopping criterion:

    % newly detected sinusoid is coarse - so we refine it to
    % imitate detection on the continuum
    for i = 1:R_s
        [omega_new, gain_new, y_r] = refineOne(y_r, omega_new, ...
            gain_new, S, sampledManifold.ant_idx, true);
    end
    % refineOne checks whether the refinement step decreases
    % the l-2 norm of the residue y_r

    % Add newly detected sinusoid to the ordered lists
    omegaList = [omegaList; omega_new]; %#ok<*AGROW>
    gainList  = [gainList; gain_new];

    % refine all frequencies detected so far
    % can be interpreted as a search for better frequency supports
    [omegaList, ~, ~] = refineAll(y_r, omegaList, gainList, S,...
        sampledManifold.ant_idx, R_s, R_c);
    % refineAll only uses refineOne to tweak parameters and the energy
    % in the residual measurements y_r can only decrease as a result

    % Solve least squares for the dictionary set [Ax(omega)] omega in
    % omegaList

    [omegaList, gainList, y_r, A] = solveLeastSquares(y , omegaList, ...
        S, sampledManifold.ant_idx);

    % ensures that for the support we have settled on our choice of
    % gains is optimal (in terms giving us the lowest residual energy)
    % though the practical value of the least-squares step is debatable
    % when the frequency support is well-conditioned (since least-squares
    % is a by-product of refineAll once the frequencies have converged)
    % we need this step for theoretical guarantees on convergence rates

    residue_new = y_r(:)' * y_r(:);
    residueList = [residueList; residue_new];
end


%model order estimation step: implemented to perform the target detection
while true
    % initialize the saving array and target index set
    y_r_det = zeros(M, T, Khat);
    Tarray_judgement = zeros(Khat, 1);
    Threshold_collect = zeros(Khat, 1);
    omegaList_save = zeros(Khat - 1, Khat);
    gainList_save = zeros(Khat - 1, T, Khat);
    tar_set = 1 : Khat;

    % calculate the criterion values and thresholds of all targets
    if Khat > 1
        for kidx = Khat : -1 : 1
            % hypothesis testing problem

            % Throwing the kidx th target and calculate the threshold
            % save all the residuals which are calculated by each
            % corresponding target added to the old residual
            y_r_det(:, :, kidx) = y_r + A(:, kidx) * gainList(kidx, :);
            y_r_det_sq = squeeze(y_r_det(:, :, kidx));
            tar_set_diff = setdiff(tar_set, kidx);

            % Refine the residual targets and solve least squares for them
            [omegaList0, ~, ~] = refineAll(y_r_det_sq,...
                omegaList(tar_set_diff), gainList(tar_set_diff,:), S,...
                sampledManifold.ant_idx, R_s, R_c);
            [omegaList_new, gainList_new, y_r0, ~] = ...
                solveLeastSquares(y, omegaList0, S,...
                sampledManifold.ant_idx);
            y_r_det(:, :, kidx) = y_r0;

            % Test the validity of the detected targets with CFAR method
            % Calculate the criterion values and thresholds in the CFAR
            % test step
            [T_judgement, Threshold_CUT] = CFAR_1D_alpha(y_r0,...
                sampledManifold, training_cells, alpha_set, omegaList_new);

            % collect the calculation result
            Tarray_judgement(kidx) = T_judgement;
            Threshold_collect(kidx) = Threshold_CUT;
            omegaList_save(:, kidx) = omegaList_new;
            gainList_save(:, :, kidx) = gainList_new;

        end
    else
        % Consider the cases where $K = 1$ or $K = 0$
        y_r_det = y_r + A * gainList;
        [T_judgement, Threshold_CUT] = CFAR_1D_alpha(y_r_det,...
            sampledManifold, training_cells, alpha_set, omegaList);

        Tarray_judgement = T_judgement;
        Threshold_collect = Threshold_CUT;
    end
    % Tarray_judgement(k): reflects the existance of the kth target: >0,
    % existence, else deactive this target
    [~, ktar_idx] = min(Tarray_judgement);
    if Tarray_judgement(ktar_idx) < 0
        % If the target whose judgement value is minimum is invalid, the
        % target will be deleted from the result list by using the
        % correspoding new omega list
        omegaList = omegaList_save(:, ktar_idx);
        gainList = gainList_save(:, :, ktar_idx);
        Khat = Khat - 1;

        if Khat == 0
            % when the Khat is zeros, set the result as empty
            omegaList = [];
            gainList = [];
            Threshold_collect = [];
            break;
        end
        y_r = squeeze(y_r_det(:, :, ktar_idx));
        % Is the refinement useful here? Can I delete them?
        [omegaList, ~, ~] = refineAll(y_r, omegaList,...
        gainList, S, sampledManifold.ant_idx, R_s, R_c);
        [omegaList, gainList, y_r, A] = solveLeastSquares(y, omegaList, ...
        S, sampledManifold.ant_idx);
    elseif Khat < K_max
        % If all the targets are valid in the result list and the number
        % of targets in the result list is smaller the K_max, we will test
        % the residaul to find whether there exists targets.
        [T_judgement, Threshold_CUT] = CFAR_1D_alpha(y_r, ...
            sampledManifold, training_cells, alpha_set, omegaList);

        if T_judgement > 0
            % detect the new target and refine it
            [omega_new, gain_new, y_r, ~] = detectNew(y_r, sampledManifold);
            for i = 1 : R_s
                [omega_new, gain_new, y_r] = refineOne(y_r, omega_new, ...
                    gain_new, S, sampledManifold.ant_idx, true);
            end

            % Collect the new target
            omegaList = [omegaList; omega_new];
            gainList  = [gainList; gain_new];
            Threshold_collect = [Threshold_collect; Threshold_CUT];

            % refine all frequencies detected so far
            % can be interpreted as a search for better frequency supports
            [omegaList, ~, ~] = refineAll(y_r, omegaList,...
                gainList, S, sampledManifold.ant_idx, R_s, R_c);
            % refineAll only uses refineOne to tweak parameters and the energy
            % in the residual measurements y_r can only decrease as a result

            % Solve least squares for the dictionary set [Ax(omega)] omega in
            % omegaList
            [omegaList, gainList, y_r, A] = solveLeastSquares(y , omegaList, ...
                S, sampledManifold.ant_idx);
        else
            break;
        end
    else
        break;
    end

end

% revert to standard notion of sinusoid:
%           exp(1j*(0:(N-1))'*omega)/sqrt(N)
gainList = bsxfun(@times,gainList,exp(1j*sampledManifold.ant_idx(1)*omegaList));

% gainList = gainList .* exp(1j*sampledManifold.ant_idx(1)*omegaList);
omegaList = wrap_2pi(omegaList);

end