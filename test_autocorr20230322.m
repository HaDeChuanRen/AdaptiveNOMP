clc; clear; close all;

MC = 100;
% Define Scenario
Nx = 256; % Length of Sinusoid
K = 16;
sigma_n = 1;              % noise variance sigma^2, instead of sigma
S_snap = 1;
SNRvec_all = 12 : 1 : 20;
length_SNR = length(SNRvec_all);

M = Nx;
Smat_com = eye(M);

% algorithm parameters
K_max = 32;
P_oe = 0.01;
CFAR_method = 'CA';
gamma_oversamping = 4;
N_r = 60;
R_c = 3;
R_s = 1;

tau_set = sigma_n * chi2inv((1 - P_oe) ^ (1 / Nx), 2 * S_snap) / 2;
alpha_set = alpha_PoebyS(P_oe, Nx, N_r, S_snap);

cmat_corr = zeros(MC, N_r);

tic;
for mc_idx = 1 : MC
    hwaibar = waitbar(mc_idx / MC);
    [y, omega_true, gain_true] = create_yvector(K, S_snap, SNRvec_all(end), ...
        sigma_n, Smat_com);

    [omegavec_CA, gainvec_CA, y_r] = ...
        MNOMP_CFAR_alpha(y, Smat_com, alpha_set, N_r, K_max);

    sampledManifold = preProcessMeasMat(Smat_com, gamma_oversamping);

    % Question: which y_r are used for the following steps?
    % initialize the parameters
    R = length(sampledManifold.coarseOmega);
    N = sampledManifold.length;
    OSR = round(R/N);
    % training_n = train_guard_cell(2);
    guard_cells = 4;


    % create the prob vector
    if sampledManifold.is_eye
        gains  = fft(y, R)/sqrt(N);   % R * T matrix
        if sampledManifold.ant_idx(1)~=0
            gains = bsxfun(@times, gains, exp(-1j * ...
                sampledManifold.coarseOmega(:) * sampledManifold.ant_idx(1)));
        end
        prob = sum(abs(gains).^2,2);
    else
        energy = sampledManifold.map_IfftMat_norm_sq.';   % R*1
        gains = bsxfun(@rdivide,sampledManifold.map_IfftMat'*y,energy); % R*T
        prob = sum(bsxfun(@times,abs(gains).^2, energy),2);
    end


    Allidx_set = 1 : N;
    Target_idx = round(omegavec_CA * N / (2 * pi) + 1);

    prob_ind = prob(1 : OSR : end);
    [res_inf_normSq_rot, peak_idx] = max(prob_ind);
    guardidx_set = (peak_idx - guard_cells) : (peak_idx + guard_cells);

    tempidx_set = setdiff(Allidx_set, Target_idx);
    noiseidx_set = setdiff(tempidx_set, guardidx_set);

    [~, peakidx_new] = min(abs(noiseidx_set - peak_idx));
    remaining_len = length(noiseidx_set);
    noiseidx_ext = repmat(noiseidx_set, [1, 3]);
    tempidx = round(N_r / 2);
    % trainingidx_set = peakidx_new - round

    noiseidx_sel = noiseidx_ext(remaining_len + peakidx_new - tempidx + (1 : N_r));
    gain_ind = gains(1 : OSR : end);
    noise_n = gain_ind(noiseidx_sel);

    [cvec_corr, lags_corr] = xcorr(noise_n);
    zero_idx = find(lags_corr == 0);
    cmat_corr(mc_idx, :) = cvec_corr(zero_idx : end) / N_r;
    % plot(10 * log10(abs(cvec_corr) / N_r))
end
toc;
delete(hwaibar)

cvec_mean = abs(mean(cmat_corr)); % Question: abs before mean?

% plot the figure
lw = 2;
fsz = 12;
msz = 8;
lag_vec = 0 : N_r - 1;

figure;
plot(lag_vec, 10 * log10(cvec_mean), '-bo', 'LineWidth', lw, 'Markersize', msz);
xlabel('lag');
ylabel('auto-correlation (dB)');

% N_r = 100;
% noise_n = sqrt(sigma_n / 2) * (randn(N_r, S_snap) + 1j*randn(N_r, S_snap));
% [c, lags] = xcorr(noise_n);
% plot(lags, 10 * log10(abs(c) / N_r))






% % Algorithm preprocessing
% sampledManifold = preProcessMeasMat(Smat_com, gamma_oversamping);

% if sampledManifold.is_eye
% 	S = []; % erase identity matrix
% end

% omegaList = [];
% gainList  = [];
% y_r = y;

% residueList = [y_r(:)' * y_r(:)];
% Khat = K_max;

% % initialization step: assumed that the number of targets $K$ is
% % upper bounded by a known constant $K_{\rm max}$
% for cnt_k = 1 : K_max

%     % keep detecting new sinusoids until power in residue
%     % becomes small; *** how small *** determined by *** tau ***

%     % detect gain and frequency of an additional sinusoid
%     [omega_new, gain_new, y_r, ~] = detectNew(y_r, sampledManifold);
%     % detecttNew removes the contribution of the newly detected
%     % from the old residue y_r(input) and reports the new residual
%     % measurement y_r (output)
%     % stopping criterion:

%     % newly detected sinusoid is coarse - so we refine it to
%     % imitate detection on the continuum
%     for i = 1:R_s
%         [omega_new, gain_new, y_r] = refineOne(y_r, omega_new, ...
%             gain_new, S, sampledManifold.ant_idx, true);
%     end
%     % refineOne checks whether the refinement step decreases
%     % the l-2 norm of the residue y_r

%     % Add newly detected sinusoid to the ordered lists
%     omegaList = [omegaList; omega_new]; %#ok<*AGROW>
%     gainList  = [gainList; gain_new];

%     % refine all frequencies detected so far
%     % can be interpreted as a search for better frequency supports
%     [omegaList, ~, ~] = refineAll(y_r, omegaList, gainList, S,...
%         sampledManifold.ant_idx, R_s, R_c);
%     % refineAll only uses refineOne to tweak parameters and the energy
%     % in the residual measurements y_r can only decrease as a result

%     % Solve least squares for the dictionary set [Ax(omega)] omega in
%     % omegaList

%     [omegaList, gainList, y_r, A] = solveLeastSquares(y, omegaList, ...
%         S, sampledManifold.ant_idx);

%     % ensures that for the support we have settled on our choice of
%     % gains is optimal (in terms giving us the lowest residual energy)
%     % though the practical value of the least-squares step is debatable
%     % when the frequency support is well-conditioned (since least-squares
%     % is a by-product of refineAll once the frequencies have converged)
%     % we need this step for theoretical guarantees on convergence rates

%     residue_new = y_r(:)' * y_r(:);
%     residueList = [residueList; residue_new];
% end


% % detection step: implemented to perform the target detection
% while true
%     % initialize the saving array and target index set
%     y_r_det = zeros(M, S_snap, Khat);
%     Tarray_judgement = zeros(Khat, 1);
%     Threshold_collect = zeros(Khat, 1);
%     omegaList_save = zeros(Khat - 1, Khat);
%     gainList_save = zeros(Khat - 1, S_snap, Khat);
%     tar_set = 1 : Khat;

%     % calculate the criterion values and thresholds of all targets
%     if Khat > 1
%         for kidx = Khat : -1 : 1
%             % hypothesis testing problem

%             % Throwing the kidx th target and calculate the threshold
%             % save all the residuals which are calculated by each
%             % corresponding target added to the old residual
%             y_r_det(:, :, kidx) = y_r + A(:, kidx) * gainList(kidx, :);
%             y_r_det_sq = squeeze(y_r_det(:, :, kidx));
%             tar_set_diff = setdiff(tar_set, kidx);

%             % Refine the residual targets and solve least squares for them
%             [omegaList0, ~, ~] = refineAll(y_r_det_sq,...
%                 omegaList(tar_set_diff), gainList(tar_set_diff,:), S,...
%                 sampledManifold.ant_idx, R_s, R_c);
%             [omegaList_new, gainList_new, y_r0, ~] = ...
%                 solveLeastSquares(y, omegaList0, S,...
%                 sampledManifold.ant_idx);
%             y_r_det(:, :, kidx) = y_r0;

%             % Test the validity of the detected targets with CFAR method
%             % Calculate the criterion values and thresholds in the CFAR
%             % test step
%             [T_judgement, Threshold_CUT] = CFAR_1D_alpha(y_r0,...
%                 sampledManifold, N_r, alpha_set, omegaList_new);

%             % collect the calculation result
%             Tarray_judgement(kidx) = T_judgement;
%             Threshold_collect(kidx) = Threshold_CUT;
%             omegaList_save(:, kidx) = omegaList_new;
%             gainList_save(:, :, kidx) = gainList_new;

%         end
%     else
%         % Consider the cases where $K = 1$ or $K = 0$
%         y_r_det = y_r + A * gainList;
%         [T_judgement, Threshold_CUT] = CFAR_1D_alpha(y_r_det,...
%             sampledManifold, N_r, alpha_set, omegaList);

%         Tarray_judgement = T_judgement;
%         Threshold_collect = Threshold_CUT;
%     end
%     % Tarray_judgement(k): reflects the existance of the kth target: >0,
%     % existence, else deactive this target
%     [~, ktar_idx] = min(Tarray_judgement);
%     if Tarray_judgement(ktar_idx) < 0
%         % If the target whose judgement value is minimum is invalid, the
%         % target will be deleted from the result list by using the
%         % correspoding new omega list
%         omegaList = omegaList_save(:, ktar_idx);
%         gainList = gainList_save(:, :, ktar_idx);
%         Khat = Khat - 1;

%         if Khat == 0
%             % when the Khat is zeros, set the result as empty
%             omegaList = [];
%             gainList = [];
%             Threshold_collect = [];
%             break;
%         end
%         y_r = squeeze(y_r_det(:, :, ktar_idx));
%         % Is the refinement useful here? Can I delete them?
%         [omegaList, ~, ~] = refineAll(y_r, omegaList,...
%         gainList, S, sampledManifold.ant_idx, R_s, R_c);
%         [omegaList, gainList, y_r, A] = solveLeastSquares(y, omegaList, ...
%         S, sampledManifold.ant_idx);
%     elseif Khat < K_max
%         % If all the targets are valid in the result list and the number
%         % of targets in the result list is smaller the K_max, we will test
%         % the residaul to find whether there exists targets.
%         [T_judgement, Threshold_CUT] = CFAR_1D_alpha(y_r, ...
%             sampledManifold, N_r, alpha_set, omegaList);

%         if T_judgement > 0 && (Khat < K_max)
%             % detect the new target and refine it
%             [omega_new, gain_new, y_r, ~] = detectNew(y_r, sampledManifold);
%             for i = 1 : R_s
%                 [omega_new, gain_new, y_r] = refineOne(y_r, omega_new, ...
%                     gain_new, S, sampledManifold.ant_idx, true);
%             end

%             % Collect the new target
%             omegaList = [omegaList; omega_new];
%             gainList  = [gainList; gain_new];
%             Threshold_collect = [Threshold_collect; Threshold_CUT];

%             % refine all frequencies detected so far
%             % can be interpreted as a search for better frequency supports
%             [omegaList, ~, ~] = refineAll(y_r, omegaList,...
%                 gainList, S, sampledManifold.ant_idx, R_s, R_c);
%             % refineAll only uses refineOne to weak the energy in the
%             % residual measurements y_r which can only decrease as a result

%             % Solve least squares for the dictionary set [Ax(omega)] omega
%             % in omegaList
%             [omegaList, gainList, y_r, A] = solveLeastSquares(y,...
%                 omegaList, S, sampledManifold.ant_idx);
%         else
%             break;
%         end
%     else
%         break;
%     end

% end

% % revert to standard notion of sinusoid:
% %           exp(1j*(0:(N-1))'*omega)/sqrt(N)
% if ~isempty(omegaList)
%     gainList = bsxfun(@times, gainList,...
%     exp(1j*sampledManifold.ant_idx(1)*omegaList));

%     % gainList = gainList .* exp(1j*sampledManifold.ant_idx(1)*omegaList);
%     omegaList = wrap_2pi(omegaList);
% end