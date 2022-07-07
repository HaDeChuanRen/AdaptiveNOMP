function [T_judgement, Threshold_CUT] = CFAR_1D_alpha(y, sampledManifold, training_cells, alpha_set, omegaList_new)
% last update date: 2022.6.30
% 2022.6.30: 1.separate the function from NOMP_CFAR_alpha.m

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
        gains = bsxfun(@times,gains,exp(-1j*sampledManifold.coarseOmega(:)*...
            sampledManifold.ant_idx(1)));
%         gains = gains.*exp(-1j*sampledManifold.coarseOmega(:)*...
%             sampledManifold.ant_idx(1));
    end
    prob = sum(abs(gains).^2,2);
else
    energy = sampledManifold.map_IfftMat_norm_sq.';   % R*1
    gains = bsxfun(@rdivide,sampledManifold.map_IfftMat'*y,energy); % R*T
%     gains = (sampledManifold.map_IfftMat'*y)./energy;
    prob = sum(bsxfun(@times,abs(gains).^2, energy),2);
%     prob = (abs(gains).^2).*energy;
end

% [~, IDX] = max(prob);
% IDX_mod = mod(IDX, OSR);
%
% if IDX_mod == 0
%     IDX_mod = 4;
% end

% [~, T_snap] = size(y);

% find the index of peak cell as the cell under test
Allidx_set = 1 : N;
Target_idx = round(omegaList_new * N / (2 * pi) + 1);
% Target_set = [Target_idx-2;Target_idx-1;Target_idx;Target_idx+1;Target_idx+2];
Targetidx_collect = [Target_idx - 1; Target_idx; Target_idx + 1];
% Targetidx_collect = setdiff(Targetidx_collect, [-N : 0, (N + 1) : (2 * N)]);

prob_ind = prob(1 : OSR : end);
[res_inf_normSq_rot, peak_idx] = max(prob_ind);
guardidx_set = (peak_idx - guard_cells) : (peak_idx + guard_cells);
% guardidx_set = setdiff(guardidx_set, [-N : 0, (N + 1) : (2 * N)]);


% select the taining cells and calculate the estimation of noise variance
tempidx_set = setdiff(Allidx_set, Targetidx_collect);
noiseidx_set = setdiff(tempidx_set, guardidx_set);

[~, peakidx_new] = min(abs(noiseidx_set - peak_idx));
remaining_len = length(noiseidx_set);
noiseidx_ext = repmat(noiseidx_set, [1, 3]);
tempidx = round(training_cells / 2);
% trainingidx_set = peakidx_new - round

noiseidx_sel = noiseidx_ext(remaining_len + peakidx_new - tempidx + (1 : training_cells));
sigma_hat = mean(prob_ind(noiseidx_sel));

T_judgement = res_inf_normSq_rot/sigma_hat/alpha_set-1;
Threshold_CUT = alpha_set * sigma_hat;


% trainingidx_set = randperm(remaining_len, training_cells);

% tempidx_set = Allidx_set;

% trainingidx_set = noise_set(1 : training_cells);
% noiseidx_sel = noise_set(trainingidx_set);

% ratio = res_inf_normSq_rot/sigma_hat;
% cfar_detector.ThresholdFactor = 'Auto';
% % cfar_detector.ProbabilityFalseAlarm = p_oe;
% cfar_detector.ThresholdOutputPort = true;
% cfar_detector.NoisePowerOutputPort = true;
% [~, ~, sigma_hat] = cfar_detector(prob_ind_ext, peak_idx + N);
% T_judgement = res_inf_normSq_rot / Threshold_CUT - 1;

% alpha_vector = linspace(1, N_r, 1000);
% Cons_ln = sum(log(1 : N - 1));
% ConsD2_ln = sum(log(1 : N/2 - 1));

% poe_vector_ln = N_r + log(N_r ./ alpha_vector) + Cons_ln - 2 * ConsD2_ln -...
% N_r * log(N/2 + N_r ./ alpha_vector) + log(N);

% [~, alpha_idx] = min(abs(log(1 - p_oe) - poe_vector_ln));
% alpha_hat = alpha_vector(alpha_idx);


% alpha_hat = alpha_Poe(p_oe, N, N_r);
% N_g = 1000;
% res_alpha = zeros(N_g, 1);
% alpha_vector = linspace(1, 40, N_g);
% for n_idx = 1 : N_g
%     res_alpha(n_idx) = Ser_sum(N, N_r, alpha_vector(n_idx));
% end
% [~, alpha_idx] = min(abs(res_alpha - p_oe));

% alpha_hat = Ser_sum(N, N_r, alpha_vector(n_idx));
% Threshold_CUT = sigma_hat * alpha_set;
% T_judgement = res_inf_normSq_rot / Threshold_CUT - 1;
% 
% if T_snap > 1
%     alpha_vector = linspace(1, N_r, 1000);
% 
%     Cons_ln = sum(log(1 : N_r * T_snap - 1));
%     TempSum_lnS = sum(log((T_snap - 1) : (T_snap + N_r * T_snap - 2))) - Cons_ln + ...
%         (T_snap - 1) * log(alpha_vector ./ (alpha_vector + N_r)) + log(T_snap);
%     % for s_idx = 0 : T_snap - 1
%     %     TempSum_S = TempSum_S + prod(T_snap * N_r - 1 : s_idx + T_snap * N_r - 1) / factorial(s_idx) .* ...
%     %     (alpha_vector ./ (N_r + alpha_vector)) .^ s_idx;
%     % end
%     pfa_vector_ln = log(p_fa_CFAR) + (T_snap * N_r) * log(1 + alpha_vector / N_r) - TempSum_lnS;
%     [~, alpha_idx] = min(abs(pfa_vector_ln));
%     alpha_hat = alpha_vector(alpha_idx);
%     Threshold_CUT = sigma_hat * alpha_hat;
%     T_judgement = res_inf_normSq_rot / Threshold_CUT - 1;
% end

%
% num_training_cell = (2 * training_n + 1) - (2 * guard_cells + 1);

% training_win_n = peak_idx + N + (- training_n : training_n);
% Y_training_win = prob_ind_ext(training_win_n);
% guard_win_n = peak_idx + N + (- guard_cells : guard_cells);
% Y_guard_win = prob_ind_ext(guard_win_n);
% sigma_hat = (sum(Y_training_win) - sum(Y_guard_win)) / num_training_cell;


end
