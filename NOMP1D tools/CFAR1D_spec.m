function [T_judgement, Threshold_CUT] = CFAR1D_spec(prob, sampledManifold, training_cells, alpha_set, omegaList_new)
% last update date: 2022.6.30
% 2022.6.30: 1.separate the function from NOMP_CFAR_alpha.m

% initialize the parameters
R = length(sampledManifold.coarseOmega);
N = sampledManifold.length;
OSR = round(R/N);
% training_n = train_guard_cell(2);
guard_cells = 4;


% create the prob vector
% if sampledManifold.is_eye
%     gains  = fft(y, R)/sqrt(N);   % R * T matrix
%     if sampledManifold.ant_idx(1)~=0
%         gains = bsxfun(@times,gains,exp(-1j*sampledManifold.coarseOmega(:)*...
%             sampledManifold.ant_idx(1)));
% %         gains = gains.*exp(-1j*sampledManifold.coarseOmega(:)*...
% %             sampledManifold.ant_idx(1));
%     end
%     prob = sum(abs(gains).^2,2);
% else
%     energy = sampledManifold.map_IfftMat_norm_sq.';   % R*1
%     gains = bsxfun(@rdivide,sampledManifold.map_IfftMat'*y,energy); % R*T
% %     gains = (sampledManifold.map_IfftMat'*y)./energy;
%     prob = sum(bsxfun(@times,abs(gains).^2, energy),2);
% %     prob = (abs(gains).^2).*energy;
% end

% [~, T_snap] = size(y);

% find the index of peak cell as the cell under test
Allidx_set = 1 : N;
Target_idx = round(omegaList_new * N / (2 * pi) + 1);

prob_ind = prob(1 : OSR : end);
[res_inf_normSq_rot, peak_idx] = max(prob_ind);
guardidx_set = (peak_idx - guard_cells) : (peak_idx + guard_cells);
% guardidx_set = setdiff(guardidx_set, [-N : 0, (N + 1) : (2 * N)]);


% select the taining cells and calculate the estimation of noise variance
% tempidx_set = setdiff(Allidx_set, Targetidx_collect);
tempidx_set = setdiff(Allidx_set, Target_idx);
noiseidx_set = setdiff(tempidx_set, guardidx_set);

[~, peakidx_new] = min(abs(noiseidx_set - peak_idx));
remaining_len = length(noiseidx_set);
noiseidx_ext = repmat(noiseidx_set, [1, 3]);
tempidx = round(training_cells / 2);
% trainingidx_set = peakidx_new - round

noiseidx_sel = noiseidx_ext(remaining_len + peakidx_new - tempidx + (1 : training_cells));
sigma_hat = mean(prob_ind(noiseidx_sel));

Threshold_CUT = alpha_set * sigma_hat;
T_judgement = log(res_inf_normSq_rot / Threshold_CUT);


end
