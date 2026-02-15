function cvec_corr = traincell_set(y, sampledManifold, N_r, omegavec_new)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
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
    Target_idx = round(omegavec_new * N / (2 * pi) + 1);

    prob_ind = prob(1 : OSR : end);
    [~, peak_idx] = max(prob_ind);
    guardidx_set = (peak_idx - guard_cells) : (peak_idx + guard_cells);

    tempidx_set = setdiff(Allidx_set, Target_idx);
    noiseidx_set = setdiff(tempidx_set, guardidx_set);

    [~, peakidx_new] = min(abs(noiseidx_set - peak_idx));
    remaining_len = length(noiseidx_set);
    noiseidx_ext = repmat(noiseidx_set, [1, 3]);
    tempidx = round(N_r / 2);

    noiseidx_sel = noiseidx_ext(remaining_len + peakidx_new - tempidx + (1 : N_r));
    gain_ind = gains(1 : OSR : end);
    noise_n = gain_ind(noiseidx_sel);

    [cvec_corr, ~] = xcorr(noise_n, 'unbiased');
    % zero_idx = find(lags_corr == 0);
    cvec_corr = cvec_corr(N_r : end);
    % plot(10 * log10(abs(cvec_corr) / N_r))
end