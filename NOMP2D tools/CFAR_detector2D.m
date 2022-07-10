function [T_judgement, Threshold_CUT] = CFAR_detector2D(y_matrix, train_guard_cell, p_fa_CFAR, overSamplingRate, CFAR_method)

    % Nx = length(y_matrix);
    [Nx, My] = size(y_matrix);
    NM_num = Nx * My;

    OSR = overSamplingRate;
    % R = Nx * overSamplingRate;
    Nx_ext = Nx * overSamplingRate(1);
    My_ext = My * overSamplingRate(2);
    guard_band = train_guard_cell(1 : 2);
    training_band = train_guard_cell(3 : 4);

    % gains  = fft(y_matrix, R) / sqrt(Nx);   % R * T matrix
    % prob = sum(abs(gains).^2,2);
    % Oversam_vec = [Nx_ext, My_ext];    % the size vector after oversampling
    % 2D-fft and find the maximum index
    Y_spec = fft2(y_matrix, Nx_ext, My_ext) / sqrt(NM_num);
    prob = abs(Y_spec) .^ 2;

    prob_ind = prob(1 : OSR : end, 1 : OSR : end);
    [res_inf_normSq_rot, peak_idx] = max(prob_ind(:));
    [peak_n, peak_m] = ind2sub([Nx, My], peak_idx);

    prob_ind_ext = repmat(prob_ind, 3, 3);
    cfar_detector = phased.CFARDetector2D('TrainingBandSize', training_band - guard_band, 'GuardBandSize', guard_band);
    cfar_detector.Method = CFAR_method;
    if strcmp(CFAR_method, 'OS')
        N_r = (2 * training_band(1) + 1) * (2 * training_band(2) + 1) - (2 * guard_band(1) + 1) * (2 * guard_band(2) + 1); 
        cfar_detector.Rank = round(N_r / 2);
    end
    cfar_detector.ThresholdFactor = 'Auto';
    cfar_detector.ProbabilityFalseAlarm = p_fa_CFAR;
    cfar_detector.ThresholdOutputPort = true;
    cfar_detector.NoisePowerOutputPort = true;
    [~, Threshold_CUT, ~] = cfar_detector(prob_ind_ext, [peak_n + Nx; peak_m + My]);
    T_judgement = res_inf_normSq_rot / Threshold_CUT - 1;

end