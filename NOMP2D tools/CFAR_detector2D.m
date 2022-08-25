function [T_judgement, Threshold_CUT] = ...
CFAR_detector2D(y_matrix, N_r, alpha_set, overSamplingRate, omegaList, guard_band)

    % Nx = length(y_matrix);
    [Nx, My] = size(y_matrix);
    NM_num = Nx * My;


    OSR_x = overSamplingRate(1);
    OSR_y = overSamplingRate(2);
    % R = Nx * overSamplingRate;
    Nx_ext = Nx * OSR_x;
    My_ext = My * OSR_y;
    % guard_band = train_guard_cell(1 : 2);
    % training_band = train_guard_cell(3 : 4);
    % N_r = (2 * training_band(1) + 1) * (2 * training_band(2) + 1) -...
    %     (2 * guard_band(1) + 1) * (2 * guard_band(2) + 1);

    % 2D-fft and find the maximum index
    Y_spec = fft2(y_matrix, Nx_ext, My_ext) / sqrt(NM_num);
    prob = abs(Y_spec) .^ 2;
    prob_ind = prob(1 : OSR_x : end, 1 : OSR_y : end);
    prob_ind_ext = repmat(prob_ind, 3, 3);

    if isempty(omegaList)
        omega_2DindexList = [];
    else
        omega_2DindexList = round((omegaList / (2 * pi)) .* [Nx, My]);
        % omega_indexList = sub2ind([Nx, My], omega_2DindexList(:, 1), omega_2DindexList(:, 2));
    end

    [res_inf_normSq_rot, peak_idx] = max(prob_ind(:));
    [peak_n, peak_m] = ind2sub([Nx, My], peak_idx);

    semi_num = 1;
    rc_dir = [1, 0];
    s_idx = 0;
    pointer_spiral = [peak_n, peak_m];
    n_idx = 0;
    sum_energy = 0;
    while n_idx < N_r
        if s_idx >= semi_num
            s_idx = 0;
            if rc_dir(2) == 1
                semi_num = semi_num + 1;
                rc_dir = [1, 0];
            else
                rc_dir = [0, 1];
            end
        end
        pointer_spiral = pointer_spiral + (-1) ^ semi_num * rc_dir;
        s_idx = s_idx + 1;

        if abs(pointer_spiral(1) - peak_n) <= guard_band(1) && ...
            abs(pointer_spiral(2) - peak_m) <= guard_band(2)
            continue;
        end

        if min(sum(abs(pointer_spiral - omega_2DindexList), 2)) == 0
            continue;
        end
        sum_energy = sum_energy + prob_ind_ext(Nx + pointer_spiral(1), My + pointer_spiral(2));
        n_idx = n_idx + 1;

    end

    sigma_hat = sum_energy / N_r;
    % training_cells = prob_ind_ext(peak_n + Nx - training_band(1) : peak_n + Nx + training_band(1),...
    % peak_m + My - training_band(2) : peak_m + My + training_band(2));
    % guard_cells = prob_ind_ext(peak_n + Nx - guard_band(1) : peak_n + Nx + guard_band(1),...
    % peak_m + My - guard_band(2) : peak_m + My + guard_band(2));
    % sigma_hat = (sum(training_cells(:)) - sum(guard_cells(:))) / N_r;

    Threshold_CUT = alpha_set * sigma_hat;
    T_judgement = res_inf_normSq_rot / Threshold_CUT - 1;

    % cfar_detector = phased.CFARDetector2D('TrainingBandSize', training_band - guard_band, 'GuardBandSize', guard_band);
    % cfar_detector.Method = CFAR_method;
    % if strcmp(CFAR_method, 'OS')
        %     cfar_detector.Rank = round(N_r / 2);
    % end
    % cfar_detector.ThresholdFactor = 'Auto';
    % cfar_detector.ProbabilityFalseAlarm = p_fa_CFAR;
    % cfar_detector.ThresholdOutputPort = true;
    % cfar_detector.NoisePowerOutputPort = true;
    % [~, Threshold_CUT, ~] = cfar_detector(prob_ind_ext, [peak_n + Nx; peak_m + My]);


end