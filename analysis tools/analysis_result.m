function results_struct = analysis_result(omega_true, gain_true, omega_hat, gainList, N, gamma_oversamping)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
    Khat = length(omega_hat);
    K = length(omega_true);
    Khat_on = 0;
    Khat_off = 0;
    % omega_min = 2 * pi / N;
    omega_trueidx = ceil(omega_true * N * gamma_oversamping / (2 * pi)) + 1;
    results_struct.Miss_Eve = 0;
    results_struct.error_ave = 0;
    results_struct.Detect_Eve = 0;
    results_struct.Equal_Eve = 0;
    results_struct.Overest_Eve = 0;
    results_struct.reconErr = 0;

    Dist_omega = zeros(K, Khat);
    if K > 0
        y_full = exp(1j * (0 : (N - 1)).' * omega_true.') / sqrt(N) * gain_true;
    else
        y_full = zeros(N, 1);
    end

    % calculate the reconstruct error
    if Khat > 0
        for k_idx = 1 : K
            for k_jdx = 1 : Khat
                Dist_omega(k_idx, k_jdx) = sum(abs(wrapToPi(omega_true(k_idx) - omega_hat(k_jdx))).^2);
            end
        end
        [i_Hun, j_Hun] = Hungarianmatch(Dist_omega);
        % Dist_judge = Dist_omega((j_Hun - 1) * K + i_Hun) > omega_min ^ 2;
        y_recon = exp(1j * (0 : (N - 1)).' * omega_hat.') / sqrt(N) * gainList;
        results_struct.reconErr = norm(y_full - y_recon) ^ 2 / norm(y_full)^2;
    else
        results_struct.reconErr = 1;
    end

    % calculate the detection event and false alarm event
    omega_hatidx = ceil(omega_hat * N * gamma_oversamping / (2 * pi)) + 1;
    vec_hatidx = zeros(N * gamma_oversamping + 2, 1);
    vec_hatidx(omega_hatidx) = 1;
    for k_kdx = 1 : K
        if sum(vec_hatidx(omega_trueidx(k_kdx) + (-1 : 1))) >= 1
            Khat_on = Khat_on + 1;
            Khat_off = Khat_off + sum(vec_hatidx(omega_trueidx(k_kdx) + (-1 : 1))) - 1;
            % vec_hatidx(omega_trueidx(k_kdx) + (-1 : 1)) = 0;
        end
    end
    results_struct.False_Eve = Khat - Khat_on;
    if Khat_on == K
        results_struct.Detect_Eve = 1;
    end
    % sum(vec_hatidx) + Khat_off;
    if Khat > K
        results_struct.Overest_Eve = 1;
        % iinj_judge = ~ismember(j_Hun, i_Hun);
        % if sum(iinj_judge) + sum(Dist_judge) == 0
        %     results_struct.Detect_Eve = 1;
        % end
    elseif Khat == K
    % calculate the estimation error
        % results_struct.Detect_Eve = 1;
        results_struct.Equal_Eve = 1;
        error_omega = 0;
        for k_ldx = 1 : K
            error_omega = error_omega + Dist_omega(i_Hun(k_ldx), j_Hun(k_ldx));
        end
        results_struct.error_ave = error_omega / K;
    else
        results_struct.Miss_Eve = K - Khat + Khat_off;
    end
end

