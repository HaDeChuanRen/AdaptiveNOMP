function results_struct = False_Detection(omega_true, gain_true, ...
    omega_hat, gainList, N)
%   this function is used to ananlysis the detection and estimation result
%   此处显示详细说明
    Khat = length(omega_hat);
    K = length(omega_true);
    omega_min = 2 * pi / N;

    results_struct.Miss_Eve = 0;
    results_struct.error_ave = 0;
    results_struct.False_Eve = 0;
    results_struct.Detect_Eve = 0;
    results_struct.Equal_Eve = 0;
    results_struct.Overest_Eve = 0;
    results_struct.reconErr = 0;

    if K > 0
        y_full = exp(1j * (0 : (N - 1)).' * omega_true.') / sqrt(N) * gain_true;
    else
        y_full = zeros(N, 1);
    end

    % calculate the reconstruct error
    if Khat > 0
        Dist_omega = abs(wrapToPi(omega_hat' - omega_true)) .^ 2;
        y_recon = exp(1j * (0 : (N - 1)).' * omega_hat.') / sqrt(N) * gainList;
        results_struct.reconErr = norm(y_full - y_recon) ^ 2 / norm(y_full) ^ 2;
    else
        results_struct.reconErr = 1;
        Dist_omega = [];
    end

    % judge the detection event
    Detect_vec = min(Dist_omega, [], 2) < (0.5 * omega_min) ^ 2;
    if sum(Detect_vec) == K
        results_struct.Detect_Eve = 1;
    end
    results_struct.Miss_Eve = K - sum(Detect_vec);

    % judge the false alarm event
    False_vec = min(Dist_omega) > (0.5 * omega_min) ^ 2;
    results_struct.False_Eve = sum(False_vec);
%         results_struct.False_Eve = 1;
%     end

    if Khat > K
        results_struct.Overest_Eve = 1;
    elseif Khat == K
        results_struct.Equal_Eve = 1;
        error_omega = 0;
        [i_Hun, j_Hun] = Hungarianmatch(Dist_omega);
        for k_ldx = 1 : K
            error_omega = error_omega + Dist_omega(i_Hun(k_ldx), j_Hun(k_ldx));
        end
        results_struct.error_ave = error_omega / K;
    end
end

