function [omega_est, ghat, y_residue_matrix] = RefineAll_2D(y_residue_matrix, omega_est, ghat, R_s, R_c)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
    [K_hat, ~] = size(omega_est);          % the number of detected frequencies
    % order = 1 : K_hat;
    for c_idx = 1 : R_c
        for k_idx = 1 : K_hat
            % l = order(j);
            omega_est_k = omega_est(k_idx, :);
            gest_kidx = ghat(k_idx);

            for s_idx = 1 : R_s
                [y_residue_matrix, omega_hat_kidx, ghat_kidx] = RefineOne_2D(y_residue_matrix, omega_est_k, gest_kidx);
            end
            omega_est(k_idx, :) = omega_hat_kidx;
            ghat(k_idx) = ghat_kidx;
        end
    end
end

