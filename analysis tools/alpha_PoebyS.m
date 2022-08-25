function alpha_hat = alpha_PoebyS(P_oe, N, N_r, S_snap)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明

    if ~exist('S_snap','var'), S_snap = 1;
    elseif isempty(S_snap),    S_snap = 1; end

    alpha_min = 1;
    alpha_max = 100;
    accuracy_result = 1e-12;

    if S_snap == 1
        while abs(alpha_max - alpha_min) > accuracy_result
            alpha_set = (alpha_max + alpha_min) / 2;
            fun_Poe = @(xvar) exp(N * log(1 - exp(- alpha_set / (2 * N_r) * xvar))+...
            (N_r - 1) * log(xvar) - xvar / 2 - N_r * log(2) - sum(log(1 : N_r - 1)));

            % fun_Poe_alpha = @(Tvar) (N_r / alpha_set) * exp(N * log(1 - exp(- Tvar)) + ...
            % (S_snap * N_r - 1) * log(N_r * Tvar / alpha_set) - N_r * Tvar / alpha_set - ...
            % sum(log(1 : (S_snap * N_r - 1))));

            int_result = integral(fun_Poe, 0, Inf);

            delta_result = 1 - int_result - P_oe;
            if delta_result > 0
                alpha_min = alpha_set;
            else
                alpha_max = alpha_set;
            end
        end
    else
        while abs(alpha_max - alpha_min) > accuracy_result
            alpha_set = (alpha_max + alpha_min) / 2;
            k_seq = (0 : S_snap - 1)';
            ln_factorialS = zeros(S_snap, 1);
            for k_idx = 1 : S_snap - 1
                ln_factorialS(k_idx + 1) = sum(log(1 : k_idx));
            end
            fun_Poe = @(Tvar) (N_r / alpha_set) *...
            exp(N * log(1 - sum(exp(- Tvar + k_seq .* log(Tvar) - ln_factorialS))) + ...
            (S_snap * N_r - 1) * log(N_r * Tvar / alpha_set) - ...
            N_r * Tvar / alpha_set - sum(log(1 : (S_snap * N_r - 1))));
            int_result = integral(fun_Poe, 0, 5e3);
            delta_result = 1 - int_result - P_oe;
            if delta_result > 0
                alpha_min = alpha_set;
            else
                alpha_max = alpha_set;
            end
        end
    end

    alpha_hat = alpha_set;

end

