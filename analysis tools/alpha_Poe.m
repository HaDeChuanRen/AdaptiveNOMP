function alpha_hat = alpha_Poe(P_oe, N, Nr)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明

    alpha_min = 1;
    alpha_max = 100;
    % delta_result = 1;
    accuracy_result = 1e-12;

    while abs(alpha_max - alpha_min) > accuracy_result
        alpha_set = (alpha_max + alpha_min) / 2;
        fun1 = @(xvar) exp(N * log(1 - exp(- alpha_set / (2 * Nr) * xvar))+...
        (Nr - 1) * log(xvar) - xvar / 2 - Nr * log(2) - sum(log(1 : Nr - 1)));
        % int_result = sum((1 - exp(- alpha_set * x_vector / N_r)) .^ N .* (x_vector) .^ (N_r - 1) .* exp(-x_vector)) *...
        % x_step / factorial(N_r - 1);
        int_result = integral(fun1, 0, Inf);

        delta_result = 1 - int_result - P_oe;
        if delta_result > 0
            alpha_min = alpha_set;
        else
            alpha_max = alpha_set;
        end
    end
    alpha_hat = alpha_set;

end

