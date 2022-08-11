% last update: 2022/7/18

clc; clear; close all;

N = 256;
N_r = 50;

alpha_all = (1 : 0.01 : 15)';
S_all = [10, 50, 100];
% S_max = 2000;
T_var = 1 : 5000;
% S_all = 1;
length_alpha = length(alpha_all);
length_S = length(S_all);
Poe_alpha = zeros(length_alpha, length_S);

Poe_NOMPS1 = chi2cdf(2 * alpha_all, 2) .^ N;

% (1 - exp(- alpha_all))


lw = 1.6;
fsz = 12;
msz = 10;

Poe_alphaS1 = zeros(length_alpha, 1);
for alpha_idx = 1 : length_alpha
    S_snap = 1;

    alpha_set = alpha_all(alpha_idx);
    fun_Poe_alpha = @(Tvar) (N_r / alpha_set) * exp(N * log(1 - exp(- Tvar)) + ...
    (S_snap * N_r - 1) * log(N_r * Tvar / alpha_set) - N_r * Tvar / alpha_set - sum(log(1 : (S_snap * N_r - 1))));
    Poe_alphaS1(alpha_idx) = sum(fun_Poe_alpha(T_var));
    % Poe_alphaS1(alpha_idx) = integral(fun_Poe_alpha, 0, Inf);
end

% .* (Tvar .^ k_seq ./ factorial(k_seq))
figure(1)
plot(alpha_all, Poe_alphaS1, 'r', 'Linewidth', lw, 'Markersize', msz);
hold on;
plot(alpha_all, Poe_NOMPS1, 'r--', 'Linewidth', lw, 'Markersize', msz);



for S_idx = 1 : length(S_all)
    S_snap = S_all(S_idx);
    k_seq = (0 : S_snap - 1)';
    ln_factorialS = zeros(S_snap, 1);
    for k_idx = 1 : S_snap - 1
        ln_factorialS(k_idx + 1) = sum(log(1 : k_idx));
    end
    % Poe_old = 0;
    % k_seq = (S_snap : S_max - 1)';
    % ln_factorialS = zeros(S_max - S_snap, 1);
    % for k_idx = S_snap : (S_max - 1)
    %     ln_factorialS(k_idx - S_snap + 1) = sum(log(1 : k_idx));
    % end
    for alpha_idx = 1 : length_alpha
        alpha_set = alpha_all(alpha_idx);
        fun_Poe_alpha = @(Tvar) (N_r / alpha_set) *...
        exp(N * log(1 - sum(exp(- Tvar + k_seq .* log(Tvar) - ln_factorialS))) + ...
        (S_snap * N_r - 1) * log(N_r * Tvar / alpha_set) - ...
        N_r * Tvar / alpha_set - sum(log(1 : (S_snap * N_r - 1))));

        Poe_temp = sum(fun_Poe_alpha(T_var));
        % Poe_temp =  integral(fun_Poe_alpha, 0, Inf);
        Poe_alpha(alpha_idx, S_idx) = Poe_temp;
%         if Poe_temp > Poe_old
%             Poe_alpha(alpha_idx, S_idx) = Poe_temp;
%             Poe_old = Poe_temp;
%         else
%             Poe_alpha(alpha_idx, S_idx) = Poe_old;
%         end
    end
end

Poe_NOMPSlow = (chi2cdf(2 * alpha_all * S_all(1), 2 * S_all(1))) .^ N;
Poe_NOMPSmedium =(chi2cdf(2 * alpha_all * S_all(2), 2 * S_all(2))) .^ N;
Poe_NOMPShigh =(chi2cdf(2 * alpha_all * S_all(3), 2 * S_all(3))) .^ N;
Poe_all = (0.001 : 0.001 : 0.999)';
% alpha_NOMPlow = chi2inv((1 - Poe_all) .^ (1 / N), 2 * S_all(1)) / 2;
% alpha_NOMPmedium = chi2inv((1 - Poe_all) .^ (1 / N), 2 * S_all(2)) / 2;

plot(alpha_all, Poe_alpha(:, 1), 'b', 'Linewidth', lw, 'Markersize', msz);
plot(alpha_all, Poe_NOMPSlow, 'b--', 'Linewidth', lw, 'Markersize', msz);
plot(alpha_all, Poe_alpha(:, 2), 'g', 'Linewidth', lw, 'Markersize', msz);
plot(alpha_all, Poe_NOMPSmedium, 'g--', 'Linewidth', lw, 'Markersize', msz);
plot(alpha_all, Poe_alpha(:, 3), 'c', 'Linewidth', lw, 'Markersize', msz);
plot(alpha_all, Poe_NOMPShigh, 'c--', 'Linewidth', lw, 'Markersize', msz);

legend('S = 1 (CA)', 'S = 1 (NOMP)', 'S = 10 (CA)', 'S = 10 (NOMP)', ...
'S = 50 (CA)', 'S = 50 (NOMP)', 'S = 100 (CA)', 'S = 100 (NOMP)');
xlabel('$\alpha$', 'Interpreter', 'latex', 'Fontsize', fsz)
ylabel('$1 - \bar{P}_{\rm FA}$', 'Interpreter', 'latex', 'Fontsize', fsz)
ylim([0 1])



% length_Poe = length(Poe_all);
% alpha_cal = zeros(length_Poe, 1);
% S = 50;
% for p_idx = 1 : length_Poe
%     Poe = Poe_all(p_idx);
%     alpha_cal(p_idx) = alpha_PoebyS(Poe, N, N_r, S);
% end
% 
% figure(2);
% plot(alpha_cal, Poe_all, 'Linewidth', lw, 'Markersize', msz)
% hold on;
% plot(alpha_all, 1- Poe_alpha(:, 2), 'Linewidth', lw, 'Markersize', msz);







% S = 50;
% k_seq = (1 : S)';
% T_var = 1 : 1000;
% alpha_set = 12;
% ln_Poe_Tvar = N * log(1 - exp(- T_var) .* (T_var .^ k_seq ./ factorial(k_seq))) + ...
% (S * N_r - 1) * log(N_r * T_var / alpha_set) - N_r * T_var / alpha_set - sum(log(1 : (S * N_r - 1)));
%
% plot(T_var, ln_Poe_Tvar);
% sum(fun_Poe(T_var))










