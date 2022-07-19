% last update: 2022/7/18

clc; clear; close all;

N = 256;
N_r = 60;

alpha_all = 1 : 0.01 : 20;
S_all = [5, 10, 15];
% S_all = 1;
length_alpha = length(alpha_all);
length_S = length(S_all);
Poe_alpha = zeros(length_alpha, length_S);

lw = 1.6;
fsz = 12;
msz = 10;

Poe_alphaS1 = zeros(length_alpha, 1);
for alpha_idx = 1 : length_alpha
    S = 1;
    k_seq = (1 : S)';
    alpha_set = alpha_all(alpha_idx);
    fun_Poe_alpha = @(Tvar) (N_r / alpha_set) * exp(N * log(1 - exp(- Tvar)) + ...
    (S * N_r - 1) * log(N_r * Tvar / alpha_set) - N_r * Tvar / alpha_set - sum(log(1 : (S * N_r - 1))));
    Poe_alphaS1(alpha_idx) = 1 - integral(fun_Poe_alpha, 0, Inf);
end

% .* (Tvar .^ k_seq ./ factorial(k_seq))
figure(1)
plot(alpha_all, Poe_alphaS1, 'Linewidth', lw, 'Markersize', msz);
hold on;

for S_idx = 1 : length(S_all)
    S = S_all(S_idx);
    k_seq = (1 : S)';
    for alpha_idx = 1 : length_alpha
        alpha_set = alpha_all(alpha_idx);
        fun_Poe_alpha = @(Tvar) (N_r / alpha_set) * exp(N * log(1 - exp(- Tvar) .* sum(Tvar .^ k_seq ./ factorial(k_seq))) + ...
        (S * N_r - 1) * log(N_r * Tvar / alpha_set) - N_r * Tvar / alpha_set - sum(log(1 : (S * N_r - 1))));
        Poe_alpha(alpha_idx, S_idx) = 1 - integral(fun_Poe_alpha, 0, Inf);
    end
end


plot(alpha_all, Poe_alpha, 'Linewidth', lw, 'Markersize', msz);
legend('S = 1', 'S = 5', 'S = 10', 'S = 15');
xlabel('$\alpha$', 'Interpreter', 'latex', 'Fontsize', fsz)
ylabel('$\bar{P}_{\rm oe}$', 'Interpreter', 'latex', 'Fontsize', fsz)
ylim([0 1])


Poe_all = 0.01 : 0.01 : 0.99;
length_Poe = length(Poe_all);
alpha_cal = zeros(length_Poe, 1);
S = 10;
for p_idx = 1 : length_Poe
    Poe = Poe_all(p_idx);
    alpha_cal(p_idx) = alpha_PoebyS(Poe, N, N_r, S);
end

figure(2);
plot(alpha_cal, Poe_all, 'Linewidth', lw, 'Markersize', msz)
hold on;
plot(alpha_all, Poe_alpha(:, 2), 'Linewidth', lw, 'Markersize', msz);


% S = 1;
% k_seq = (1 : S)';
% T_var = 1 : 1000;
% alpha_set = 12;
% ln_Poe_Tvar = N * log(1 - exp(- T_var) .* (T_var .^ k_seq ./ factorial(k_seq))) + ...
% (S * N_r - 1) * log(N_r * T_var / alpha_set) - N_r * T_var / alpha_set - sum(log(1 : (S * N_r - 1)));
%
% plot(T_var, ln_Poe_Tvar);
% sum(fun_Poe(T_var))










