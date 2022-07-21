clc; clear; close all;

N = 256;
N_r = 60;

S_snap = 50;

Tvar = 1 : 1000;
alpha_set = 6;

k_seq = (1 : S_snap)';
ln_factorialS = zeros(S_snap, 1);
for k_idx = 2 : S_snap
    ln_factorialS(k_idx + 1) = sum(log(1 : k_idx));
end

fun_Poe_alpha = @(Tvar) (N_r / alpha_set) *...
exp(N * log(1 - sum(exp(- Tvar + k_seq .* log(Tvar) - ln_factorialS))) + ...
(S_snap * N_r - 1) * log(N_r * Tvar / alpha_set) - ...
N_r * Tvar / alpha_set - sum(log(1 : (S_snap * N_r - 1))));


% S_max = 2000;
% k_seq = (S_snap : S_max - 1)';
% ln_factorialS = zeros(S_max - S_snap, 1);
% for k_idx = S_snap : (S_max - 1)
%     ln_factorialS(k_idx - S_snap + 1) = sum(log(1 : k_idx));
% end

% fun_Poe_alpha = @(Tvar) (N_r / alpha_set) *...
% exp(N * log(sum(exp(- Tvar + k_seq .* log(Tvar) - ln_factorialS))) + ...
% (S_snap * N_r - 1) * log(N_r * Tvar / alpha_set) - ...
% N_r * Tvar / alpha_set - sum(log(1 : (S_snap * N_r - 1))));

% S_max = 2000;
% k_seq = (S_snap : (S_max - 1))';
% ln_factorialS = zeros(S_max - S_snap, 1);
% for k_idx = (S_snap + 1) : (S_max - 1)
%     ln_factorialS(k_idx - S_snap + 1) = sum(log(S_snap + 1 : k_idx));
% end
% 
% ln_Poe_Tvar = N * (- Tvar + S_snap * log(Tvar) - sum(log(1 : S_snap))) +...
% N * log(sum(exp((k_seq - S_snap) .* log(Tvar) - ln_factorialS))) + ...
% (S_snap * N_r - 1) * log(N_r * Tvar / alpha_set) - ...
% N_r * Tvar / alpha_set - sum(log(1 : (S_snap * N_r - 1)));


% fun_Poe_alpha = @(Tvar) (N_r / alpha_set) *...
% exp(N * (- Tvar + S_snap * log(Tvar) - sum(log(1 : S_snap))) +...
% N * log(sum(exp((k_seq - S_snap) .* log(Tvar) - ln_factorialS))) + ...
% (S_snap * N_r - 1) * log(N_r * Tvar / alpha_set) - ...
% N_r * Tvar / alpha_set - sum(log(1 : (S_snap * N_r - 1))));
% 
% 
% % exp(- Tvar) .* sum(Tvar .^ k_seq ./ factorial(k_seq))
% 
% Poe_alpha1 = integral(fun_Poe_alpha, 0, 1e4)
% % Poe_alpha2 = integral(fun_Poe_alpha, 0, 1e5)
% Poe_alpha3 = integral(fun_Poe_alpha, 0, inf)


% figure;
% plot(T_var, ln_Poe_Tvar);
% sum((N_r / alpha_set) * exp(ln_Poe_Tvar))

% k_seq = (1 : S_snap)';
% ln_Poe_Tvar = N * log(1 - exp(- Tvar) .* sum(Tvar .^ k_seq ./ factorial(k_seq))) + ...
% (S_snap * N_r - 1) * log(N_r * Tvar / alpha_set) - N_r * Tvar / alpha_set -...
% sum(log(1 : (S_snap * N_r - 1)));
