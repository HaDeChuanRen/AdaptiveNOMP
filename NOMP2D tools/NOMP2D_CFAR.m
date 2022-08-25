function [omegaList, gainList, y_residue_matrix, Threshold_collect] =...
NOMP2D_CFAR(y_matrix, alpha_set, N_r, K_max, guard_band, overSamplingRate, R_s, R_c)
%NOMP2D_CFAR 此处显示有关此函数的摘要
%   此处显示详细说明
    % if ~exist('CFAR_method', 'var'), CFAR_method = 'CA';
    % elseif isempty(CFAR_method), CFAR_method = 'CA'; end

    if ~exist('guard_band','var'), guard_band = [3, 3];
    elseif isempty(guard_band), guard_band = [3, 3]; end

    if ~exist('overSamplingRate','var'), overSamplingRate = [4; 4];
    elseif isempty(overSamplingRate), overSamplingRate = [4; 4]; end

    if ~exist('R_s','var'), R_s = 1;
    elseif isempty(R_s),    R_s = 1; end

    if ~exist('R_c','var'), R_c = 3;
    elseif isempty(R_c),    R_c = 3; end

    [Nx, My] = size(y_matrix);
    ant_idx_Nx = (0 : (Nx - 1))' - (Nx - 1) / 2;
    ant_idx_My = (0 : (My - 1))' - (My - 1) / 2;

    K_est = 0;
    omegaList = zeros(K_max, 2);
    gainList = zeros(K_max, 1);
    y_residue_matrix = y_matrix;
    % y_fftabs = abs(fft2(y_residue_matrix)); figure; imagesc(y_fftabs);

    while (K_est < K_max)
        % [omega_k, gain_k, y_residue_matrix] = detectnew_1D_serial(y_residue_matrix, overSamplingRate);
        [omega_k, gain_k, y_residue_matrix] = DetectNew_2D(y_residue_matrix, overSamplingRate);

        for refione_idx = 1 : R_s
            [y_residue_matrix, omega_k, gain_k] = RefineOne_2D(y_residue_matrix, omega_k, gain_k);
        end

        K_est = K_est + 1;
        omegaList(K_est, :) = omega_k;
        gainList(K_est) = gain_k;

        [omega_est, ~, ~] = RefineAll_2D(y_residue_matrix, omegaList(1 : K_est, :), gainList(1 : K_est), R_s, R_c);
        [gain_est, y_residue_matrix, A_all_omega] = LeastSquares_2D(y_matrix, omega_est);
        omegaList(1 : K_est, :) = omega_est;
        gainList(1 : K_est) = gain_est;

        % [omegaList, gainList, y_residue_matrix] = refineall_1D_serial(y_residue_matrix, omegaList, gainList, R_s, R_c);
        % [gainList, y_residue_matrix, A_all_omega] = LeastSquares_1D(y_matrix, omegaList, ant_idx);
    end

    % Threshold_collect = zeros(K_est, 1);
    Khat = length(gainList);

    while true
        y_r_det = zeros(Nx, My, Khat);
        % alpha_hat = zeros(Khat, 1);
        Tarray_judgement = zeros(Khat, 1);
        Threshold_collect = zeros(Khat, 1);
        omegaList_save = zeros(Khat - 1, 2, Khat);
        gainList_save = zeros(Khat - 1, Khat);
        tar_set = 1 : Khat;
        if Khat > 1
            for kidx = Khat : -1 : 1
                % hypothesis testing problem
                % throwing the kidx th target and calculate the threshold

                % y_r_det(:, :, kidx) = y_residue_matrix + gainList(kidx) * reshape(A_all_omega(:, kidx), [Nx, My]);
                y_r_det_sq = squeeze(y_r_det(:, :, kidx));
                tar_set_diff = setdiff(tar_set, kidx);

                [omegaList_temp, ~, ~] = RefineAll_2D(y_r_det_sq, omegaList(tar_set_diff, :), gainList(tar_set_diff,:), R_s, R_c);
                [gainList_temp, y_test, ~] = LeastSquares_2D(y_matrix, omegaList_temp);
                y_r_det(:, :, kidx) = y_test;

                [T_judgement, Threshold_CUT] = CFAR_detector2D(y_test, N_r, ...
                alpha_set, overSamplingRate, omegaList, guard_band);
                % alpha_hat(kidx) = res_inf_normSq_rot / sigma_hat;
                omegaList_save(:, :, kidx) = omegaList_temp;
                gainList_save(:, kidx) = gainList_temp;
                Tarray_judgement(kidx) = T_judgement;
                Threshold_collect(kidx) = Threshold_CUT;
            end
        else
            ymat_test = y_residue_matrix + reshape(A_all_omega, Nx, My) * gainList;
            [T_judgement, Threshold_CUT] = CFAR_detector2D(ymat_test, N_r, ...
            alpha_set, overSamplingRate, omegaList, guard_band);
            % alpha_hat0 = res_inf_normSq_rot / sigma_hat;
            % Tarray_judgement = alpha_hat0/tau-1;
            Tarray_judgement = T_judgement;
            Threshold_collect = Threshold_CUT;
        end
        %Tarray_judgement(k): reflects the existance of the kth target: >0, existence,
        %else deactive this target
        [~, ktar_idx] = min(Tarray_judgement);
        if Tarray_judgement(ktar_idx) < 0
            omegaList = omegaList_save(:, :, ktar_idx);
            gainList = gainList_save(:, ktar_idx);
            % omegaList(ktar_idx, :) = [];
            % gainList(ktar_idx, :) = [];
            Khat = Khat - 1;
            if Khat == 0
                omegaList = [];
                gainList = [];
                Threshold_collect = [];
                break;
            end
            y_residue_matrix = squeeze(y_r_det(:, :, ktar_idx));
            [omegaList, ~, ~] = RefineAll_2D(y_residue_matrix, omegaList, gainList, R_s, R_c);
            [gainList, y_residue_matrix, A_all_omega] = LeastSquares_2D(y_matrix, omegaList);
        else
            [T_judgement, ~] = CFAR_detector2D(y_residue_matrix, ...
            N_r, alpha_set, overSamplingRate, omegaList, guard_band);
            % alpha_hat0 = res_inf_normSq_rot / sigma_hat;
            % Tarray_judgement0 = alpha_hat0/tau-1;
            if (T_judgement > 0) && (Khat < K_max)
                % detect
                [omega_new, gain_new, y_residue_matrix] = DetectNew_2D(y_residue_matrix, overSamplingRate);
                for i = 1 : R_s
                    [y_residue_matrix, omega_new, gain_new] = RefineOne_2D(y_residue_matrix, omega_new, gain_new);
                end
                omegaList = [omegaList; omega_new];
                gainList  = [gainList; gain_new];
                % Threshold_collect = [Threshold_collect; Threshold_CUT];
                Khat = Khat + 1;

                % refine all frequencies detected so far
                % can be interpreted as a search for better frequency supports
                [omegaList, gainList, y_residue_matrix] = RefineAll_2D(y_residue_matrix, omegaList, gainList, R_s, R_c);
                % refineAll only uses refineOne to tweak parameters and the energy
                % in the residual measurements y_r can only decrease as a result

                % Solve least squares for the dictionary set [Ax(omega)] omega in
                % omegaList
                [gainList, y_residue_matrix, A_all_omega] = LeastSquares_2D(y_matrix, omegaList);
            else
                break;
            end
        end
    end

    if ~isempty(omegaList)
        gainList = bsxfun(@times, gainList, exp(1j * (ant_idx_Nx(1) * omegaList(:, 1) + ant_idx_My(1) * omegaList(:, 2))));

        % gainList = gainList .* exp(1j*sampledManifold.ant_idx(1)*omegaList);
        omegaList =  wrapTo2Pi(omegaList);
    end


end

