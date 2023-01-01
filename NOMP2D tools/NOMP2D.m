function [omegaList, gainList, y_residue_matrix] =...
NOMP2D(y_matrix, tau, K_max, overSamplingRate, R_s, R_c)
%NOMP2D_CFAR 此处显示有关此函数的摘要
%   此处显示详细说明


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
    omegaList = [];
    gainList = [];
    y_residue_matrix = y_matrix;
    % y_fftabs = abs(fft2(y_residue_matrix)); figure; imagesc(y_fftabs);

    while (1)
        % [omega_k, gain_k, y_residue_matrix] = detectnew_1D_serial(y_residue_matrix, overSamplingRate);
        [omega_k, gain_k, y_residue_matrix] = DetectNew_2D(y_residue_matrix, overSamplingRate);

        if abs(gain_k) ^ 2 < tau
            break;
        end

        for refione_idx = 1 : R_s
            [y_residue_matrix, omega_k, gain_k] = RefineOne_2D(y_residue_matrix, omega_k, gain_k);
        end

        K_est = K_est + 1;
        omegaList = [omegaList; omega_k];
        gainList = [gainList; gain_k];

        [omega_est, ~, ~] = RefineAll_2D(y_residue_matrix, omegaList, gainList, R_s, R_c);
        [gain_est, y_residue_matrix, ~] = LeastSquares_2D(y_matrix, omega_est);
        omegaList(1 : K_est, :) = omega_est;
        gainList(1 : K_est) = gain_est;

        % [omegaList, gainList, y_residue_matrix] = refineall_1D_serial(y_residue_matrix, omegaList, gainList, R_s, R_c);
        % [gainList, y_residue_matrix, A_all_omega] = LeastSquares_1D(y_matrix, omegaList, ant_idx);
    end

    if ~isempty(omegaList)
        gainList = bsxfun(@times, gainList, exp(1j * (ant_idx_Nx(1) * omegaList(:, 1) + ant_idx_My(1) * omegaList(:, 2))));
        % gainList = gainList .* exp(1j*sampledManifold.ant_idx(1)*omegaList);
        omegaList =  wrapTo2Pi(omegaList);
    end
end

