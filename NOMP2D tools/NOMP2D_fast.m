function [omegaList, gainList, Ymat_res, Threshold_collect] = ...
    NOMP2D_fast(Ymat, alpha_set, N_r, K_max, guard_band, ...
    overSamplingRate, R_s, R_c)

    if ~exist('guard_band','var'), guard_band = [3, 3];
    elseif isempty(guard_band), guard_band = [3, 3]; end

    if ~exist('overSamplingRate','var'), overSamplingRate = [4; 4];
    elseif isempty(overSamplingRate), overSamplingRate = [4; 4]; end

    if ~exist('R_s','var'), R_s = 1;
    elseif isempty(R_s),    R_s = 1; end

    if ~exist('R_c','var'), R_c = 3;
    elseif isempty(R_c),    R_c = 3; end

    [Nx, My] = size(Ymat);
    antvec_Nx = (0 : (Nx - 1))' - (Nx - 1) / 2;
    antvec_My = (0 : (My - 1))' - (My - 1) / 2;

    K_est = 0;
    omegaList = [];
    gainList = [];
    Threshold_collect = [];
    Ymat_res = Ymat;
    % y_fftabs = abs(fft2(y_residue_matrix)); figure; imagesc(y_fftabs);
    flag_cfar = 0;

    while (K_est < K_max) && (flag_cfar < 3)
        [T_judgement, Threshold_CUT] = CFAR_detector2D(Ymat_res, N_r, ...
        alpha_set, overSamplingRate, omegaList, guard_band);

        [omega_k, gain_k, Ymat_res] = DetectNew_2D(Ymat_res, ...
            overSamplingRate);
        for refione_idx = 1 : R_s
            [Ymat_res, omega_k, gain_k] = RefineOne_2D(Ymat_res, ...
            omega_k, gain_k);
        end

        if T_judgement < 0
            flag_cfar = flag_cfar + 1;
            continue;
        else
            flag_cfar = 0;
            K_est = K_est + 1;
            omegaList = [omegaList; omega_k];
            gainList = [gainList; gain_k];
            Threshold_collect = [Threshold_collect; Threshold_CUT];
        end
    end

    if ~isempty(omegaList)
        gainList = bsxfun(@times, gainList, exp(1j * (antvec_Nx(1) * ...
        omegaList(:, 1) + antvec_My(1) * omegaList(:, 2))));

        % gainList = gainList .* exp(1j*sampledManifold.ant_idx(1)*omegaList);
        omegaList =  wrapTo2Pi(omegaList);
    end

end










