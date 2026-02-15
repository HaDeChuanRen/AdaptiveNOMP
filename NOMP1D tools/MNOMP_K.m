function [omegaList, gainList, y_r] = MNOMP_K(y, K, S, ...
    overSamplingRate, R_s , R_c)

if ~exist('overSamplingRate','var'), overSamplingRate = 4;
elseif isempty(overSamplingRate), overSamplingRate = 4; end

if ~exist('R_s','var'), R_s = 1;
elseif isempty(R_s),    R_s = 1; end

if ~exist('R_c','var'), R_c = 3;
elseif isempty(R_c),    R_c = 3; end

% Algorithm preprocessing
sampledManifold = preProcessMeasMat(S, overSamplingRate);

if sampledManifold.is_eye
	S = []; % erase identity matrix
end

omegaList = [];
gainList  = [];
y_r = y;

k_idx = 0;

while k_idx < K
    [omega_new, gain_new, y_r, ~] = detectNew(y_r, sampledManifold);

    for i = 1:R_s
        [omega_new, gain_new, y_r] = refineOne(y_r, omega_new, ...
            gain_new, S, sampledManifold.ant_idx, true);
    end

    k_idx = k_idx + 1;
    omegaList = [omegaList; omega_new];
    gainList  = [gainList; gain_new];

    [omegaList, ~, y_r] = refineAll(y_r, omegaList, gainList, S, ...
        sampledManifold.ant_idx, R_s, R_c);

    [omegaList, gainList, y_r] = solveLeastSquares(y, omegaList, ...
        S, sampledManifold.ant_idx);
end

gainList = bsxfun(@times,gainList,exp(1j*sampledManifold.ant_idx(1)*omegaList));

omegaList = wrap_2pi(omegaList);

end