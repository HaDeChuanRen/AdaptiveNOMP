function [omega, gain, y_r, res_inf_normSq_rot] = detectNew(y, sampledManifold)
% last update date: 2022.6.30

% SUMMARY:
% 
% 	detects a new sinusoid on the coarse grid
% 
% INPUT:
% 	y - measurements
% 	sampledManifold - when measurements are compressive,
% 		  we precompute (using IFFT operation) a
% 		  **dictionary** of responses for sinusoids 
%  		  corresponding to frequencies on a coarse grid 
% 		  and store them in the MATLAB **structure** 
% 		  sampledManifold
%
% OUTPUT:
% 	omega - frequency on [0,2*pi) which best explains
% 			the measurements y
% 	gain  - corresponding complex gain
% 	y_r   - after removing the detected sinusoid from the
%         	measurements y, y_r is the ***residual measurement***
%          res_inf_normSq_rot - max energy among DFT directions - needed
%          for stopping criterion

R = length(sampledManifold.coarseOmega);
N = sampledManifold.length;
OSR = round(R/N);

if sampledManifold.is_eye
    gains  = fft(y, R)/sqrt(N);   % R *T matrix
    
    if sampledManifold.ant_idx(1)~=0
        gains = bsxfun(@times,gains,exp(-1j*sampledManifold.coarseOmega(:)*...
            sampledManifold.ant_idx(1)));
%         gains = gains.*exp(-1j*sampledManifold.coarseOmega(:)*...
%             sampledManifold.ant_idx(1));
    end
    prob = sum(abs(gains).^2,2);
else
    energy = sampledManifold.map_IfftMat_norm_sq.';   % R*1
    gains = bsxfun(@rdivide,sampledManifold.map_IfftMat'*y,energy); % R*T
%     gains = (sampledManifold.map_IfftMat'*y)./energy;
    prob = sum(bsxfun(@times,abs(gains).^2, energy),2);
%     prob = (abs(gains).^2).*energy;
end

[~,IDX] = max(prob);

omega = sampledManifold.coarseOmega(IDX);
gain = gains(IDX,:);

% compute the response corresponding to the
% current estimate of the sinusoid
if sampledManifold.is_eye
    x = exp(1j*sampledManifold.ant_idx * omega)...
        /sqrt(sampledManifold.length);
else
    x = sampledManifold.map_IfftMat(:,IDX);
end

% residual measurements after subtracting
% out the newly detected sinusoid
y_r = y - x*gain;

% For stopping criterion
% we check only DFT frequencies - gives us a handle
% the false alarm rate (a measure of how often we would 
% over estimate the size of the support set)

res_inf_normSq_rot = max(prob(1:OSR:end));

end