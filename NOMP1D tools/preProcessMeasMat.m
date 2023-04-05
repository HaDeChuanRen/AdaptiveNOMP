function sampledManifold = preProcessMeasMat(S, overSamplingRate)
% SUMMARY:
%   compute overSamplingRate*M IFFTs once when measurement matrix
%   is not equal to eye(N); avoids repetitive IFFT computations 
%   later on
% INPUT:
%   S - M times N measurement matrix
%       M - number of measurements
%       N - length of sinusoid
%   examples: eye(N) - normal measurement matrix
%             diag(hamming(N)) - normal measurements with hamming
%                                weights
%             randn(M,N)/sqrt(N) - compressive measurement matrix
%             randn(M,N)/sqrt(N) * diag(hamming(N)) - compressive
%                       measurements with hamming weights matrix
%    overSamplingRate (optional) - how fine (in terms of multiples
%       of the FFT grid) do we want the coarse grid to be?
%       number of grid points = overSamplingRate * N
%       Default value is 3
%
% OUTPUT:
%       data structure with sinusoid responses 
%       precomputed using IFFTs when S is not equal to eye(N)


M = size(S,1);
N = size(S,2);

sampledManifold.length = N;
R = round(overSamplingRate * N);

sampledManifold.coarseOmega = 2*pi*(0:(R-1))/R;  % omegaCoarse

% definition of a sinusoid: 
%              exp(-1j*omega((N-1)/2):((N-1)/2))/sqrt(N)

ant_idx = 0:(N-1);
ant_idx = ant_idx - (N-1)/2;

% IFFT definition of a sinusoid(omega) takes the following form:
% 	sinusoid    = @(omega) exp(1j*(0:(N-1)).'*omega);
% To reiterate, we assume that a sinusoid is given by
%	sinusoid    = @(omega) exp(1j*ant_idx.'*omega)/sqrt(N);
% So we store this information in sampledManifold container
sampledManifold.ant_idx = ant_idx(:);

% Check if the measurement matrix is the identity matrix
if M == N
    is_eye = norm(eye(N) - S, 'fro') == 0;
else
    is_eye = false;
end
sampledManifold.is_eye = is_eye;

% WHEN THE MEASUREMENT MATRIX NOT THE IDENTITY MATRIX
% we compute the Fourier transforms of sensing weights 
% and store them for use in the coarse stage (can use 
% IFFT to speed up these steps)

% IN THE FOLLOWING COMPUTATIONS WE COMPENSATE FOR THE 
% DIFFERENCES BETWEEN NOTION OF A SINUSOID THAT WE TAKE 
% AND THAT TAKEN BY "IFFT"

if ~is_eye
    
    % S times x(omegaCoarse)
    sampledManifold.map_IfftMat = R/sqrt(N)*ifft(S,R,2) * ...
        sparse(1:R, 1:R, exp(1j * sampledManifold.coarseOmega...
        * ant_idx(1)), R, R);
    
    % norm square of S times x(omegaCoarse)
    % energy as a function of frequency
    sampledManifold.map_IfftMat_norm_sq = ...
        sum(abs(sampledManifold.map_IfftMat).^2, 1);
    
end
end
