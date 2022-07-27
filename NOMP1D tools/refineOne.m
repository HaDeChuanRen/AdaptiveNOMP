function [omega, gain, y_r] = refineOne(y_r, omega, gain, S,...
			 ant_idx, isOrth)
% SUMMARY:
%   Refines parameters (gain and frequency) of a single sinusoid
%   and updates the residual measurement vector y_r to reflect
%   the refinement -- This function applies one Newton step only.
% INPUT:
% 	y_r - residual measurement (all detected sinusoids removed)
%	omega - current estimate of frequency of sinusoid we want to
% 			refine
%	gain - current estimate of gain of sinusoid we want to refine
% 	S - measurement matrix
%           if [], measurements are direct i.e S = eye(N),
% 	    where N is length of sinusoid	
% 	ant_idx - translating indexes to phases in definition of sinusoid
%	isOrth - binary flag - is y_r orthogonal to x(omega) 
%	       - default - false
% OUTPUT:
%       refined versions of omega, gain and y_r
%       (see INPUT for definitions)

if ~exist('isOrth', 'var'), isOrth = false; end

if isempty(S)
    is_eye = true;
else
    is_eye = false;
end

N = length(ant_idx);
x_theta  = exp(1j*ant_idx*omega)/sqrt(N);
dx_theta = 1j * ant_idx .* x_theta;
d2x_theta = -(ant_idx.^2) .* x_theta;

if ~is_eye
    x_theta   = S * x_theta;
    dx_theta  = S * dx_theta;
    d2x_theta = S * d2x_theta;
end

% add the current estimate of the sinusoid to residue
y = y_r + x_theta*gain;


% UPDATE GAIN
% recompute gain and residue to ensure that 
% y_r is orthogonal to x_theta - this property
% is lost when we refine other sinusoids
if ~isOrth
    if is_eye
        energy = 1;
    else
        energy = x_theta'*x_theta;
    end
    gain = (x_theta'*y)/energy;
    y_r = y - x_theta * gain;
end

% der1 = -2*real(gain * y_r'*dx_theta);
% der2 = -2*real(gain * y_r'*d2x_theta) +...
%     2*abs(gain)^2*(dx_theta'*dx_theta);
der1 = -2*real(gain*(y_r'*dx_theta));
% der2 = -2*real(gain*(y'*d2x_theta));
der2 = -2*real(gain * (y_r'*d2x_theta)) +...
    2*(gain*gain')*(dx_theta'*dx_theta);
% UPDATE OMEGA
if der2 > 0
    omega_next = omega - der1/der2;
else
    omega_next = omega - sign(der1)*(1/4)*(2*pi/N) * 0.25;
end

% COMPUTE x_theta for omega_next so that we can compute 
% gains_next and y_r_next
x_theta  = exp(1j*ant_idx*omega_next)/sqrt(N);
if is_eye
    energy = 1;
else
    x_theta = S * x_theta;
    energy = (x_theta'*x_theta);
end

% UPDATE GAIN
gain_next = (x_theta'*y)/energy;

% UPDATE RESIDUE
y_r_next = y - x_theta*gain_next;

% check for decrease in residue -  needed as a result of 
% non-convexity of residue (even when the cost surface 
% is locally convex); This is the same as checking whether 
% |<y, x_theta>|^2/<x_theta,x_theta> improves as we ensured 
% that y - gain*x_theta is perp to x_theta by recomputing gain

if (y_r_next(:)'*y_r_next(:)) <= (y_r(:)'*y_r(:))
    % commit only if the residue decreases
    omega = omega_next;
    gain = gain_next;
    y_r = y_r_next;
end

end