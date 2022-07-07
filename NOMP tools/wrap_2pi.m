function omega_prime= wrap_2pi(omega)
% SUMMARY: Restricts frequencies to [0, 2*pi)
% INPUT: A vector of spatial frequencies
% OUTPUT: A vector of coresponding spatial frequencies in [0, 2*pi)

omega_prime = angle(exp(1j*omega));
omega_prime(omega_prime < 0) = omega_prime(omega_prime < 0) + 2*pi;

end