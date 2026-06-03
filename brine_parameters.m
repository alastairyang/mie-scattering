function [eps_p, eps_pp] = brine_parameters(T, frequency)
%BRINEPARAMETERS Compute the dielectric parameters for brine
    sigma = brine_conductivity(T);
    eps_brine_debye = brine_permittivity(T,frequency);
    [eps_p, eps_pp] = decompose_permittivity(eps_brine_debye, sigma, frequency);
end


function [eps_prime, eps_double_prime] = decompose_permittivity(eps_debye, sigma, frequency)
    % Decompose complex permittivity into real and imaginary parts
    % by combining Debye relaxation and ionic conductivity contributions.
    %
    % Parameters
    % ----------
    % eps_debye  : complex scalar
    %     Complex relative permittivity from Debye single relaxation model
    %     (already computed, e.g. eps_inf + (eps_s - eps_inf)/(1 + 1j*omega*tau))
    % sigma      : float, S/m
    %     Ionic conductivity of brine
    % frequency  : float, Hz
    %     Linear frequency
    %
    % Returns
    % -------
    % eps_prime       : float, real part of total permittivity (epsilon')
    % eps_double_prime: float, loss part of total permittivity (epsilon'')
    %                   defined as POSITIVE, per convention: eps = eps' - j*eps''

    eps0  = 8.854187817e-12;   % F/m
    omega = 2 * pi * frequency;

    % Combine Debye + conductivity into one complex permittivity
    eps_total = eps_debye - 1j * sigma / (omega * eps0);

    % Extract parts — eps'' is defined positive, hence the minus sign
    eps_prime        =  real(eps_total);
    eps_double_prime = -imag(eps_total);
end
