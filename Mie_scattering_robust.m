function [Es, Ea, Ee, Eb] = Mie_scattering_robust(r, f, epsp, epsb)
% Mie scattering efficiency factors for spherical inclusions.
% Uses Lentz continued-fraction for A_1(zeta), then stable downward
% recurrence for A_l, l>=2. Riccati-Bessel W via Ulaby upward recurrence.
%
% Source: Ulaby & Long (2014), stabilised per Bohren & Huffman (1983)
% and Wiscombe (1980).

c         = 3e8;
epsb_real = real(epsb);

np = sqrt(epsp);
nb = sqrt(epsb);
n  = np / nb;

lambda = c / f;
chi    = (2*pi*r / lambda) * sqrt(epsb_real);
zeta   = n * chi;

% Wiscombe stopping criterion
l_max = max(ceil(abs(chi)  + 4*abs(chi)^(1/3)  + 2), ...
            ceil(abs(zeta) + 4*abs(zeta)^(1/3) + 2));
l_max = max(l_max, 20);
l_min = 5;

% =========================================================
% A_1(zeta) via Lentz continued fraction (Bohren & Huffman 1983)
% This is stable for ALL zeta including small |zeta|
%
%   A_1(zeta) = 1/zeta * CF
%   where CF = 1 / (2/zeta - 1/(3/zeta - 1/(4/zeta - ...)))
% =========================================================
A1 = lentz_A1(zeta);

% =========================================================
% Downward recurrence for A_l, l = 2..l_max
% A_{l-1} = l/zeta - 1/(A_l + l/zeta)
% Seed: A_1 from Lentz above
% =========================================================
A_down    = zeros(1, l_max + 1);
A_down(1) = A1;
for l = 1 : l_max - 1
    % This recurs UPWARD from A_1 — wait, we need downward.
    % Store all via downward from a HIGH seed, but correct A_1 via Lentz.
    % Simplest: use Lentz for A_1, then recur UPWARD is UNSTABLE.
    % So: use downward from l_max, but REPLACE A_down(1) with Lentz value.
end

% Correct approach: downward from l_max, then override A_down(1) with Lentz
A_down        = zeros(1, l_max + 2);
A_down(l_max+2) = 0;
for l = l_max+1 : -1 : 1
    A_down(l) = (l)/zeta - 1.0 / (A_down(l+1) + (l)/zeta);
end
% Override l=1 with the accurate Lentz value
A_down(1) = A1;

% =========================================================
% Pre-compute Riccati-Bessel W_l via Ulaby upward recurrence
% W_l = (2l-1)/chi * W_{l-1} - W_{l-2}
% W_0  = sin(chi) + i*cos(chi)
% W_{-1} = cos(chi) - i*sin(chi)
% =========================================================
W = zeros(1, l_max + 2);   % W(l+1) = W_l
W(1) = sin(chi) + 1i*cos(chi);    % W_0
Wm1  = cos(chi) - 1i*sin(chi);    % W_{-1}
W(2) = (1/chi)*W(1) - Wm1;        % W_1
for l = 2 : l_max
    W(l+1) = (2*l-1)/chi * W(l) - W(l-1);
end

% =========================================================
% Accumulate Mie sums — Ulaby (8.33), (8.32), (8.40)
% =========================================================
Qs_sum = 0;  oldQs = 0;
Qe_sum = 0;  oldQe = 0;
Qb_sum = 0;

for l = 1 : l_max

    Wl   = W(l+1);
    Wlm1 = W(l);
    An   = A_down(l);

    a_l = ((An/n  + l/chi)*real(Wl) - real(Wlm1)) / ...
          ((An/n  + l/chi)*Wl       - Wlm1);

    b_l = ((n*An  + l/chi)*real(Wl) - real(Wlm1)) / ...
          ((n*An  + l/chi)*Wl       - Wlm1);

    Qs_sum = Qs_sum + (2*l+1) * (abs(a_l)^2 + abs(b_l)^2);
    Qe_sum = Qe_sum + (2*l+1) * real(a_l + b_l);
    Qb_sum = Qb_sum + (-1)^l  * (2*l+1) * (a_l - b_l);

    if l >= l_min
        pdiff_s = 1; pdiff_e = 1;
        if abs(Qs_sum) > 1e-30
            pdiff_s = abs((Qs_sum - oldQs) / Qs_sum) * 100;
        end
        if abs(Qe_sum) > 1e-30
            pdiff_e = abs((Qe_sum - oldQe) / Qe_sum) * 100;
        end
        if pdiff_s < 0.001 && pdiff_e < 0.001
            break;
        end
    end
    oldQs = Qs_sum;
    oldQe = Qe_sum;

end

Es = (2 / chi^2) * Qs_sum;
Ee = (2 / chi^2) * Qe_sum;
Eb = (1 / chi^2) * abs(Qb_sum)^2;
Ea = Ee - Es;

if Ea < 0 && abs(Ea) < 1e-8 * abs(Es)
    Ea = 0;
end

end


function A1 = lentz_A1(zeta)
% Compute A_1(zeta) = D_1(zeta) via Lentz continued fraction.
% Bohren & Huffman eq. 4.91:
%   D_1 = -1/zeta + 1/(3/zeta - 1/(5/zeta - 1/(7/zeta - ...)))
%
% Equivalently via modified Lentz method (Numerical Recipes):
%   D_1(zeta) = CF where CF is evaluated by forward recurrence.

    TINY = 1e-30;
    max_iter = 10000;
    tol = 1e-12;

    % Continued fraction: A_1 = 1/(zeta*(2/zeta - 1/(3/zeta - ...)))
    % Use the recurrence from B&H eq 4.91 directly:
    % nu = 1.5 (start), increment by 1 each step
    % f_j = (2*nu/zeta) - 1/f_{j-1}

    nu  = 1.5;
    f   = 2*nu/zeta;
    if abs(f) < TINY; f = TINY; end
    C = f;
    D = 0;

    for j = 1 : max_iter
        nu  = nu + 1;
        a_j = -1;
        b_j = 2*nu/zeta;

        D = b_j + a_j * D;
        if abs(D) < TINY; D = TINY; end
        C = b_j + a_j / C;
        if abs(C) < TINY; C = TINY; end

        D = 1/D;
        delta = C * D;
        f = f * delta;

        if abs(delta - 1) < tol
            break;
        end
    end

    % A_1 = 1/zeta - 1/f  (B&H eq. 4.91 rearranged)
    A1 = 1/zeta - 1/f;

end
