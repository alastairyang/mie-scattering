function [Es, Ea, Ee, Eb] = Mie_scattering_new(r, f, epsp, epsb)
% Calculates the Mie scattering efficiency factors resulting from spherical
% inclusions embedded in a medium.
%
% Syntax:
% [Es Ea Ee Eb] = Mie_scattering(r, f, epsp, epsb)
%
% Inputs:
% r         Particle radius (m)
% f         Frequency (Hz)
% epsp      Particle relative permittivity
% epsb      Background relative permittivity
%
% Outputs:
% Es        Scattering efficiency factor
% Ea        Absorption efficiency factor
% Ee        Extinction efficiency
% Eb        Backscattering efficiency
%
% Source:
% Ulaby and Long (2014)
%
% Author:
% Ulaby and Long (2014) modified by Natalie Wolfenbarger
% nswolfen@gmail.com
%
%%
c = 3e8; % m/s
epsb_real = real(epsb);

np = sqrt(epsp); % index of refraction of spherical particle
nb = sqrt(epsb); % index of refraction of background medium

n = np./nb; % relative index of refraction (8.31a)



lambda = c/f;
chi = (2*pi*r/lambda)*sqrt(epsb_real); % normalized circumference in background (8.31b)

% Values of W0 and W-1
W_1 = sin(chi)+1i*cos(chi); % (8.35a)
W_2 = cos(chi)-1i*sin(chi); % (8.35b)

% Value of A0
A = cot(n*chi); % (8.37)

oldSum = 0;
pdiff = 1;
l = 1;
while pdiff>=0.001
    W = (2*l-1)/chi*W_1-W_2; % (8.34)

    A = -l/(n*chi)+(l/(n*chi)-A)^-1; % (8.36)

    a = ((A/n+l/chi)*real(W)-real(W_1))/((A/n+l/chi)*W-W_1); % (8.33a)
    b = ((n*A+l/chi)*real(W)-real(W_1))/((n*A+l/chi)*W-W_1); % (8.33b)

    sigma = (2*l+1)*(abs(a)^2+abs(b)^2); % (8.32a)
    newSum = oldSum+sigma;

    l = l+1;
    W_2 = W_1;
    W_1 = W;

    pdiff=abs((newSum-oldSum)/newSum)*100;
    oldSum = newSum;

end
Es = 2/(chi)^2*newSum; % (8.32a)

% Values of W0 and W-1
W_1 = sin(chi)+1i*cos(chi); % (8.35a)
W_2 = cos(chi)-1i*sin(chi); % (8.35b)

% Value of A0
A = cot(n*chi); % (8.37)

oldSum = 0;
pdiff = 1;
l = 1;
while pdiff>=0.001
    W = (2*l-1)/chi*W_1-W_2; % (8.34)

    A = -l/(n*chi)+(l/(n*chi)-A)^-1; % (8.36)

    a = ((A/n+l/chi)*real(W)-real(W_1))/((A/n+l/chi)*W-W_1); % (8.33a)
    b = ((n*A+l/chi)*real(W)-real(W_1))/((n*A+l/chi)*W-W_1); % (8.33b)

    sigma = (2*l+1)*real(a+b); % (8.32b)
    newSum = oldSum+sigma;

    l = l+1;
    W_2=W_1;
    W_1=W;

    pdiff=abs((newSum-oldSum)/newSum)*100;
    oldSum = newSum;

end
Ee = 2/(chi)^2*newSum; % (8.32b)

% Values of W0 and W-1
W_1 = sin(chi)+1i*cos(chi); % (8.35a)
W_2 = cos(chi)-1i*sin(chi); % (8.35b)

% Value of A0
A = cot(n*chi); % (8.37)

oldSum = 0;
pdiff = 1;
l = 1;
while pdiff>=0.001
    W = (2*l-1)/chi*W_1-W_2; % (8.34)

    A = -l/(n*chi)+(l/(n*chi)-A)^-1; % (8.36)

    a = ((A/n+l/chi)*real(W)-real(W_1))/((A/n+l/chi)*W-W_1); % (8.33a)
    b = ((n*A+l/chi)*real(W)-real(W_1))/((n*A+l/chi)*W-W_1); % (8.33b)

    sigma = (-1)^l*(2*l+1)*(a-b); % (8.32b)
    newSum = oldSum+sigma;

    l = l+1;
    W_2=W_1;
    W_1=W;

    pdiff=abs((newSum-oldSum)/newSum)*100;
    oldSum = newSum;

end
Eb = 1/(chi)^2*abs(newSum)^2; % (8.40)

Ea = Ee-Es; % (8.28b)

end

