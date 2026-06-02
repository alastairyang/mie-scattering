% Aglyamov2017.m
clear all; close all; clc

%% Figure 1B
epsp = 1;
epsb = 1.78^2;
c = 3e8;

f = 50e6; % Hz
lambda = c/f;
r = linspace(0,10);
Es_vec = zeros(length(r),1);
Ee_vec = Es_vec;
Ea_vec = Es_vec;
for n = 1:length(r)
    [Es, Ea, Ee, Eb] = Mie_scattering(r(n), f, epsp, epsb);
    Es_vec(n) = Es;
    Ee_vec(n) = Ee;
    Ea_vec(n) = Ea;
end

figure
plot(r,Es_vec,'b')
hold on
plot(r,Ee_vec,'r')
plot(r,Ea_vec,'k')

ax = gca;
ax.XLabel.String = 'Pore Radius (m)';
ax.YLabel.String = 'Scattering Efficiency, $\xi$';
legend('$\xi_s$','$\xi_e$','$\xi_a$','Location','NorthWest')

%% Table 2
f = 60e6;

epsp = 1;
epsb = 1.78^2;

r =  1e-3*[1 2 4.5 10 22 46 100 220 460 1e3]';
Ee_vec = zeros(length(r),1);
for n = 1:length(r)
    [~, ~, Ee, ~] = Mie_scattering(r(n), f, epsp, epsb);
    Ee_vec(n) = Ee;
end
phi = [0.01 0.02 0.05 0.1 0.22];
Qe = pi*r.^2.*Ee_vec;
V = 4/3*pi*r.^3;
N = phi./V;

d = 1e3*[2.9 3.3 3.5 3.5 3.1;...
    1.9 2.3 2.7 2.9 2.7;...
    1.2 1.5 1.9 2.2 2.2;...
    0.7 0.9 1.2 1.5 1.6;...
    0.4 0.6 0.7 0.9 1.1;...
    0.4 0.5 0.6 0.7 0.8;...
    0.5 0.5 0.6 0.7 0.8;...
    0.6 0.7 0.8 0.9 1;...
    0.8 1.0 1.1 1.3 1.4;...
    1.0 1.2 1.5 1.8 2.0];

tau = (N.*Qe).*d;
L = -10*log10(exp(-2*tau));
L(L>100) = 100;
disp(round(L,1))

Na_check = zeros(size(N));
for m = 1:length(r)
    for n = 1:length(phi)
        [alpha,Na] = EMscattering(r(m),f,epsp,epsb,phi(n));
        Na_check(m,n) = Na;
    end
end

L_check = 2*Na_check.*d;
L_check(L_check>100) = 100;

%% Table 3
f = 9e6;

epsp = 1;
epsb = 1.78^2;

r =  1e-3*[1 2 4.5 10 22 46 100 220 460 1e3]';
Ee_vec = zeros(length(r),1);
for n = 1:length(r)
    [~, ~, Ee, ~] = Mie_scattering(r(n), f, epsp, epsb);
    Ee_vec(n) = Ee;
end
phi = [0.01 0.02 0.05 0.1 0.22];
Qe = pi*r.^2.*Ee_vec;
V = 4/3*pi*r.^3;
N = phi./V;

d = 1e3*[2.9 3.3 3.5 3.5 3.1;...
    1.9 2.3 2.7 2.9 2.7;...
    1.2 1.5 1.9 2.2 2.2;...
    0.7 0.9 1.2 1.5 1.6;...
    0.4 0.6 0.7 0.9 1.1;...
    0.4 0.5 0.6 0.7 0.8;...
    0.5 0.5 0.6 0.7 0.8;...
    0.6 0.7 0.8 0.9 1;...
    0.8 1.0 1.1 1.3 1.4;...
    1.0 1.2 1.5 1.8 2.0];

tau = (N.*Qe).*d;
L = -10*log10(exp(-2*tau));
L(L>100) = 100;
disp(round(L,1))


%% REASON VHF Scattering Loss
f = 60e6;

epsp = 1;
epsb = 1.78^2;

r = logspace(-3,0,1e3)';
Ee_vec = zeros(length(r),1);
for n = 1:length(r)
    [~, ~, Ee, ~] = Mie_scattering(r(n), f, epsp, epsb);
    Ee_vec(n) = Ee;
end
phi = logspace(-3,0);
d = 1e3; % m
Qe = pi*r.^2.*Ee_vec;
V = 4/3*pi*r.^3;
N = phi./V;
tau = N.*Qe.*d;
L = -10*log10(exp(-2*tau));

x = repmat(phi,length(r),1);
y = repmat(r,1,length(phi));

figure
p = pcolor(x,y,L);
p.LineStyle = 'none';
p.FaceColor = 'interp';
clim([0 100])
axis tight

ax = gca;
ax.XLabel.String = 'Porosity, $\phi$';
ax.YLabel.String = 'Pore Radius, $r$ (m)';
ax.Title.String = 'Volume Scattering Loss, $f=60$ MHz, $d$ = 1 km';

ax.XScale = 'log';
ax.YScale = 'log';
ax.Layer = 'Top';

cb = colorbar;
cb.TickLabelInterpreter = 'latex';
cb.Title.String = {'dB'};
cb.Title.Interpreter = 'latex';

hold on
[C,h] = contour(x,y,L,[1e-3 1e-2 1e-1 1e0],'LineColor','w');
clabel(C,h,'Interpreter','latex','Color','w')

%% REASON HF Scattering Loss
f = 9e6;

epsp = 1;
epsb = 1.78^2;

r = logspace(-3,0,1e3)';
Ee_vec = zeros(length(r),1);
for n = 1:length(r)
    [~, ~, Ee, ~] = Mie_scattering(r(n), f, epsp, epsb);
    Ee_vec(n) = Ee;
end
phi = logspace(-3,0);
d = 1e3; % m
Qe = pi*r.^2.*Ee_vec;
V = 4/3*pi*r.^3;
N = phi./V;
tau = N.*Qe.*d;
L = -10*log10(exp(-2*tau));

x = repmat(phi,length(r),1);
y = repmat(r,1,length(phi));

figure
p = pcolor(x,y,L);
p.LineStyle = 'none';
p.FaceColor = 'interp';
clim([0 100])
axis tight

ax = gca;
ax.XLabel.String = 'Porosity, $\phi$';
ax.YLabel.String = 'Pore Radius, $r$ (m)';
ax.Title.String = 'Volume Scattering Loss, $f=9$ MHz, $d$ = 1 km';

ax.XScale = 'log';
ax.YScale = 'log';
ax.Layer = 'Top';

cb = colorbar;
cb.TickLabelInterpreter = 'latex';
cb.Title.String = {'dB'};
cb.Title.Interpreter = 'latex';

hold on
[C,h] = contour(x,y,L,[1e-3 1e-2 1e-1 1e0],'LineColor','w');
clabel(C,h,'Interpreter','latex','Color','w')
