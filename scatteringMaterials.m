% Mie Scattering Analysis - MATLAB Version
% Author: Donglai Yang

%% Define Material Properties

% Real permittivity for common materials
eps_r_ice   = 3.15;      % real permittivity of ice
eps_r_dust  = 8.8;       % real permittivity of dust (high dielectric rock)
eps_r_water = 80.1;      % real permittivity of water at room temperature
eps_r_air   = 1.00054;   % real permittivity of air

sigma_ice   = 10e-6;    % S/m, conductivity of pure ice
sigma_dust  = 1.5e-3;   % S/m, conductivity of dust
sigma_water = 400e-3;   % S/m, conductivity of water
sigma_air   = 3e-15;    % S/m, conductivity of air

% Brine temperature
T_brine_K = 250;        % Kelvin

%% Example 1: Scattering Efficiency vs Pore Radius

epsp = 100;
epsb = 1.78^2;
c = 3e8;
m = sqrt(epsp/epsb);
fprintf("The refractive index is: %d \n", m)

f = 9e6;  % Hz
lambda_wave = c / f;
r = linspace(1, 30, 500);  % 100 points
Es_vec = zeros(length(r), 1);
Ee_vec = zeros(length(r), 1);
Ea_vec = zeros(length(r), 1);

for n = 1:length(r)
    [Es, Ea, Ee, Eb] = Mie_scattering(r(n), f, epsp, epsb);
    Es_vec(n) = Es;
    Ee_vec(n) = Ee;
    Ea_vec(n) = Ea;
end

% Create figure
figure('Position', [100, 100, 800, 600]);
plot(r, Es_vec, 'b', 'LineWidth', 1.5, 'DisplayName', '\xi_s');
hold on;
plot(r, Ee_vec, 'r', 'LineWidth', 1.5, 'DisplayName', '\xi_e');
plot(r, Ea_vec, 'k', 'LineWidth', 1.5, 'DisplayName', '\xi_a');
hold off;

xlabel('Pore Radius (m)');
ylabel('Scattering Efficiency, \xi');
legend('Location', 'northwest');
grid on;
title(['Scattering Efficiency vs Pore Radius at m = ' num2str(m)]);
exportgraphics(gcf,['scattering_efficiency_m=' num2str(m) '.png'], 'Resolution',300)

% Display values
fprintf('Es_vec(2:10):\n');
disp(Es_vec(2:10)');
fprintf('\nEe_vec(2:10):\n');
disp(Ee_vec(2:10)');
fprintf('\nEa_vec(2:10):\n');
disp(Ea_vec(2:10)');

fprintf("The refractive index is: %d \n", m)

%
r_sample_N = 1000;
r = linspace(1e-3, 1, r_sample_N);
% Porosity
phi = [0.01, 0.02, 0.05, 0.1, 0.20];

% Frequencies
f_1   = 1e6;    % Hz
f_10  = 10e6;   % Hz
f_100 = 100e6;  % Hz
frequencies = [f_1, f_10, f_100];
freq_labels  = {'1 MHz', '10 MHz', '100 MHz'};
linestyles   = {'-', '--', ':'};

% ---- Material tables ----
% For brine: eps_r and sigma are placeholders (computed per-frequency inside loop)
% sentinel value NaN flags "use brine_parameters()" instead of static values
materials_epsr = {NaN,         'Brine';
                  eps_r_water, 'Water';
                  eps_r_dust,  'Dust';
                  eps_r_air,   'Air'};
materials_cond = {NaN,         'Brine';
                  sigma_water, 'Water';
                  sigma_dust,  'Dust';
                  sigma_air,   'Air'};

n_materials = size(materials_epsr, 1);

%% Example 2: Attenuation Rate for Different Materials at Three Frequencies
% ---- Plot mode ----
% 'radius'  → x-axis is scatterer radius in cm
% 'xparam'  → x-axis is size parameter x = 2*pi*r/lambda
plot_mode = 'radius';   % <-- change this to 'xparam' to switch

% ========================================================

figure('Position', [100, 100, 1400, 500]);
colors = lines(length(phi));

for mat_idx = 1:n_materials
    ax = subplot(1, n_materials, mat_idx);
    title_str = materials_epsr{mat_idx, 2};

    Na_list    = cell(1, 3);
    xaxis_list = cell(1, 3);   % store x-axis vector per frequency

    for freq_idx = 1:numel(frequencies)
        if isnan(materials_epsr{mat_idx, 1})
            [eps_p, eps_pp] = brine_parameters(T_brine_K, frequencies(freq_idx));
            epsp = eps_p - 1j * eps_pp;
        else
            epsp = materials_epsr{mat_idx,1} + ...
                   compute_imag_permittivity(frequencies(freq_idx), materials_cond{mat_idx,1});
        end
        epsb = eps_r_ice + compute_imag_permittivity(frequencies(freq_idx), sigma_ice);

        fprintf("Material: %s | f = %.0f Hz\n", title_str, frequencies(freq_idx));
        fprintf("  EPSP: %s\n", num2str(epsp));
        fprintf("  EPSB: %s\n", num2str(epsb));
        m_ref = sqrt(epsp / epsb);
        fprintf("  Refractive index: %s\n", num2str(m_ref));

        [Na_list{freq_idx}, ~] = compute_atten_rate(phi, r, frequencies(freq_idx), epsb, epsp);

        % Build x-axis vector for this frequency
        if strcmp(plot_mode, 'xparam')
            lambda = 3e8 / frequencies(freq_idx);   % wavelength in vacuum (m)
            xaxis_list{freq_idx} = 2 * pi * r / lambda;
        else
            xaxis_list{freq_idx} = r * 1e2;         % radius in cm
        end

        if freq_idx == 1
            fprintf("  Na(1): %s\n", num2str(Na_list{freq_idx}(1)));
        end
    end

    % --- Plot ---
    hold(ax, 'on');
    for j = 1:length(phi)
        for k = 1:3
            Na    = Na_list{k};
            xaxis = xaxis_list{k};

            if k == 1
                plot(ax, xaxis, Na(:, j), 'Color', colors(j,:), ...
                     'LineStyle', linestyles{k}, 'LineWidth', 1.5, ...
                     'DisplayName', sprintf('\\phi=%.2f', phi(j)));
            else
                plot(ax, xaxis, Na(:, j), 'Color', colors(j,:), ...
                     'LineStyle', linestyles{k}, 'LineWidth', 1.5, ...
                     'HandleVisibility', 'off');
            end
        end
    end

    % --- Axes config ---
    ylim(ax, [1e-15, 1e4]);
    if strcmp(plot_mode, 'xparam')
        xlabel(ax, 'Size Parameter, x = 2\pir/\lambda');
    else
        xlabel(ax, 'Scatterer Radius (cm)');
    end
    if mat_idx == 1
        ylabel(ax, 'Attenuation Rate (dB/km)');
    end
    grid(ax, 'on');
    ax.GridLineStyle = '-';
    title(ax, title_str);
    set(ax, 'TickDir', 'in');

    % Porosity legend (first call — only phi entries are in legend so far)
    legend(ax, 'Location', 'southwest');

    % Frequency dummy lines
    for k = 1:3
        plot(ax, nan, nan, 'k', 'LineStyle', linestyles{k}, ...
             'LineWidth', 1.5, 'DisplayName', freq_labels{k});
    end
    legend(ax, 'Location', 'southeast');
    set(gca,'XScale','log','YScale','log')

    hold(ax, 'off');
end

% Build filename dynamically so you don't overwrite the wrong one
if strcmp(plot_mode, 'xparam')
    exportgraphics(gcf, 'mie_scattering_attenuation_xparam.png', 'Resolution', 300);
else
    exportgraphics(gcf, 'mie_scattering_attenuation.png', 'Resolution', 300);
end

%% ============================================================
%  Normalized Size-Parameter Plot  (standalone block)
%  X-axis: xi = 2*pi*r/lambda, swept from 1e-3 to 10
%  For each frequency, physical r is back-computed from xi
% =============================================================
normalize_attenuation = false;   % true  → N_a · r  (dB/km · m)
                                 % false → N_a      (dB/km)

% ---- Sweep in size parameter space ----
xi_sample_N = 1000;
xi          = logspace(-3, 1, xi_sample_N);
phi         = [0.01, 0.02, 0.05, 0.1, 0.20];

figure('Position', [1000, 500, 1600, 400]);
colors_xi = lines(length(phi));

for mat_idx = 1:n_materials
    ax        = subplot(1, n_materials, mat_idx);
    title_str = materials_epsr{mat_idx, 2};

    Na_xi_list     = cell(1, numel(frequencies));
    lambda_list    = zeros(1, numel(frequencies));   % store lambda per freq

    for freq_idx = 1:numel(frequencies)

        % --- Back-compute physical r from xi ---
        lambda_list(freq_idx) = 3e8 / frequencies(freq_idx);
        r_xi = xi * lambda_list(freq_idx) / (2 * pi);

        % --- Permittivities (real-only special case) ---
        if isnan(materials_epsr{mat_idx, 1})
            [eps_p, ~] = brine_parameters(T_brine_K, frequencies(freq_idx));
            epsp = eps_p;
        else
            epsp = materials_epsr{mat_idx, 1};
        end
        epsb = eps_r_ice;

        fprintf("[xi-plot] Material: %s | f = %.0f Hz\n", title_str, frequencies(freq_idx));
        fprintf("  EPSP: %s\n",  num2str(epsp));
        fprintf("  EPSB: %s\n",  num2str(epsb));
        fprintf("  m   : %s\n",  num2str(sqrt(epsp / epsb)));
        fprintf("  r range: [%.3e, %.3e] m\n", r_xi(1), r_xi(end));

        [Na_xi_list{freq_idx}, ~] = compute_atten_rate(phi, r_xi, frequencies(freq_idx), epsb, epsp);
    end

    % --- Plot ---
    hold(ax, 'on');
    for j = 1:length(phi)
        for k = 1:numel(frequencies)
            Na      = Na_xi_list{k};
            r_xi_k  = xi * lambda_list(k) / (2 * pi);   % (m), column after (:)

            if normalize_attenuation
                y_data = Na(:, j) .* r_xi_k(:);   % N_a · r  (dB/km · m)
            else
                y_data = Na(:, j);                 % N_a      (dB/km)
            end

            if k == 1
                plot(ax, xi, y_data, ...
                     'Color',            colors_xi(j, :), ...
                     'LineStyle',        linestyles{k}, ...
                     'LineWidth',        1.5, ...
                     'DisplayName',      sprintf('\\phi=%.2f', phi(j)));
            else
                plot(ax, xi, y_data, ...
                     'Color',            colors_xi(j, :), ...
                     'LineStyle',        linestyles{k}, ...
                     'LineWidth',        1.5, ...
                     'HandleVisibility', 'off');
            end
        end
    end

    % --- Axes config ---
    set(ax, 'XScale', 'log', 'YScale', 'log', 'TickDir', 'in');
    xlim(ax, [1e-3, 10]);
    if normalize_attenuation
        ylim(ax, [1e-5, 1e4]);
    else
        ylim(ax, [1e-5, 1e4]);
    end
    grid(ax, 'on');
    ax.GridLineStyle = '-';

    xline(ax, 1, '--k', '\xi = 1', ...
          'LabelVerticalAlignment', 'bottom', ...
          'HandleVisibility',       'off');

    xlabel(ax, 'Size Parameter  \xi = 2\pir/\lambda');
    if mat_idx == 1
        if normalize_attenuation
            ylabel(ax, 'Normalized Attenuation  N_a \cdot r  (dB/km \cdot m)');
        else
            ylabel(ax, 'Attenuation Rate  N_a  (dB/km)');
        end
    end
    title(ax, title_str);

    % Porosity legend entries (already in legend from k==1 plots above)
    legend(ax, 'Location', 'northwest');

    % Frequency dummy lines
    for k = 1:numel(frequencies)
        plot(ax, nan, nan, 'k', ...
             'LineStyle',  linestyles{k}, ...
             'LineWidth',  1.5, ...
             'DisplayName', freq_labels{k});
    end
    legend(ax, 'Location', 'northwest');

    hold(ax, 'off');
end

% --- Export with descriptive filename ---
if normalize_attenuation
    exportgraphics(gcf, 'mie_xiparam_normalized.png', 'Resolution', 300);
else
    exportgraphics(gcf, 'mie_xiparam_raw.png',        'Resolution', 300);
end



%% --- Difference Plot: Brine minus Water ---

% Indices in your materials table
idx_brine = 1;
idx_water = 2;

figure('Position', [100, 100, 600, 500]);
colors_diff = lines(length(phi));

% Recompute Na for brine and water across all frequencies
Na_brine = cell(1, 3);
Na_water = cell(1, 3);

for freq_idx = 1:numel(frequencies)

    % --- Brine ---
    [eps_p, eps_pp] = brine_parameters(T_brine_K, frequencies(freq_idx));
    epsp_brine = eps_p - 1j * eps_pp;
    epsb = eps_r_ice + compute_imag_permittivity(frequencies(freq_idx), sigma_ice);
    Na_brine{freq_idx} = compute_atten_rate(phi, r, frequencies(freq_idx), epsb, epsp_brine);

    % --- Water ---
    epsp_water = eps_r_water + compute_imag_permittivity(frequencies(freq_idx), sigma_water);
    Na_water{freq_idx} = compute_atten_rate(phi, r, frequencies(freq_idx), epsb, epsp_water);

end

ax_diff = axes;
hold(ax_diff, 'on');

for j = 1:length(phi)
    for k = 1:3
        delta_Na = Na_brine{k}(:, j) - Na_water{k}(:, j);

        if k == 1
            plot(ax_diff, r * 1e2, delta_Na, ...
                 'Color', colors_diff(j,:), ...
                 'LineStyle', linestyles{k}, 'LineWidth', 2.5, ...
                 'DisplayName', sprintf('\\phi=%.2f', phi(j)));
        else
            plot(ax_diff, r * 1e2, delta_Na, ...
                 'Color', colors_diff(j,:), ...
                 'LineStyle', linestyles{k}, 'LineWidth', 2.5, ...
                 'HandleVisibility', 'off');
        end
    end
end

% Frequency dummy lines
for k = 1:3
    plot(ax_diff, nan, nan, 'k', 'LineStyle', linestyles{k}, ...
         'LineWidth', 1.5, 'DisplayName', freq_labels{k});
end

hold(ax_diff, 'off');

set(ax_diff, 'XScale', 'log', 'TickDir', 'in');
grid(ax_diff, 'on');
ax_diff.GridLineStyle = '-';

% ylim([0,50])
xlabel(ax_diff, 'Scatterer Radius (cm)');
ylabel(ax_diff, '\Delta Attenuation Rate (dB/km)');
title(ax_diff, 'Attenuation Rate Difference: Brine - Water');
legend(ax_diff, 'Location', 'northwest');

exportgraphics(gcf, 'mie_diff_brine_water.png', 'Resolution', 300);

%% Just Plot Scattering Efficiency at 100 MHz

high_freq  = 100e6;
r_sample_N = 200;
r = linspace(1e-3, 1, r_sample_N);

Es_list = zeros(size(materials_epsr,1), numel(r));
Ea_list = zeros(size(materials_epsr,1), numel(r));
Ee_list = zeros(size(materials_epsr,1), numel(r));
Eb_list = zeros(size(materials_epsr,1), numel(r));

colors_e       = ["r","b","g","k"];
materials_legend = [];

% Background ice at 100 MHz (needed for Mie call)
epsb_100 = eps_r_ice + compute_imag_permittivity(high_freq, sigma_ice);

figure;
for mat_idx = 1:size(materials_epsr, 1)
    materials_legend = [materials_legend string(materials_epsr{mat_idx, 2})];

    % Build complex permittivity at 100 MHz
    if isnan(materials_epsr{mat_idx, 1})
        [eps_p, eps_pp] = brine_parameters(T_brine_K, high_freq);
        mat_eps = eps_p - 1j * eps_pp;
    else
        mat_eps = materials_epsr{mat_idx, 1} + ...
                  compute_imag_permittivity(high_freq, materials_cond{mat_idx, 1});
    end

    for kk = 1:numel(r)
        [Es_list(mat_idx,kk), Ea_list(mat_idx,kk), ...
         Ee_list(mat_idx,kk), Eb_list(mat_idx,kk)] = ...
            Mie_scattering(r(kk), high_freq, mat_eps, epsb_100);
    end

    subplot(1,2,1)
    loglog(r .* 1e2, Es_list(mat_idx, :), 'Color', colors_e(mat_idx), 'LineWidth', 1.5)
    hold on;
    % set(gca, 'XScale', 'log', 'YScale', 'log');
    ylim([0, 20]);
    xlabel('Scatterer Radius (cm)');
    if mat_idx == 1; ylabel('Scattering efficiency'); end
    grid on;
    ax = gca;
    ax.GridLineStyle = '-';
    title("...");
    %set(ax, 'TickDir', 'in');

    if mat_idx == size(materials_epsr,1)
        legend(materials_legend, 'Location', 'northwest');
    end

    subplot(1,2,2)
    loglog(r .* 1e2, Es_list(mat_idx, :), 'Color', colors_e(mat_idx), 'LineWidth', 1.5)
    hold on;
    % set(gca, 'XScale', 'log');
    ylim([0, 20]);
    xlabel('Scatterer Radius (cm)');
    if mat_idx == 1; ylabel('Scattering efficiency'); end
    grid on;
    ax = gca;
    ax.GridLineStyle = '-';
    title("...");
    %set(ax, 'TickDir', 'in');

end
sgtitle("Scattering efficiency at 100 MHz");
hold off;
exportgraphics(gcf, 'scattering-efficiency-100MHz.png', 'Resolution', 300);

%% Example 3: Single Plot with All Frequencies

% Brine permittivity at each frequency (for Example 3, use f_1 as reference display)
% epsp is rebuilt per-frequency inside Na computation via compute_atten_rate,
% but Example 3 uses a fixed epsp — so we pick 10 MHz as a representative value.
[eps_p_ex3, eps_pp_ex3] = brine_parameters(T_brine_K, f_10);
epsp = eps_p_ex3 - 1j * eps_pp_ex3;
epsb = 1.78^2;

r   = 1e-3 * [1, 2, 4.5, 10, 22, 46, 100, 220, 460, 1e3];
phi = [0.01, 0.02, 0.05, 0.1, 0.20];

Na_1   = compute_atten_rate(phi, r, f_1,   epsb, epsp);
Na_10  = compute_atten_rate(phi, r, f_10,  epsb, epsp);
Na_100 = compute_atten_rate(phi, r, f_100, epsb, epsp);

Na_list    = {Na_1, Na_10, Na_100};
linestyles = {'-', '--', ':'};
freq_labels = {'1 MHz', '10 MHz', '100 MHz'};

figure('Position', [100, 100, 900, 700]);
colors = lines(length(phi));
hold on;
for j = 1:length(phi)
    color = colors(j, :);
    for k = 1:length(Na_list)
        Na = Na_list{k};
        if k == 1
            loglog(r * 1e2, Na(:, j), 'Color', color, ...
                 'LineStyle', linestyles{k}, 'LineWidth', 1.5, ...
                 'DisplayName', sprintf('\\phi=%.2f', phi(j)));
        else
            loglog(r * 1e2, Na(:, j), 'Color', color, ...
                 'LineStyle', linestyles{k}, 'LineWidth', 1.5, ...
                 'HandleVisibility', 'off');
        end
    end
end
hold off;

% set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Scatterer Radius (cm)');
ylabel('Attenuation Rate (dB/km)');
grid on;
title('Attenuation Rate for Brine Inclusions — Different Frequencies and Porosities');
%set(gca, 'TickDir', 'in');

leg1 = legend('Location', 'northwest');
title(leg1, 'Porosity (\phi)');

hold on;
for k = 1:3
    loglog(nan, nan, 'k', 'LineStyle', linestyles{k}, ...
         'LineWidth', 1.5, 'DisplayName', freq_labels{k});
end
leg2 = legend('Location', 'northeast');
title(leg2, 'Frequency');
hold off;


%% Helper Functions

function eps_imag = compute_imag_permittivity(frequency, conductivity)
    eps0    = 8.854e-12;
    omega   = 2 * pi * frequency;
    eps_imag = -1j * conductivity / (omega * eps0);
end

function [Na, Es_vec] = compute_atten_rate(phi, r, f, epsb, epsp)
    if nargin < 4; epsb = 1.78^2; end
    if nargin < 5; epsp = 1;      end

    Es_vec = zeros(1, length(r));
    for n = 1:length(r)
        [Es, ~, ~, ~] = Mie_scattering(r(n), f, epsp, epsb);
        Es_vec(n) = Es;
    end

    V_particle = (4/3) * pi * r.^3;
    N          = phi ./ V_particle(:);
    Qs         = Es_vec .* pi .* r.^2;
    Qs_array   = repmat(Qs(:), 1, length(phi));
    alpha_s    = N .* Qs_array;
    Na         = 10 * log10(exp(1)) * alpha_s * 1e3;
end

function [eps_p, eps_pp] = brine_parameters(T, frequency)
    sigma          = brine_conductivity(T);
    eps_brine_debye = brine_permittivity(T, frequency);
    [eps_p, eps_pp] = decompose_permittivity(eps_brine_debye, sigma, frequency);
end

function [eps_prime, eps_double_prime] = decompose_permittivity(eps_debye, sigma, frequency)
    eps0    = 8.854187817e-12;
    omega   = 2 * pi * frequency;
    eps_total        = eps_debye - 1j * sigma / (omega * eps0);
    eps_prime        =  real(eps_total);
    eps_double_prime = -imag(eps_total);
end
