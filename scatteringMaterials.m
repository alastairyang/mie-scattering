% Mie Scattering Analysis - MATLAB Version
% Author: Donglai Yang

%% Define Material Properties

% Real permittivity for common materials
eps_r_ice   = 3.15;      % real permittivity of ice
eps_r_dust  = 8.8;       % real permittivity of dust (high dielectric rock)
eps_r_water = 80.1;      % real permittivity of water at room temperature
eps_r_air   = 1.00054;   % real permittivity of air
eps_r_brine = 42;        % real permittivity of brine (deduced from Debye model)

sigma_ice   = 10e-6;    % S/m, conductivity of pure ice
sigma_dust  = 1.5e-3;   % S/m, conductivity of dust
sigma_water = 400e-3;   % S/m, conductivity of water
sigma_air   = 3e-15;    % S/m, conductivity of air
sigma_brine = 11.67;    % S/m, conductivity of brine (effective sigma, deduced from the Debye model)

% Brine temperature
T_brine_K = 250;        % Kelvin

% ---- Material tables ----
% For brine: eps_r and sigma are placeholders (computed per-frequency inside loop)
% sentinel value NaN flags "use brine_parameters()" instead of static values
materials_epsr = {eps_r_brine, 'Brine';
                  eps_r_water, 'Water';
                  eps_r_dust,  'Dust';
                  eps_r_air,   'Air'};
materials_cond = {sigma_brine, 'Brine';
                  sigma_water, 'Water';
                  sigma_dust,  'Dust';
                  sigma_air,   'Air'};

n_materials = size(materials_epsr, 1);

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

%% Specifying parameteric sweep
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

%% ACTUAL PLOT IN THE MANUSCRIPT
% =========================================================
% Aesthetics / style config (set once before the mat loop)
% =========================================================
plot_scattering_flag = true;
plot_absorption_flag = true;
plot_extinction_flag = true;

pub_colors = [
    0.122, 0.471, 0.706;   % blue
    0.890, 0.102, 0.110;   % red
    0.992, 0.553, 0.235;   % orange
    0.416, 0.239, 0.604;   % purple
    0.180, 0.627, 0.173;   % green
];
pub_linestyles    = {'-', '--', ':'};
rayleigh_colors   = [
    0.122, 0.471, 0.706;
    0.890, 0.102, 0.110;
    0.180, 0.627, 0.173;
];

set(groot, 'DefaultAxesFontName',  'Helvetica');
set(groot, 'DefaultTextFontName',  'Helvetica');
set(groot, 'DefaultAxesFontSize',  11);
set(groot, 'DefaultAxesLineWidth', 0.8);

% =========================================================
% Figure + tiled layout
% =========================================================
n_rows = plot_scattering_flag + plot_absorption_flag + plot_extinction_flag;

figure('Position', [1000, 500, 1400, 400 * n_rows]);
t = tiledlayout(n_rows, n_materials, ...
    'TileSpacing', 'none', ...
    'Padding',     'compact');

colors_xi = pub_colors;   % reuse pub palette for phi lines

% =========================================================
% Main loop
% =========================================================
for mat_idx = 1:n_materials
    title_str = materials_epsr{mat_idx, 2};

    % --- Allocate storage ---
    Na_Es_xi_list = cell(1, numel(frequencies));
    Na_Ea_xi_list = cell(1, numel(frequencies));
    Na_Ee_xi_list = cell(1, numel(frequencies));

    Es_xi_list = cell(1, numel(frequencies));
    Ea_xi_list = cell(1, numel(frequencies));
    Ee_xi_list = cell(1, numel(frequencies));

    xi_rayleigh   = zeros(1, numel(frequencies));

    % --- Compute Mie quantities for every frequency ---
    for freq_idx = 1:numel(frequencies)
        lambda = 3e8 / frequencies(freq_idx);
        r_xi   = xi * lambda / (2 * pi);

%         if isnan(materials_epsr{mat_idx, 1})
%             [eps_p, eps_pp] = brine_parameters(T_brine_K, frequencies(freq_idx));
%             epsp = eps_p - 1j * eps_pp;
%         else
        epsp = materials_epsr{mat_idx, 1} + ...
                   compute_imag_permittivity(frequencies(freq_idx), materials_cond{mat_idx, 1});
%         end
        epsb = eps_r_ice + compute_imag_permittivity(frequencies(freq_idx), sigma_ice);

        % Rayleigh boundary: |n| * xi = 0.5  =>  xi_R = 0.5 / |n|
        n_mag                 = abs(sqrt(epsp / epsb));
        xi_rayleigh(freq_idx) = 0.5 / n_mag;

        fprintf("[xi-plot] Material: %s | f = %.0f Hz\n", title_str, frequencies(freq_idx));
        fprintf("  EPSP         : %s\n",   num2str(epsp));
        fprintf("  EPSB         : %s\n",   num2str(epsb));
        fprintf("  n            : %s\n",   num2str(sqrt(epsp / epsb)));
        fprintf("  |n|          : %.4f\n", n_mag);
        fprintf("  xi_Rayleigh  : %.4f\n", xi_rayleigh(freq_idx));
        fprintf("  r range      : [%.3e, %.3e] m\n", r_xi(1), r_xi(end));

        [Na_Es_xi_list{freq_idx}, ...
         Na_Ea_xi_list{freq_idx}, ...
         Na_Ee_xi_list{freq_idx}] = ...
            compute_atten_rate(phi, r_xi, frequencies(freq_idx), epsb, epsp);

        % FOR DEBUGGING
        [Es_xi_list{freq_idx},...
         Ea_xi_list{freq_idx},...
         Ee_xi_list{freq_idx}] = ...
            compute_scat_crossection(r_xi, frequencies(freq_idx), epsb, epsp);

    end

    % --- Build row list based on active flags ---
    row_data   = {};
    row_labels = {};
    if plot_scattering_flag
        row_data{end+1}   = Na_Es_xi_list;
        row_labels{end+1} = 'Scattering (dB km^{-1})';
    end
    if plot_absorption_flag
        row_data{end+1}   = Na_Ea_xi_list;
        row_labels{end+1} = 'Absorption (dB km^{-1})';
    end
    if plot_extinction_flag
        row_data{end+1}   = Na_Ee_xi_list;
        row_labels{end+1} = 'Extinction (dB km^{-1})';
    end

    % --- Draw each row for this material column ---
    for row = 1:n_rows
        tile_idx = (row - 1) * n_materials + mat_idx;
        ax = nexttile(tile_idx);
        hold(ax, 'on');

        is_last_row  = (row == n_rows);
        is_first_col = (mat_idx == 1);
        is_first_row = (row == 1);
        is_last_col  = (mat_idx == n_materials);

        data_list = row_data{row};
        
        grid(ax, 'on');
        ax.GridColor          = [0.85 0.85 0.85];
        ax.GridLineStyle      = '-';
        ax.GridAlpha          = 1.0;
        ax.MinorGridColor     = [0.92 0.92 0.92];
        ax.MinorGridLineStyle = ':';
        ax.MinorGridAlpha     = 1.0;
        set(ax, 'XMinorGrid', 'on', 'YMinorGrid', 'on');

        % --- Data lines (phi x frequency) ---
        for j = 1:length(phi)
            for k = 1:numel(frequencies)
                Na = data_list{k};
                if k == 1
                    plot(ax, xi, Na(:, j), ...
                         'Color',     pub_colors(j, :), ...
                         'LineStyle', pub_linestyles{k}, ...
                         'LineWidth', 1.4, ...
                         'DisplayName', sprintf('\\phi = %.2f', phi(j)));
                else
                    plot(ax, xi, Na(:, j), ...
                         'Color',     pub_colors(j, :), ...
                         'LineStyle', pub_linestyles{k}, ...
                         'LineWidth', 1.4, ...
                         'HandleVisibility', 'off');
                end
            end
        end

        for k = 1:numel(frequencies)
            if is_first_row && is_first_col
                label_str = sprintf('\\xi_R %s', freq_labels{k});
            else
                label_str = '';   % line still drawn, just no text
            end
        
            xline(ax, xi_rayleigh(k), ...
                  pub_linestyles{k}, ...
                  'Color',                    'k', ...
                  'LineWidth',                2.0, ...
                  'Alpha',                    0.55, ...
                  'Label',                    label_str, ...
                  'LabelOrientation',         'aligned', ...
                  'LabelVerticalAlignment',   'top', ...
                  'LabelHorizontalAlignment', 'left', ...
                  'FontSize',                 12, ...
                  'HandleVisibility',         'off');
        end


        % --- Ice background attenuation reference ---
        yline(ax, 5, '--', '5 dB km^{-1}', ...
              'Color',                   [0.4 0.4 0.4], ...
              'LineWidth',               2.0, ...
              'LabelVerticalAlignment',  'bottom', ...
              'LabelHorizontalAlignment','right', ...
              'FontSize',                12, ...
              'HandleVisibility',        'off');

        % --- Axes config ---
        set(ax, 'XScale',      'log', ...
                'YScale',      'log', ...
                'TickDir',     'out', ...
                'TickLength',  [0.015 0.015], ...
                'Box',         'on', ...
                'Layer',       'top');
        xlim(ax, [1e-3, 10]);
        ylim(ax, [1e-5, 1e4]);

        % --- X ticks: first col full range, others drop 10^-3 ---
        if is_first_col
            set(ax, 'XTick', [1e-3, 1e-2, 1e-1, 1e0, 1e1]);
        else
            set(ax, 'XTick', [1e-2, 1e-1, 1e0, 1e1]);
        end

        % --- X label: last row only ---
        if is_last_row
            xlabel(ax, 'k_s = 2\pir/\lambda', 'FontSize', 12);
        else
            set(ax, 'XTickLabel', {});
        end

        % --- Y label: first column only ---
        if is_first_col
            ylabel(ax, row_labels{row}, 'FontSize', 12);
        else
            set(ax, 'YTickLabel', {});
        end

        % --- Title: first row only ---
        if is_first_row
            title(ax, title_str, 'FontSize', 18, 'FontWeight','normal');
        end

        % --- Legend: first row, last column only ---
        if is_first_row && is_last_col
            % Frequency linestyle dummy entries
            for k = 1:numel(frequencies)
                plot(ax, nan, nan, ...
                     'Color',       [0.15 0.15 0.15], ...
                     'LineStyle',   pub_linestyles{k}, ...
                     'LineWidth',   1.4, ...
                     'DisplayName', freq_labels{k});
            end
            lg = legend(ax, ...
                        'Location',  'northwest', ...
                        'FontSize',  10, ...
                        'Box',       'on', ...
                        'EdgeColor', [0.7 0.7 0.7]);
            lg.ItemTokenSize = [18, 9];
        end

        hold(ax, 'off');
    end
end

% --- Shared x-axis label via tiledlayout ---
xlabel(t, 'Size Parameter  k_s = 2\pir/\lambda', 'FontSize', 13);

% =========================================================
% Export
% =========================================================
if plot_scattering_flag && plot_absorption_flag && plot_extinction_flag
    exportgraphics(gcf, 'mie_scattering_absorption_extinction_xiparam.png', 'Resolution', 200);
elseif plot_scattering_flag && plot_absorption_flag
    exportgraphics(gcf, 'mie_scattering_absorption_xiparam.png', 'Resolution', 300);
elseif plot_scattering_flag && plot_extinction_flag
    exportgraphics(gcf, 'mie_scattering_extinction_xiparam.png', 'Resolution', 300);
elseif plot_absorption_flag && plot_extinction_flag
    exportgraphics(gcf, 'mie_absorption_extinction_xiparam.png', 'Resolution', 300);
elseif plot_scattering_flag
    exportgraphics(gcf, 'mie_scattering_xiparam.png', 'Resolution', 300);
elseif plot_absorption_flag
    exportgraphics(gcf, 'mie_absorption_xiparam.png', 'Resolution', 300);
else
    exportgraphics(gcf, 'mie_extinction_xiparam.png', 'Resolution', 300);
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

%% Effective dielectrics for brine
T = 253.15; 
sigma = brine_conductivity(T);
freqs = [1e6,10e6,100e6];

eps_all    = zeros(3,1);
eps_p_all  = zeros(3,1);
eps_pp_all = zeros(3,1);
sigma_eff_all = zeros(3,1); 

for ii = 1:numel(freqs)
    [eps_p, eps_pp, eps_brine_debye] = brine_parameters(T, freqs(ii));
    eps_all(ii) = eps_brine_debye;
    [eps_p_all(ii), eps_pp_all(ii)] = decompose_permittivity(eps_all(ii), sigma, freqs(ii));
    % compute the effective sigma
    sigma_eff_all(ii) = imagi_eps_to_sigma(eps_pp_all(ii), freqs(ii));
end

%

%% Helper Functions

function eps_imag = compute_imag_permittivity(frequency, conductivity)
    eps0    = 8.854e-12;
    omega   = 2 * pi * frequency;
    eps_imag = -1j * conductivity / (omega * eps0);
end

function [Na_s, Na_a, Na_e] = compute_atten_rate(phi, r, f, epsb, epsp)
    if nargin < 4; epsb = 1.78^2; end
    if nargin < 5; epsp = 1;      end

    Es_vec = zeros(1, length(r));
    Ea_vec = zeros(1, length(r));
    Ee_vec = zeros(1, length(r));
 
    for n = 1:length(r)
        [Es, Ea, Ee, ~] = Mie_scattering_new(r(n), f, epsp, epsb);
        Es_vec(n) = Es;
        Ea_vec(n) = Ea;
        Ee_vec(n) = Ee;
    end

    V_particle = (4/3) * pi * r.^3;
    N          = phi ./ V_particle(:);
    Qs         = Es_vec .* pi .* r.^2;
    Qa         = Ea_vec .* pi .* r.^2;
    Qe         = Ee_vec .* pi .* r.^2;
    Qs_array   = repmat(Qs(:), 1, length(phi));
    Qa_array   = repmat(Qa(:), 1, length(phi));
    Qe_array   = repmat(Qe(:), 1, length(phi));
    alpha_s    = N .* Qs_array;
    alpha_a    = N .* Qa_array;
    alpha_e    = N .* Qe_array;
    Na_s       = 10 * log10(exp(1)) * alpha_s * 1e3;
    Na_a       = 10 * log10(exp(1)) * alpha_a * 1e3;
    Na_e       = 10 * log10(exp(1)) * alpha_e * 1e3;
end

function [Es_vec, Ea_vec, Ee_vec] = compute_scat_crossection(r, f, epsb, epsp)
    if nargin < 4; epsb = 1.78^2; end
    if nargin < 5; epsp = 1;      end

    Es_vec = zeros(1, length(r));
    Ea_vec = zeros(1, length(r));
    Ee_vec = zeros(1, length(r));
 
    for n = 1:length(r)
        [Es, Ea, Ee, ~] = Mie_scattering_new(r(n), f, epsp, epsb);
        Es_vec(n) = Es;
        Ea_vec(n) = Ea;
        Ee_vec(n) = Ee;
    end
end


function [eps_p, eps_pp, eps_brine_debye] = brine_parameters(T, frequency)
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

function sigma_eff = imagi_eps_to_sigma(eps_pp, frequency)
    eps0    = 8.854187817e-12;
    omega   = 2 * pi * frequency;
    sigma_eff = eps_pp * (omega * eps0);
end

