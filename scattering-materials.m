% Mie Scattering Analysis - MATLAB Version
% Author: Donglai Yang

%% Define Material Properties

% Real permittivity for common materials
eps_r_ice = 3.15;      % real permittivity of ice
eps_r_dust = 8.8;      % real permittivity of dust (high dielectric rock)
eps_r_ldr = 4.0;       % real permittivity of Low Dielectric Rock (LDR)
eps_r_water = 80.1;    % real permittivity of water at room temperature
eps_r_air = 1.00054;   % real permittivity of air

sigma_ice = 10e-6;  % S/m, conductivity of pure ice
sigma_ldr = 100e-6;  % S/m, conductivity of Low Dielectric Rock (LDR)
sigma_dust = 1.5e-3;  % S/m, conductivity of dust
sigma_water = 400e-3; % S/m, conductivity of water
sigma_air = 3e-15; % S/m, conductivity of air
%% Example 1: Scattering Efficiency vs Pore Radius

epsp = 1;
epsb = 1.78^2;
c = 3e8;

f = 9e6;  % Hz
lambda_wave = c / f;
r = linspace(1, 30, 100);  % 100 points
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
title('Scattering Efficiency vs Pore Radius');

% Display values
fprintf('Es_vec(2:10):\n');
disp(Es_vec(2:10)');
fprintf('\nEe_vec(2:10):\n');
disp(Ee_vec(2:10)');
fprintf('\nEa_vec(2:10):\n');
disp(Ea_vec(2:10)');

%% Example 2: Attenuation Rate for Different Materials at three frequencies

% Common relative permittivity

% scatterer radius
r_sample_N = 1000;
r = linspace(1e-3, 1, r_sample_N);
% r = 1e-3 * [1, 2, 4.5, 10, 22, 46, 100, 220, 460, 1e3];
% porosity
phi = [0.01, 0.02, 0.05, 0.1, 0.20];

% Compute attenuation rate
f_1 = 1e6;      % Hz
f_10 = 10e6;    % Hz
f_100 = 100e6;  % Hz

% Materials and frequencies
materials_epsr = {eps_r_ldr, 'Low dielectric rock'; 
                  eps_r_water, 'Water'; 
                  eps_r_dust, 'Dust';
                  eps_r_air,'Air'};
materials_cond = {sigma_ldr, 'Low dielectric rock';
                  sigma_water, 'Water';
                  sigma_dust, 'Dust';
                  sigma_air, 'Air'};
frequencies = [f_1, f_10, f_100];
freq_labels = {'1 MHz', '10 MHz', '100 MHz'};
linestyles = {'-', '--', ':'}; % one for each frequency
n_materials = size(materials_epsr, 1);

% Create figure with subplots
figure('Position', [100, 100, 1400, 500]);
colors = lines(length(phi));

for mat_idx = 1:n_materials
    subplot(1, n_materials, mat_idx);
    title_str = materials_epsr{mat_idx, 2};
    
    % Compute attenuation for three frequencies
    Na_list = cell(1, 3);
    for freq_idx = 1:numel(frequencies)
        % scatterer dielectric
        epsp = materials_epsr{mat_idx,1} + compute_imag_permittivity(frequencies(freq_idx), materials_cond{mat_idx,1});
        % background ice dielectric
        epsb = eps_r_ice + compute_imag_permittivity(frequencies(freq_idx), sigma_ice);
        fprintf("EPSP: " + num2str(epsp) + "\n")
        fprintf("EPSB: " + num2str(epsb) + "\n")
        Na_list{freq_idx} = compute_atten_rate(phi, r, frequencies(freq_idx), epsb, epsp);
        if freq_idx == 1
            Na_sample = Na_list{freq_idx};
            fprintf("Na: " + num2str(Na_sample(1)) + "\n");
        end
    end
    
    % Plot for each porosity
    hold on;
    for j = 1:length(phi)
        for k = 1:3
            Na = Na_list{k};
            if k == 1
                plot(r * 1e2, Na(:, j), 'Color', colors(j,:), ...
                     'LineStyle', linestyles{k}, 'LineWidth', 1.5, ...
                     'DisplayName', sprintf('\\phi=%.2f', phi(j)));
                set(gca, 'XScale', 'log', 'YScale', 'log');
            else
                plot(r * 1e2, Na(:, j), 'Color', colors(j,:), ...
                     'LineStyle', linestyles{k}, 'LineWidth', 1.5, ...
                     'HandleVisibility', 'off');
                set(gca, 'XScale', 'log', 'YScale', 'log');
            end
        end
    end
    hold off;
    
    % set(gca, 'XScale', 'log', 'YScale', 'log');
    ylim([1e-15, 1e4]);
    xlabel('Scatterer Radius (cm)');
    if mat_idx == 1
        ylabel('Attenuation Rate (dB/km)');
    end
    grid on;
    ax = gca;
    ax.GridLineStyle = '-';
    title(title_str);
    set(gca, 'TickDir', 'in');
    
    % Add legends
    if mat_idx == n_materials
        leg1 = legend('Location', 'southwest');
    else
        leg1 = legend('Location', 'northwest');
    end
    
    % Add frequency legend
    hold on;
    for k = 1:3
        plot(nan, nan, 'k', 'LineStyle', linestyles{k}, ...
             'LineWidth', 1.5, 'DisplayName', freq_labels{k});
    end
    if mat_idx == 1
        leg2 = legend('Location', 'northwest');
    else
        leg2 = legend('Location', 'southeast');
    end
    hold off;
end

% Save figure
exportgraphics(gcf, 'mie_scattering_attenuation.png','Resolution',300);

%% Just plot the scattering efficiency at high frequency (here 100 MHz)
high_freq = 100e6; % 100 MHz
r_sample_N = 200;
r = linspace(1e-3, 1, r_sample_N);
Es_list = zeros(size(materials_epsr,1),numel(r));
Ea_list = zeros(size(materials_epsr,1),numel(r));
Ee_list = zeros(size(materials_epsr,1),numel(r));
Eb_list = zeros(size(materials_epsr,1),numel(r));
colors_e = ["r","b","g","k"]; % each represents a material
materials_legend = [];

figure;
for mat_idx = 1:size(materials_epsr, 1)
    mat_eps = materials_epsr{mat_idx, 1};
    materials_legend = [materials_legend string(materials_epsr{mat_idx, 2})];
    
    % Compute attenuation for three frequencies
    for kk = 1:numel(r)
        [Es_list(mat_idx, kk), Ea_list(mat_idx,kk), Ee_list(mat_idx,kk), Eb_list(mat_idx,kk)] = Mie_scattering(r(kk), high_freq, mat_eps, epsb);
    end
    
    subplot(1,2,1)
    plot(r .* 1e2, Es_list(mat_idx, :), 'Color', colors_e(mat_idx), ...
         'LineWidth',1.5)

    hold on;
    
    set(gca, 'XScale', 'log','YScale','log');
    % set(gca, 'XScale', 'log');    
    ylim([0, 20]);
    xlabel('Scatterer Radius (cm)');
    if mat_idx == 1
        ylabel('Scattering efficiency');
    end
    grid on;
    ax = gca;
    ax.GridLineStyle = '-';
    title("Scattering efficiency (log-space)");
    set(gca, 'TickDir', 'in');
    
    if mat_idx == 3
    legend(materials_legend, 'Location', 'northwest');
    end
   
    % -------
    subplot(1,2,2)
    plot(r .* 1e2, Es_list(mat_idx, :), 'Color', colors_e(mat_idx), ...
         'LineWidth',1.5)

    hold on;
    
    %set(gca, 'XScale', 'log','YScale','log');
    set(gca, 'XScale', 'log');    
    ylim([0, 20]);
    xlabel('Scatterer Radius (cm)');
    if mat_idx == 1
        ylabel('Scattering efficiency');
    end
    grid on;
    ax = gca;
    ax.GridLineStyle = '-';
    title("Scattering efficiency (linear space)");
    set(gca, 'TickDir', 'in');
    
end
sgtitle("Scattering efficiency at 100 MHz")
hold off;

exportgraphics(gcf,'scattering-efficiency-100MHz.png','Resolution',300)

%% Example 3: Single Plot with All Frequencies

% Common relative permittivity
epsp_water = 88;                    % water, at 0 degree C
epsp_dust = 8.8 - 1.7e-2i;         % dust
epsp_salty_ice = 3.5 - 2.7e-20i;   % salty ice

epsp = epsp_salty_ice;
epsb = 1.78^2;
c = 3e8;

r = 1e-3 * [1, 2, 4.5, 10, 22, 46, 100, 220, 460, 1e3];
phi = [0.01, 0.02, 0.05, 0.1, 0.20];

% Compute attenuation rate
f_1 = 1e6;      % Hz
f_10 = 10e6;    % Hz
f_100 = 100e6;  % Hz
Na_1 = compute_atten_rate(phi, r, f_1, epsb, epsp);
Na_10 = compute_atten_rate(phi, r, f_10, epsb, epsp);
Na_100 = compute_atten_rate(phi, r, f_100, epsb, epsp);

% Plot all frequencies on one axes
Na_list = {Na_1, Na_10, Na_100};
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
            plot(r * 1e2, Na(:, j), 'Color', color, ...
                 'LineStyle', linestyles{k}, 'LineWidth', 1.5, ...
                 'DisplayName', sprintf('\\phi=%.2f', phi(j)));
        else
            plot(r * 1e2, Na(:, j), 'Color', color, ...
                 'LineStyle', linestyles{k}, 'LineWidth', 1.5, ...
                 'HandleVisibility', 'off');
        end
    end
end
hold off;

set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Scatterer Radius (cm)');
ylabel('Attenuation Rate (dB/km)');
grid on;
title('Attenuation Rate for Different Frequencies and Porosities');
set(gca, 'TickDir', 'in');

% Legend for porosity
leg1 = legend('Location', 'northwest');
title(leg1, 'Porosity (\phi)');

% Add frequency legend
hold on;
for k = 1:3
    plot(nan, nan, 'k', 'LineStyle', linestyles{k}, ...
         'LineWidth', 1.5, 'DisplayName', freq_labels{k});
end
leg2 = legend('Location', 'northeast');
title(leg2, 'Frequency');
hold off;


%% Helper Functions

function eps_imag = compute_imag_permittivity(frequency, conductivity)
    % Compute the relative and imaginary part of the permittivity from first
    % principles.
    %
    % Parameters
    % ----------  
    % frequency : float
    %     Frequency in Hz.
    % conductivity : float
    %     Conductivity in S/m.
    
    eps_imag = 1j * conductivity / (2 * pi * frequency);
end

function Na = compute_atten_rate(phi, r, f, epsb, epsp)
    % Compute the attenuation rate in dB/km due to Mie scattering 
    %
    % Parameters
    % ---------
    % phi: array
    %     porosity values
    % r: array
    %     scatterer radius in m 
    % f: float
    %     frequency in Hz
    % epsb: float, default=1.78^2
    %     background relative permittivity
    % epsp: float or complex, default=1
    %     scatterer relative permittivity
    %
    % Return
    % ----------
    % Na: matrix
    %     attenuation rate in dB/km (size: length(r) x length(phi))
    
    if nargin < 4
        epsb = 1.78^2;
    end
    if nargin < 5
        epsp = 1;
    end
    
    Es_vec = zeros(1,length(r));
    for n = 1:length(r)
        [Es, ~, ~, ~] = Mie_scattering(r(n), f, epsp, epsb);
        Es_vec(n) = Es;
    end
    % fprintf("Min scattering efficiency: " + num2str(min(Es_vec)) + "\n")
    
    % Compute particle volume
    V_particle = (4/3) * pi * r.^3;
    
    % Make N a (len(r), len(phi)) matrix: porosity divided by particle volume
    N = phi ./ V_particle(:);
    
    % Compute scattering cross-section
    Qs = Es_vec .* pi .* r.^2;
    Qs_array = repmat(Qs(:), 1, length(phi));
    
    % Compute scattering attenuation rate
    alpha_s = N .* Qs_array;
    Na = 10 * log10(exp(1)) * alpha_s * 1e3;
end
