%% ========================================================================
%  CW preheating for Gruneisen nonlinear photoacoustics
%  ------------------------------------------------------------------------
%  Purpose
%  This script supports the weak-nonlinearity compensation analysis in the
%  manuscript. It simulates continuous-wave (CW) preheating followed by a
%  nanosecond probe pulse and evaluates the resulting nonlinear PA gain.
%
%  Modeling workflow
%  1. Run one k-Wave reference simulation to obtain the acoustic impulse
%     response, ref_wave.
%  2. Compute the temperature rise caused by CW heating.
%  3. Use the temperature-modulated Gruneisen parameter to generate the
%     probe-pulse source term.
%  4. Plot the temperature rise and nonlinear gain versus heating duration.
%
%  Requirements
%  MATLAB with k-Wave installed. The lowpass filter is applied only when the
%  lowpass function is available.
%% ========================================================================
clear; clc; close all;

verbose = true;
export_results = true;
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end
output_dir = fullfile(script_dir, 'outputs');
if export_results && ~exist(output_dir, 'dir'), mkdir(output_dir); end

set(groot,'defaultTextInterpreter','tex');
set(groot,'defaultAxesTickLabelInterpreter','tex');
set(groot,'defaultLegendInterpreter','tex');
set(groot,'defaultAxesFontName','Arial');

%% 1. Physical and numerical parameters
Gamma0      = 0.12;
dGamma_dT   = 0.01;
rho         = 1000;
Cp          = 3600;
c0          = 1500;
alpha_th    = 1.3e-7;
eta_th      = 1;

F_probe_mJcm2 = 2;      % Probe fluence [mJ/cm^2]
mu_a_cm       = 240;    % Absorption coefficient [1/cm]
w0_um         = 5;      % Heating length scale [um]
tauL_ns       = 10;     % Probe-pulse FWHM [ns]

I_cw_W_cm2 = 100;       % CW heating intensity [W/cm^2]
heating_durations = linspace(0, 300e-6, 50);
probe_at_end_of_heating = true;

CFL = 0.05;
t_end = max(heating_durations) + 50e-6;

%% 2. k-Wave reference impulse response
if verbose, disp('[1/4] Running k-Wave reference simulation...'); end

Nx = 128;
dx = 50e-6;
PML = 20;
sensor_row = PML + 12;

kgrid = kWaveGrid(Nx, dx, Nx, dx);
dt = CFL * dx / c0;
Nt = ceil(t_end / dt);
kgrid.setTime(Nt, dt);
t = kgrid.t_array;

medium.sound_speed = c0;
medium.density = rho;

[xc, yc] = deal(round(Nx/2), round(Nx/2));
src_mask = makeDisc(Nx, Nx, xc, yc, 3);

sensor.mask = zeros(Nx, Nx);
sensor.mask(sensor_row, yc) = 1;
sensor.record = {'p'};

imp = zeros(1, Nt);
imp(1) = 1;
source.p_mask = src_mask;
source.p = repmat(imp, nnz(src_mask), 1);
args = {'PMLSize', PML, 'PlotPML', false, 'Smooth', true, 'DataCast', 'single'};

sdat = kspaceFirstOrder2D(kgrid, medium, source, sensor, args{:});
if isstruct(sdat) && isfield(sdat, 'p')
    ref_wave = double(sdat.p(:)');
else
    ref_wave = double(sdat(:)');
end

fs = 1 / dt;
f_max_expected = c0 / (3 * dx);
ref_wave = optional_lowpass(ref_wave, f_max_expected, fs);
clear sdat source imp;

%% 3. CW preheating scan
if verbose, disp('[2/4] Running CW preheating scan...'); end

F_probe_Jm2 = F_probe_mJcm2 * 10;
I_cw_W_m2   = I_cw_W_cm2 * 1e4;
mu_a_m = mu_a_cm * 100;
w0_m   = w0_um * 1e-6;
tauL_s = tauL_ns * 1e-9;
tau_th = w0_m^2 / alpha_th;

pp_values_kPa = zeros(size(heating_durations));
max_temp_values = zeros(size(heating_durations));

for i = 1:numel(heating_durations)
    current_heating_duration = heating_durations(i);
    if probe_at_end_of_heating
        t_probe = current_heating_duration;
    else
        t_probe = max(heating_durations);
    end
    [pinc_probe, T_at_probe] = build_pincs_cw_heat(t, dt, ...
        current_heating_duration, t_probe, I_cw_W_m2, F_probe_Jm2, ...
        tauL_s, mu_a_m, eta_th, rho, Cp, Gamma0, dGamma_dT, tau_th);
    p_probe_conv = conv(pinc_probe, ref_wave, 'full');
    pp_values_kPa(i) = (max(p_probe_conv) - min(p_probe_conv)) / 1e3;
    max_temp_values(i) = T_at_probe;
end

%% 4. Threshold analysis and figures
if verbose, disp('[3/4] Computing thresholds and generating figures...'); end

pa_linear_ref_kPa = pp_values_kPa(1);
pa_thresh_1pct_kPa = pa_linear_ref_kPa * 1.01;
pa_thresh_5pct_kPa = pa_linear_ref_kPa * 1.05;

[t_reach_1pct_s, T_reach_1pct_K] = find_threshold_crossing(heating_durations, pp_values_kPa, ...
    max_temp_values, pa_thresh_1pct_kPa);
[t_reach_5pct_s, T_reach_5pct_K] = find_threshold_crossing(heating_durations, pp_values_kPa, ...
    max_temp_values, pa_thresh_5pct_kPa);

if verbose
    fprintf('  Thermal relaxation time: %.1f us\n', tau_th * 1e6);
    fprintf('  Probe fluence: %.2f mJ/cm^2\n', F_probe_mJcm2);
    fprintf('  Absorption coefficient: %d 1/cm\n', mu_a_cm);
    fprintf('  CW heating intensity: %d W/cm^2\n', I_cw_W_cm2);
    fprintf('  Linear reference PA signal: %.3f kPa\n', pa_linear_ref_kPa);
    if ~isnan(t_reach_1pct_s)
        fprintf('  1%% gain: %.2f us heating, %.2f K temperature rise.\n', t_reach_1pct_s * 1e6, T_reach_1pct_K);
    else
        fprintf('  1%% gain was not reached within the simulated heating window.\n');
    end
    if ~isnan(t_reach_5pct_s)
        fprintf('  5%% gain: %.2f us heating, %.2f K temperature rise.\n', t_reach_5pct_s * 1e6, T_reach_5pct_K);
    else
        fprintf('  5%% gain was not reached within the simulated heating window.\n');
    end
end

figure('Color', 'w', 'Name', 'Temperature rise versus CW heating duration');
hold on;
plot(heating_durations * 1e6, max_temp_values, '-o', 'LineWidth', 1.5, ...
    'Color', [0.85 0.33 0.1], 'DisplayName', 'Temperature rise at probe time');
if ~isnan(T_reach_1pct_K)
    yline(T_reach_1pct_K, '--', '1% gain threshold', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1.5);
end
if ~isnan(T_reach_5pct_K)
    yline(T_reach_5pct_K, '--', '5% gain threshold', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5);
end
hold off;
grid on; box on;
xlabel('CW heating duration [us]');
ylabel('Temperature rise at probe time [K]');
title('Temperature rise during CW preheating');
legend('Location', 'southeast');

rel_gain_percent = ((pp_values_kPa / pa_linear_ref_kPa) - 1) * 100;
figure('Color', 'w', 'Name', 'Relative nonlinear gain versus CW heating duration');
hold on;
plot(heating_durations * 1e6, rel_gain_percent, '-o', 'LineWidth', 1.5, ...
    'DisplayName', 'Relative nonlinear gain');
yline(1, '--', '1% threshold', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1.5);
yline(5, '--', '5% threshold', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5);
hold off;
grid on; box on;
xlabel('CW heating duration [us]');
ylabel('Relative nonlinear gain, (PA - PA_{linear})/PA_{linear} [%]');
title('Relative nonlinear gain during CW preheating');
legend('Location', 'southeast');
ylim([0, max(max(rel_gain_percent), 6)]);

if export_results
    results = struct();
    results.parameters = struct('F_probe_mJcm2', F_probe_mJcm2, 'mu_a_cm', mu_a_cm, ...
        'w0_um', w0_um, 'tauL_ns', tauL_ns, 'I_cw_W_cm2', I_cw_W_cm2, 'tau_th_s', tau_th);
    results.heating_durations_s = heating_durations;
    results.temperature_rise_K = max_temp_values;
    results.pp_values_kPa = pp_values_kPa;
    results.relative_gain_percent = rel_gain_percent;
    results.t_reach_1pct_s = t_reach_1pct_s;
    results.T_reach_1pct_K = T_reach_1pct_K;
    results.t_reach_5pct_s = t_reach_5pct_s;
    results.T_reach_5pct_K = T_reach_5pct_K;
    save(fullfile(output_dir, 'Fig6_CW_preheating_source_data.mat'), 'results');
    save_all_open_figures(output_dir, 'Fig6_CW_preheating');
end

if verbose, disp('[4/4] CW preheating simulation completed.'); end

%% Local functions
function [p_inc_probe, T_at_probe] = build_pincs_cw_heat(t, dt, ...
    heat_duration, t_probe, I_cw, F_probe, tauL_probe, mu, eta, rho, Cp, G0, dGdT, tau_th)
    q_cw_train = zeros(size(t));
    heating_indices = t >= 0 & t <= heat_duration;
    q_cw_train(heating_indices) = mu * eta * I_cw;
    T_accum = 0;
    T_heat_hist = zeros(size(t));
    for n = 1:numel(t)
        T_heat_hist(n) = T_accum;
        T_accum = T_accum * exp(-dt / tau_th) + q_cw_train(n) * dt / (rho * Cp);
    end
    I_probe = gaussI(F_probe, t, t_probe, tauL_probe, dt);
    q_probe = mu * eta * I_probe;
    [~, probe_idx] = min(abs(t - t_probe));
    T_at_probe = T_heat_hist(probe_idx);
    p_inc_probe = (G0 + dGdT * T_at_probe) .* q_probe * dt;
end

function I = gaussI(F, t, t0, tauFWHM, dt)
    sigma = tauFWHM / (2*sqrt(2*log(2)));
    pulse = exp(-0.5 * ((t - t0) / sigma).^2);
    norm_factor = sum(pulse) * dt;
    if norm_factor > 0
        I = F * pulse / norm_factor;
    else
        I = zeros(size(t));
    end
end

function [t_cross_s, T_cross_K] = find_threshold_crossing(heating_durations, pp_values_kPa, temp_values_K, threshold_kPa)
    idx = find(pp_values_kPa >= threshold_kPa, 1, 'first');
    if ~isempty(idx) && idx > 1
        p1 = pp_values_kPa(idx - 1);
        p2 = pp_values_kPa(idx);
        t1 = heating_durations(idx - 1);
        t2 = heating_durations(idx);
        t_cross_s = t1 + (t2 - t1) * (threshold_kPa - p1) / (p2 - p1);
        T_cross_K = interp1(heating_durations, temp_values_K, t_cross_s);
    else
        t_cross_s = NaN;
        T_cross_K = NaN;
    end
end

function y = optional_lowpass(x, f_cutoff, fs)
    if exist('lowpass', 'file') == 2
        y = lowpass(x, f_cutoff, fs);
    else
        y = x;
        warning('lowpass function not found. Reference impulse response was not filtered.');
    end
end

function save_all_open_figures(output_dir, prefix)
    figs = findall(0, 'Type', 'figure');
    [~, order] = sort([figs.Number]);
    figs = figs(order);
    for i = 1:numel(figs)
        fig = figs(i);
        name = get(fig, 'Name');
        if isempty(name), name = sprintf('figure_%02d', i); end
        safe_name = regexprep(name, '[^A-Za-z0-9_\-]+', '_');
        base = fullfile(output_dir, sprintf('%s_%02d_%s', prefix, i, safe_name));
        exportgraphics(fig, [base '.png'], 'Resolution', 300);
        savefig(fig, [base '.fig']);
    end
end
