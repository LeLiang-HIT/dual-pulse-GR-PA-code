%% ========================================================================
%  Multi-pulse heat accumulation for Gruneisen nonlinear photoacoustics
%  ------------------------------------------------------------------------
%  Purpose
%  This script supports the weak-nonlinearity compensation analysis in the
%  manuscript. It simulates residual heat accumulation from a train of
%  preheating pulses and evaluates the nonlinear enhancement of a subsequent
%  probe pulse.
%
%  Modeling workflow
%  1. Run one k-Wave reference simulation to obtain a bipolar acoustic
%     impulse response, ref_wave.
%  2. Build the initial-pressure source term for each pulse sequence.
%  3. Compute detected PA waveforms by convolving the source term with
%     ref_wave.
%  4. Plot the nonlinear gain pathway versus the number of heating pulses.
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

F_mJcm2     = 2;       % Fluence per pulse [mJ/cm^2]
mu_a_cm     = 240;     % Absorption coefficient [1/cm]
w0_um       = 5;       % Heating length scale [um]
tauL_ns     = 10;      % Optical-pulse FWHM [ns]

t_start     = 10e-6;   % Time of the first heating pulse [s]
t_rep_s     = 30e-6;   % Pulse interval: 30 us, corresponding to 33.3 kHz
N_max       = 15;      % Maximum number of heating pulses

CFL         = 0.1;
t_end       = t_start + (N_max + 1) * t_rep_s;

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
source.p_mode = 'additive';
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

%% 3. Multi-pulse heat-accumulation simulation
if verbose, disp('[2/4] Running multi-pulse heat-accumulation scan...'); end

F_Jm2  = F_mJcm2 * 10;
mu_a_m = mu_a_cm * 100;
w0_m   = w0_um * 1e-6;
tauL_s = tauL_ns * 1e-9;
tau_th = w0_m^2 / alpha_th;

p_probe_all = cell(1, N_max);
pp_values_kPa = zeros(1, N_max);
delta_pa_path = zeros(1, N_max);

[pinc_linear_ref, ~] = build_pincs_cumulative(t, dt, 0, ...
    F_Jm2, t_rep_s, F_Jm2, t_start, tauL_s, mu_a_m, eta_th, rho, Cp, ...
    Gamma0, dGamma_dT, tau_th);
p_linear_ref = conv(pinc_linear_ref, ref_wave, 'full');
PA_linear_ref_val = max(p_linear_ref) - min(p_linear_ref);
t_long_linear = (0:numel(p_linear_ref)-1) * dt;

for n = 1:N_max
    [~, pinc_probe_n] = build_pincs_cumulative(t, dt, n, ...
        F_Jm2, t_rep_s, F_Jm2, t_start, tauL_s, mu_a_m, eta_th, rho, Cp, ...
        Gamma0, dGamma_dT, tau_th);
    p_probe_n_conv = conv(pinc_probe_n, ref_wave, 'full');
    p_probe_all{n} = p_probe_n_conv;
    pp_values_kPa(n) = (max(p_probe_n_conv) - min(p_probe_n_conv)) / 1e3;
    delta_pa_path(n) = (max(p_probe_n_conv) - min(p_probe_n_conv)) - PA_linear_ref_val;
end

rel_gain_path = delta_pa_path / PA_linear_ref_val;
[F_weak_mJcm2, F_strong_mJcm2] = F_thresholds(mu_a_m, Gamma0, dGamma_dT, rho, Cp, t_rep_s, tau_th);
N_strong_idx = find(rel_gain_path >= 0.05, 1, 'first');

if verbose
    fprintf('  Thermal relaxation time: %.1f us\n', tau_th * 1e6);
    fprintf('  1%% fluence threshold: %.2f mJ/cm^2\n', F_weak_mJcm2);
    fprintf('  5%% fluence threshold: %.2f mJ/cm^2\n', F_strong_mJcm2);
    if ~isempty(N_strong_idx)
        fprintf('  At F = %.1f mJ/cm^2, 5%% gain is reached after %d heating pulses.\n', F_mJcm2, N_strong_idx);
    end
end

%% 4. Figures
if verbose, disp('[3/4] Generating figures...'); end

t_long_final = (0:numel(p_probe_all{end})-1) * dt;
[~, idx_final_peak] = max(p_probe_all{N_max});
t_final_arrival_us = t_long_final(idx_final_peak) * 1e6;
plot_end_time_us = t_final_arrival_us + 50;

figure('Color','w', 'Name', 'Multi-pulse signal evolution');
ax1 = axes;
hold(ax1, 'on');
plot(ax1, t_long_linear * 1e6, p_linear_ref / 1e3, '--', 'Color', [0.5 0.5 0.5], ...
    'LineWidth', 2, 'DisplayName', 'Linear reference (N=0)');
colors = interp1([1, N_max], [0.1 0.1 0.8; 0.8 0.1 0.1], 1:N_max);
for n = 1:N_max
    t_long_n = (0:numel(p_probe_all{n})-1) * dt;
    if n == 1 || n == round(N_max/2) || n == N_max
        plot(ax1, t_long_n * 1e6, p_probe_all{n} / 1e3, '-', 'Color', colors(n,:), ...
            'LineWidth', 1.5, 'DisplayName', sprintf('Probe signal (N=%d)', n));
    else
        plot(ax1, t_long_n * 1e6, p_probe_all{n} / 1e3, '-', 'Color', colors(n,:), ...
            'LineWidth', 1.5, 'HandleVisibility', 'off');
    end
end
hold(ax1, 'off');
grid(ax1, 'on'); box(ax1, 'on');
xlabel(ax1, 'Time [us]');
ylabel(ax1, 'Pressure [kPa]');
title(ax1, 'Signal enhancement by multi-pulse heat accumulation');
colormap(ax1, colors);
cb = colorbar(ax1);
cb.Label.String = 'Number of heating pulses, N';
cb.Ticks = [0, 0.5, 1];
cb.TickLabels = {1, round(N_max/2), N_max};
legend(ax1, 'show', 'Location', 'northeast');
xlim(ax1, [-10, plot_end_time_us]);

figure('Color','w', 'Name', 'Multi-pulse nonlinear gain pathway');
hold on;
h_main = plot(1:N_max, rel_gain_path * 100, '-o', 'LineWidth', 2, ...
    'MarkerFaceColor', 'auto', 'DisplayName', 'Relative nonlinear gain');
h_1pct = yline(1, '--', '1% threshold', 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 1.5);
h_5pct = yline(5, '--', '5% threshold', 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 1.5);
hold off;
grid on; box on;
xlabel('Number of heating pulses, N');
ylabel('Relative nonlinear gain, \DeltaPA/PA_{linear} [%]');
legend([h_main, h_1pct, h_5pct], 'Location', 'southeast');
xticks(0:max(1,floor(N_max/5)):N_max);
ylim([0, max(max(rel_gain_path)*100, 6)]);

final_heated_signal = p_probe_all{N_max};
[~, idx_linear_peak] = max(p_linear_ref);
t_peak_linear_us = t_long_linear(idx_linear_peak) * 1e6;
time_shift_us = t_final_arrival_us - t_peak_linear_us;
t_final_shifted_us = t_long_final * 1e6 - time_shift_us;

figure('Color', 'w', 'Name', 'Linear and final heated signal comparison');
hold on;
plot(t_long_linear * 1e6, p_linear_ref / 1e3, '-', 'LineWidth', 1.5, 'DisplayName', 'Linear reference (N=0)');
plot(t_final_shifted_us, final_heated_signal / 1e3, '-', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Final heated signal (N=%d)', N_max));
hold off;
grid on; box on;
zoom_window_us = 1.0;
xlim([t_peak_linear_us - zoom_window_us, t_peak_linear_us + zoom_window_us]);
xlabel('Aligned time [us]');
ylabel('Pressure [kPa]');
title('Waveform comparison after multi-pulse accumulation');
legend('Location', 'northwest');

figure('Color', 'w', 'Name', 'Peak-to-peak signal versus heating-pulse number');
plot(1:N_max, pp_values_kPa, '-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'auto');
grid on; box on;
xlabel('Number of heating pulses, N');
ylabel('Peak-to-peak signal [kPa]');
title('Peak-to-peak signal versus number of heating pulses');
xlim([0, N_max + 1]);

if export_results
    results = struct();
    results.parameters = struct('F_mJcm2', F_mJcm2, 'mu_a_cm', mu_a_cm, 'w0_um', w0_um, ...
        'tauL_ns', tauL_ns, 't_rep_s', t_rep_s, 'N_max', N_max, 'tau_th_s', tau_th);
    results.pp_values_kPa = pp_values_kPa;
    results.rel_gain_path = rel_gain_path;
    results.F_threshold_1pct_mJcm2 = F_weak_mJcm2;
    results.F_threshold_5pct_mJcm2 = F_strong_mJcm2;
    results.N_to_5pct = N_strong_idx;
    save(fullfile(output_dir, 'Fig6_multi_pulse_source_data.mat'), 'results');
    save_all_open_figures(output_dir, 'Fig6_multi_pulse');
end

if verbose, disp('[4/4] Multi-pulse simulation completed.'); end

%% Local functions
function [p_inc_linear, p_inc_probe] = build_pincs_cumulative(...
        t, dt, N_heat, F_heat, t_rep, F_probe, t_start, tauL, mu, eta, rho, Cp, G0, dGdT, tau_th)
    q_heat_train = zeros(size(t));
    if N_heat > 0
        for i = 1:N_heat
            t_pulse_i = t_start + (i-1) * t_rep;
            I_heat_i = gaussI(F_heat, t, t_pulse_i, tauL, dt);
            q_heat_train = q_heat_train + mu * eta * I_heat_i;
        end
    end
    T_accum = 0;
    T_heat_hist = zeros(size(t));
    for n = 1:numel(t)
        T_heat_hist(n) = T_accum;
        T_accum = T_accum * exp(-dt / tau_th) + q_heat_train(n) * dt / (rho * Cp);
    end
    t_probe = t_start + N_heat * t_rep;
    I_probe = gaussI(F_probe, t, t_probe, tauL, dt);
    q_probe = mu * eta * I_probe;
    p_inc_probe = (G0 + dGdT * T_heat_hist) .* q_probe * dt;
    if N_heat == 0
        p_inc_linear = G0 * q_probe * dt;
    else
        p_inc_linear = [];
    end
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

function [F1_mJcm2, F5_mJcm2] = F_thresholds(mu_m1, G0, dGdT, rho, Cp, Dt, tau_th)
    eta1 = 0.01;
    eta5 = 0.05;
    fac = (G0 * rho * Cp) / (dGdT * mu_m1) * exp(Dt / tau_th);
    F1_mJcm2 = (eta1 * fac) / 10;
    F5_mJcm2 = (eta5 * fac) / 10;
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
