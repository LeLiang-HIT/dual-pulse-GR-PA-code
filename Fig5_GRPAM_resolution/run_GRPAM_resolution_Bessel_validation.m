%% GR-PAM resolution and Bessel sidelobe validation
% This script supports the application-validation part of the manuscript:
%   1) GR-PAM lateral resolution enhancement.
%   2) GR-PAM axial optical sectioning.
%   3) Noise-limited effective resolution for three representative cases.
%   4) GR-Bessel sidelobe suppression under the same three-case logic.
%
% The script is intentionally kept compact for GitHub release. Exploratory
% phantom, RBC-style contrast, and parameter-map sections from the working
% script have been removed.

clear; clc; close all;
rng(4);

%% Output settings
scriptDir = fileparts(mfilename('fullpath'));
outDir = fullfile(scriptDir, 'outputs');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
export_results = true;

%% Physical parameters
Physics.rho      = 1000;        % kg/m^3
Physics.Cp       = 3600;        % J/(kg*K)
Physics.G0       = 0.13;        % baseline Grueneisen parameter
Physics.dGdT     = 0.01;        % 1/K
Physics.alpha_th = 1.4e-7;      % m^2/s
Physics.lambda   = 532e-9;      % m
Physics.eta_th   = 1.0;

%% OR-PAM baseline and thermal-memory scale
P.FWHM_lat_or_um = 5.0;         % OR-PAM lateral FWHM
P.d_um           = 5.0;         % heated-region size used for tau_th

tau_th_us = (P.d_um * 1e-6)^2 / Physics.alpha_th * 1e6;

%% Axial baseline and optical-gating limit
Axial.FWHM_acous_um = 45.0;     % acoustic axial FWHM
Axial.NA            = 0.63;

w0_m = 0.61 * Physics.lambda / Axial.NA;
zR_m = pi * w0_m^2 / Physics.lambda;
Axial.zR_um = zR_m * 1e6;
Axial.FWHM_opt_um = 2 * Axial.zR_um * sqrt(sqrt(2) - 1);
Gain_axial_max = Axial.FWHM_acous_um / Axial.FWHM_opt_um;


%% Three representative operating cases
caseNames = {'Global best', 'Thermal mismatch', 'Weak nonlinearity'};
F_cases_mJcm2 = [10.0, 10.0, 3.0];
mu_cases_cm1  = [400, 400, 400];
beta_target   = [0.03, 1.00, 0.03];

nCases = numel(caseNames);
beta_list  = zeros(1, nCases);
eta_list   = zeros(1, nCases);
dt_list_us = zeros(1, nCases);

for k = 1:nCases
    dt_list_us(k) = beta_target(k) * tau_th_us;
    [beta_list(k), eta_list(k)] = calc_beta_eta( ...
        F_cases_mJcm2(k), mu_cases_cm1(k), P.d_um, dt_list_us(k), Physics);
end

%% Spatial axes and ideal PSFs
x_range = linspace(-3 * P.FWHM_lat_or_um, 3 * P.FWHM_lat_or_um, 600);
z_range = linspace(-80, 80, 600);

[Lat_OR_all, Lat_GR_ideal, LatInfo] = calc_lateral_psf(x_range, P);
[Ax_OR_all,  Ax_GR_ideal,  AxInfo]  = calc_axial_psf(z_range, Axial);

Lat_GR_cases = repmat({Lat_GR_ideal}, 1, nCases);
Ax_GR_cases  = repmat({Ax_GR_ideal},  1, nCases);

%% Noise model and effective-resolution calculations
signal_scale_GR = F_cases_mJcm2 .* mu_cases_cm1 .* eta_list;
signal_scale_GR = signal_scale_GR / signal_scale_GR(1);

signal_scale_OR = F_cases_mJcm2 .* mu_cases_cm1;
signal_scale_OR = signal_scale_OR / signal_scale_OR(1);

SNR_peak_case1_GR = 40;
sigma_GR = 1 / SNR_peak_case1_GR;
sigma_OR = sigma_GR / sqrt(2);
SNR_min_det = 6;

SNR_GR_peak = signal_scale_GR / sigma_GR;
SNR_OR_peak = signal_scale_OR / sigma_OR;

dsep_grid_um = linspace(0.3 * P.FWHM_lat_or_um, 2.5 * P.FWHM_lat_or_um, 60);
dz_grid_or_um = linspace(0.3 * Axial.FWHM_acous_um, 2.5 * Axial.FWHM_acous_um, 80);
dz_grid_gr_um = linspace(0.3 * Axial.FWHM_opt_um, 2.5 * Axial.FWHM_acous_um, 200);

Lat_eff_GR_um = nan(1, nCases);
Lat_eff_OR_um = nan(1, nCases);
Ax_eff_GR_um  = nan(1, nCases);
Ax_eff_OR_um  = nan(1, nCases);

for k = 1:nCases
    Lat_eff_GR_um(k) = eval_lateral_eff( ...
        Lat_GR_cases{k}, signal_scale_GR(k), sigma_GR, dsep_grid_um, x_range, SNR_min_det);
    Lat_eff_OR_um(k) = eval_lateral_eff( ...
        Lat_OR_all, signal_scale_OR(k), sigma_OR, dsep_grid_um, x_range, SNR_min_det);

    Ax_eff_GR_um(k) = eval_axial_eff( ...
        Ax_GR_cases{k}, signal_scale_GR(k), sigma_GR, dz_grid_gr_um, z_range, SNR_min_det);
    Ax_eff_OR_um(k) = eval_axial_eff( ...
        Ax_OR_all, signal_scale_OR(k), sigma_OR, dz_grid_or_um, z_range, SNR_min_det);
end


%% Figure 1: GR-PAM lateral and axial validation
fig1 = figure('Color', 'w', 'Name', 'GR-PAM resolution validation', ...
    'Position', [100, 100, 1180, 780]);
tiledlayout(fig1, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

colors = [0.85 0.33 0.10; 0.00 0.45 0.74; 0.47 0.67 0.19];

ax1 = nexttile; hold(ax1, 'on'); box(ax1, 'on');
plot(ax1, x_range, Lat_OR_all, 'k--', 'LineWidth', 1.6, 'DisplayName', 'OR-PAM');
plot(ax1, x_range, Lat_GR_ideal, 'r-', 'LineWidth', 2.0, 'DisplayName', 'GR-PAM');
xlabel(ax1, 'Lateral x [um]');
ylabel(ax1, 'Normalized amplitude');
title(ax1, sprintf('Lateral PSF: %.2f -> %.2f um', LatInfo.fwhm_lin, LatInfo.fwhm_gr));
legend(ax1, 'Location', 'best');
axis_style(ax1);

ax2 = nexttile; hold(ax2, 'on'); box(ax2, 'on');
fill(ax2, [x_range fliplr(x_range)], ...
    [sigma_GR * ones(size(x_range)) zeros(size(x_range))], ...
    [0.88 0.88 0.88], 'EdgeColor', 'none', 'DisplayName', 'Noise floor');
for k = 1:nCases
    prof = signal_scale_GR(k) * Lat_GR_cases{k};
    plot(ax2, x_range, prof, '-', 'Color', colors(k,:), 'LineWidth', 1.8, ...
        'DisplayName', caseNames{k});
end
xlabel(ax2, 'Lateral x [um]');
ylabel(ax2, 'Signal amplitude');
title(ax2, 'Lateral GR signal with noise floor');
legend(ax2, 'Location', 'best');
axis_style(ax2);

ax3 = nexttile; hold(ax3, 'on'); box(ax3, 'on');
plot(ax3, z_range, Ax_OR_all, 'b--', 'LineWidth', 1.6, 'DisplayName', 'OR-PAM');
plot(ax3, z_range, Ax_GR_ideal, 'r-', 'LineWidth', 2.0, 'DisplayName', 'GR-PAM');
xlabel(ax3, 'Axial z [um]');
ylabel(ax3, 'Normalized amplitude');
title(ax3, sprintf('Axial response: %.1f -> %.2f um', AxInfo.fwhm_or, AxInfo.fwhm_gr));
legend(ax3, 'Location', 'best');
axis_style(ax3);

ax4 = nexttile; hold(ax4, 'on'); box(ax4, 'on');
fill(ax4, [z_range fliplr(z_range)], ...
    [sigma_GR * ones(size(z_range)) zeros(size(z_range))], ...
    [0.88 0.88 0.88], 'EdgeColor', 'none', 'DisplayName', 'Noise floor');
for k = 1:nCases
    prof_z = signal_scale_GR(k) * Ax_GR_cases{k};
    plot(ax4, z_range, prof_z, '-', 'Color', colors(k,:), 'LineWidth', 1.8, ...
        'DisplayName', caseNames{k});
end
xlabel(ax4, 'Axial z [um]');
ylabel(ax4, 'Signal amplitude');
title(ax4, 'Axial GR signal with noise floor');
legend(ax4, 'Location', 'best');
axis_style(ax4);

%% Figure 2: effective resolution under noise
fig2 = figure('Color', 'w', 'Name', 'Noise-limited effective resolution', ...
    'Position', [150, 150, 1050, 420]);
tiledlayout(fig2, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

caseLabels = {'Global', 'Mismatch', 'Weak'};

ax5 = nexttile; hold(ax5, 'on'); box(ax5, 'on');
bar(ax5, [Lat_eff_OR_um(:), Lat_eff_GR_um(:)], 'grouped');
set(ax5, 'XTick', 1:nCases, 'XTickLabel', caseLabels);
ylabel(ax5, 'Effective lateral resolution [um]');
title(ax5, sprintf('Lateral effective resolution, SNR_{min}=%.0f', SNR_min_det));
legend(ax5, {'OR-PAM', 'GR-PAM'}, 'Location', 'northwest');
axis_style(ax5);

ax6 = nexttile; hold(ax6, 'on'); box(ax6, 'on');
bar(ax6, [Ax_eff_OR_um(:), Ax_eff_GR_um(:)], 'grouped');
set(ax6, 'XTick', 1:nCases, 'XTickLabel', caseLabels);
ylabel(ax6, 'Effective axial resolution [um]');
title(ax6, sprintf('Axial effective resolution, SNR_{min}=%.0f', SNR_min_det));
legend(ax6, {'OR-PAM', 'GR-PAM'}, 'Location', 'northwest');
axis_style(ax6);

%% Bessel sidelobe-suppression validation
% This section generates Figure 3 and saves GR_Bessel_sidelobe_suppression.png.
Bessel.Rmax_um        = 50;
Bessel.Nr             = 2000;
Bessel.kr_m           = 2.5e6;
Bessel.mu_a_cm1       = 200;
Bessel.FWHM_acous_um  = 30;

BesselNoise.SNR_ref_main = 40;
[r_um, psf_GR_case1_clean, psf_lin_case1_clean] = build_bessel_gr_psf( ...
    F_cases_mJcm2(1), beta_list(1), P.d_um, Bessel, Physics);

amp_scale_ref = (Bessel.mu_a_cm1^2) * (F_cases_mJcm2(1)^2) * exp(-beta_list(1));
BesselNoise.sigma = (amp_scale_ref * max(psf_GR_case1_clean)) / BesselNoise.SNR_ref_main;

Metrics = repmat(struct('FWHM_um', NaN, 'PSLR_clean_dB', NaN, 'PSLR_eff_dB', NaN, ...
    'SNR_main', NaN, 'SNR_side', NaN, 'DeltaT_max', NaN), 1, nCases);

for k = 1:nCases
    BCase.F_mJcm2 = F_cases_mJcm2(k);
    BCase.beta = beta_list(k);
    BCase.d_heat_um = P.d_um;
    Metrics(k) = eval_bessel_gr_case(BCase, Bessel, Physics, BesselNoise);
end


%% Figure 3: Bessel sidelobe-suppression images
% This figure shows the Bessel-beam image panels used for the sidelobe-suppression validation.
imgAxis_um = linspace(-48, 48, 420);
dimg_um = imgAxis_um(2) - imgAxis_um(1);
[Xb, Yb] = meshgrid(imgAxis_um, imgAxis_um);
Rb = hypot(Xb, Yb);

besselFwhm_um = 6.5;
kB = 1.126 / (besselFwhm_um / 2);
finiteAperture = exp(-(Rb / 42).^8);
BesselLinear2D = besselj(0, kB * Rb).^2 .* finiteAperture;
BesselLinear2D = BesselLinear2D / max(BesselLinear2D(:));
BesselLinear2D = blur2d(BesselLinear2D, 0.45, dimg_um);
BesselLinear2D = BesselLinear2D / max(BesselLinear2D(:));

sigma_img_OR = sigma_OR / sqrt(4);
sigma_img_GR = sigma_GR / sqrt(16);
linearBesselImage = max(BesselLinear2D + sigma_img_OR * randn(size(BesselLinear2D)), 0);
linearBesselSLR_dB = sidelobe_rejection_db(imgAxis_um, BesselLinear2D, besselFwhm_um);

BesselImages = cell(1, nCases);
BesselShapes = cell(1, nCases);
BesselImageSLR_dB = nan(1, nCases);
alpha_um2_us = Physics.alpha_th * 1e6;

for k = 1:nCases
    sigma_th_um = sqrt(2 * alpha_um2_us * dt_list_us(k));
    heatBessel = blur2d(BesselLinear2D, sigma_th_um, dimg_um);
    shapeBessel = BesselLinear2D .* heatBessel;
    shapeBessel = shapeBessel / max(shapeBessel(:));

    BesselShapes{k} = shapeBessel;
    BesselImages{k} = max(signal_scale_GR(k) * shapeBessel + ...
        sigma_img_GR * randn(size(shapeBessel)), 0);

    if SNR_GR_peak(k) >= SNR_min_det
        BesselImageSLR_dB(k) = sidelobe_rejection_db(imgAxis_um, shapeBessel, besselFwhm_um);
    end
end

fig3 = figure('Color', 'w', 'Name', 'GR-Bessel sidelobe suppression', ...
    'Position', [180, 120, 1200, 900]);
tiledlayout(fig3, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

ax7 = nexttile;
show_pam_image(ax7, imgAxis_um, linearBesselImage, pam_colormap(256), [0 1], ...
    sprintf('Linear Bessel  |  SLR %.1f dB', linearBesselSLR_dB));

ax8 = nexttile;
show_pam_image(ax8, imgAxis_um, BesselImages{1}, pam_colormap(256), [0 1], ...
    format_bessel_label(caseNames{1}, BesselImageSLR_dB(1)));

ax9 = nexttile;
show_pam_image(ax9, imgAxis_um, BesselImages{2}, pam_colormap(256), [0 1], ...
    format_bessel_label(caseNames{2}, BesselImageSLR_dB(2)));

ax10 = nexttile;
show_pam_image(ax10, imgAxis_um, BesselImages{3}, pam_colormap(256), [0 1], ...
    format_bessel_label(caseNames{3}, BesselImageSLR_dB(3)));

%% Save outputs
caseName = string(caseNames(:));
F_mJcm2 = F_cases_mJcm2(:);
mu_cm1 = mu_cases_cm1(:);
beta = beta_list(:);
dt_us = dt_list_us(:);
eta = eta_list(:);
SNR_GR = SNR_GR_peak(:);
SNR_OR = SNR_OR_peak(:);
Lat_eff_OR = Lat_eff_OR_um(:);
Lat_eff_GR = Lat_eff_GR_um(:);
Ax_eff_OR = Ax_eff_OR_um(:);
Ax_eff_GR = Ax_eff_GR_um(:);
Bessel_FWHM_um = [Metrics.FWHM_um].';
Bessel_PSLR_eff_dB = [Metrics.PSLR_eff_dB].';
Bessel_SNR_main = [Metrics.SNR_main].';
Bessel_DeltaT_max = [Metrics.DeltaT_max].';
Bessel_Image_SLR_dB = BesselImageSLR_dB(:);

metricsTable = table(caseName, F_mJcm2, mu_cm1, beta, dt_us, eta, ...
    SNR_GR, SNR_OR, Lat_eff_OR, Lat_eff_GR, Ax_eff_OR, Ax_eff_GR, ...
    Bessel_FWHM_um, Bessel_PSLR_eff_dB, Bessel_SNR_main, Bessel_DeltaT_max, Bessel_Image_SLR_dB);

if export_results
    exportgraphics(fig1, fullfile(outDir, 'GRPAM_resolution_profiles.png'), ...
        'Resolution', 300, 'BackgroundColor', 'white');
    exportgraphics(fig2, fullfile(outDir, 'GRPAM_effective_resolution.png'), ...
        'Resolution', 300, 'BackgroundColor', 'white');
    exportgraphics(fig3, fullfile(outDir, 'GR_Bessel_sidelobe_suppression.png'), ...
        'Resolution', 300, 'BackgroundColor', 'white');

    savefig(fig1, fullfile(outDir, 'GRPAM_resolution_profiles.fig'));
    savefig(fig2, fullfile(outDir, 'GRPAM_effective_resolution.fig'));
    savefig(fig3, fullfile(outDir, 'GR_Bessel_sidelobe_suppression.fig'));

    writetable(metricsTable, fullfile(outDir, 'GRPAM_Bessel_metrics.csv'));
    save(fullfile(outDir, 'GRPAM_Bessel_source_data.mat'), ...
        'Physics', 'P', 'Axial', 'caseNames', 'F_cases_mJcm2', 'mu_cases_cm1', ...
        'beta_list', 'eta_list', 'dt_list_us', 'tau_th_us', ...
        'x_range', 'z_range', 'Lat_OR_all', 'Lat_GR_cases', 'Ax_OR_all', 'Ax_GR_cases', ...
        'LatInfo', 'AxInfo', 'signal_scale_GR', 'signal_scale_OR', ...
        'SNR_GR_peak', 'SNR_OR_peak', 'sigma_GR', 'sigma_OR', 'SNR_min_det', ...
        'Lat_eff_OR_um', 'Lat_eff_GR_um', 'Ax_eff_OR_um', 'Ax_eff_GR_um', ...
        'Bessel', 'BesselNoise', 'r_um', 'psf_lin_case1_clean', 'Metrics', ...
        'imgAxis_um', 'BesselLinear2D', 'linearBesselImage', 'linearBesselSLR_dB', ...
        'BesselImages', 'BesselShapes', 'BesselImageSLR_dB', 'metricsTable');
end


%% Local functions
function [Lin, GR, Info] = calc_lateral_psf(x_um, P)
    k = 4 * log(2);
    I = exp(-k * (x_um.^2) / (P.FWHM_lat_or_um^2));
    I = I / max(I);

    Lin = I;
    GR = I.^2;
    GR = GR / max(GR);

    Info.fwhm_lin = calc_fwhm_numeric(x_um, Lin);
    Info.fwhm_gr = calc_fwhm_numeric(x_um, GR);
    Info.fwhm_ideal = P.FWHM_lat_or_um / sqrt(2);
    Info.gain = Info.fwhm_lin / Info.fwhm_gr;
    Info.gain_ideal = Info.fwhm_lin / Info.fwhm_ideal;
end

function [Ax_OR, Ax_GR, Info] = calc_axial_psf(z_um, Axial)
    k = 4 * log(2);
    Ax_OR = exp(-k * (z_um.^2) / (Axial.FWHM_acous_um^2));
    Ax_OR = Ax_OR / max(Ax_OR);

    zR_um = Axial.zR_um;
    Ax_GR = 1 ./ (1 + (z_um ./ zR_um).^2).^2;
    Ax_GR = Ax_GR / max(Ax_GR);

    Info.fwhm_or = calc_fwhm_numeric(z_um, Ax_OR);
    Info.fwhm_gr = calc_fwhm_numeric(z_um, Ax_GR);
    Info.fwhm_opt = Axial.FWHM_opt_um;
    Info.gain = Info.fwhm_or / Info.fwhm_gr;
    Info.gain_max = Info.fwhm_or / Info.fwhm_opt;
end

function f = calc_fwhm_numeric(x, y)
    y = y / max(y);
    idx = find(y >= 0.5);
    if numel(idx) < 2
        f = NaN;
    else
        f = x(idx(end)) - x(idx(1));
    end
end

function [beta, eta] = calc_beta_eta(F_mJcm2, mu_a_cm1, d_um, dt_us, Phys)
    F_J_m2 = F_mJcm2 * 1e-3 / 1e-4;
    mu_a_m1 = mu_a_cm1 * 100;
    d_m = d_um * 1e-6;
    dt_s = dt_us * 1e-6;

    tau_th = d_m^2 / Phys.alpha_th;
    beta = dt_s / tau_th;

    DeltaT0 = Phys.eta_th * mu_a_m1 * F_J_m2 / (Phys.rho * Phys.Cp);
    DeltaT = DeltaT0 * exp(-beta);
    PA1 = Phys.G0 * Phys.eta_th * mu_a_m1 * F_J_m2;
    DeltaPA = Phys.dGdT * DeltaT * Phys.eta_th * mu_a_m1 * F_J_m2;
    eta = DeltaPA / PA1;
end

function d_eff = eval_lateral_eff(psf_norm, Ak, sigma_n, dsep_grid_um, x_range, SNR_min_det)
    psf_norm = psf_norm / max(psf_norm);
    d_eff = NaN;
    for id = 1:numel(dsep_grid_um)
        dsep = dsep_grid_um(id);
        half = dsep / 2;
        psf1 = interp1(x_range, psf_norm, x_range - half, 'linear', 0);
        psf2 = interp1(x_range, psf_norm, x_range + half, 'linear', 0);
        prof = Ak * (psf1 + psf2);
        centerMask = abs(x_range) < (dsep / 4);
        leftMask = x_range < -(dsep / 4);
        rightMask = x_range > (dsep / 4);
        if ~any(centerMask) || ~any(leftMask) || ~any(rightMask)
            continue;
        end
        leftPeak = max(prof(leftMask));
        rightPeak = max(prof(rightMask));
        peakVal = min(leftPeak, rightPeak);
        valleyVal = min(prof(centerMask));
        peakSNR = peakVal / sigma_n;
        contrastSNR = (peakVal - valleyVal) / sigma_n;
        if peakSNR >= SNR_min_det && contrastSNR >= SNR_min_det
            d_eff = dsep;
            break;
        end
    end
end

function d_eff = eval_axial_eff(psf_norm, Ak, sigma_n, dsep_grid_um, z_range, SNR_min_det)
    psf_norm = psf_norm / max(psf_norm);
    d_eff = NaN;
    for id = 1:numel(dsep_grid_um)
        dsep = dsep_grid_um(id);
        half = dsep / 2;
        psf1 = interp1(z_range, psf_norm, z_range - half, 'linear', 0);
        psf2 = interp1(z_range, psf_norm, z_range + half, 'linear', 0);
        prof = Ak * (psf1 + psf2);
        centerMask = abs(z_range) < (dsep / 4);
        leftMask = z_range < -(dsep / 4);
        rightMask = z_range > (dsep / 4);
        if ~any(centerMask) || ~any(leftMask) || ~any(rightMask)
            continue;
        end
        leftPeak = max(prof(leftMask));
        rightPeak = max(prof(rightMask));
        peakVal = min(leftPeak, rightPeak);
        valleyVal = min(prof(centerMask));
        peakSNR = peakVal / sigma_n;
        contrastSNR = (peakVal - valleyVal) / sigma_n;
        if peakSNR >= SNR_min_det && contrastSNR >= SNR_min_det
            d_eff = dsep;
            break;
        end
    end
end

function [r_um, psf_GR_clean, psf_lin_clean] = build_bessel_gr_psf(F_mJcm2, beta, d_heat_um, Bessel, Phys) %#ok<INUSD>
    N = Bessel.Nr;
    R = Bessel.Rmax_um;
    r_full_um = linspace(-R, R, 2 * N - 1);
    r_full_m = r_full_um * 1e-6;
    idx0 = N;
    r_um = r_full_um(idx0:end);

    I_full = besselj(0, Bessel.kr_m * abs(r_full_m)).^2;
    I_full = I_full / max(I_full);

    k_gauss = 4 * log(2);
    h_acous_full = exp(-k_gauss * (r_full_um.^2) / (Bessel.FWHM_acous_um^2));
    h_acous_full = h_acous_full / sum(h_acous_full);

    psf_lin_full = conv(I_full, h_acous_full, 'same');
    psf_lin_full = psf_lin_full / max(psf_lin_full);
    psf_lin_clean = psf_lin_full(idx0:end);
    psf_lin_clean = psf_lin_clean / max(psf_lin_clean);

    beta_clamped = max(beta, 0);
    FWHM_th = d_heat_um * sqrt(1 + beta_clamped);
    G_th_full = exp(-k_gauss * (r_full_um.^2) / (FWHM_th^2));
    G_th_full = G_th_full / sum(G_th_full);

    DeltaT_full = conv(I_full, G_th_full, 'same');
    DeltaT_full = DeltaT_full / max(DeltaT_full);

    src_GR_full = I_full .* DeltaT_full;
    src_GR_full = src_GR_full / max(src_GR_full);

    psf_GR_full = conv(src_GR_full, h_acous_full, 'same');
    psf_GR_full = psf_GR_full / max(psf_GR_full);
    psf_GR_clean = psf_GR_full(idx0:end);
    psf_GR_clean = psf_GR_clean / max(psf_GR_clean);
end

function Metrics = eval_bessel_gr_case(Case, Bessel, Phys, Noise)
    [~, psf_GR_clean, ~] = build_bessel_gr_psf( ...
        Case.F_mJcm2, Case.beta, Case.d_heat_um, Bessel, Phys);

    Metrics.FWHM_um = calc_fwhm_numeric(linspace(0, Bessel.Rmax_um, Bessel.Nr), psf_GR_clean);
    [A_main, A_side] = find_main_sidelobe(psf_GR_clean);

    amp_scale = (Bessel.mu_a_cm1^2) * (Case.F_mJcm2^2) * exp(-Case.beta);
    sigma_n = Noise.sigma;

    Metrics.SNR_main = (amp_scale * A_main) / sigma_n;
    Metrics.SNR_side = (amp_scale * A_side) / sigma_n;

    if A_side <= 0
        Metrics.PSLR_clean_dB = Inf;
    else
        Metrics.PSLR_clean_dB = 20 * log10(A_main / A_side);
    end

    noise_rel = sigma_n / max(amp_scale * A_main, eps);
    A_side_eff = max(A_side, noise_rel);
    Metrics.PSLR_eff_dB = 20 * log10(A_main / A_side_eff);

    Metrics.DeltaT_max = calc_two_pulse_deltaT_max( ...
        Case.F_mJcm2, Bessel.mu_a_cm1, Case.d_heat_um, Case.beta, Phys);
end

function [A_main, A_side] = find_main_sidelobe(psf)
    psf = psf(:).';
    A_main = psf(1);
    if numel(psf) <= 3
        A_side = 0;
        return;
    end
    y = psf(2:end);
    peakMask = [false, y(2:end-1) > y(1:end-2) & y(2:end-1) > y(3:end), false];
    pks = y(peakMask);
    if isempty(pks)
        A_side = 0;
    else
        A_side = max(pks);
    end
end

function DeltaT_max = calc_two_pulse_deltaT_max(F_mJcm2, mu_a_cm1, d_heat_um, beta, Phys)
    F_J_m2 = F_mJcm2 * 1e-3 / 1e-4;
    mu_a_m1 = mu_a_cm1 * 100;
    DeltaT0 = Phys.eta_th * mu_a_m1 * F_J_m2 / (Phys.rho * Phys.Cp);
    DeltaTres = DeltaT0 * exp(-beta);
    DeltaT_max = DeltaT0 + DeltaTres;
end

function axis_style(ax)
    set(ax, 'FontName', 'Arial', 'FontSize', 11, 'LineWidth', 1.0, ...
        'TickDir', 'out', 'Box', 'on');
    grid(ax, 'off');
end

function y = blur2d(x, sigma_um, dx_um)
    if sigma_um <= 0.05
        y = x;
        return;
    end
    halfWidth = max(3, ceil(5 * sigma_um / dx_um));
    ax = (-halfWidth:halfWidth) * dx_um;
    ker = exp(-0.5 * (ax / sigma_um).^2);
    ker = ker / sum(ker);
    y = conv2(conv2(x, ker, 'same'), ker.', 'same');
end

function db = sidelobe_rejection_db(axis_um, profile, fwhm_um)
    if ~isvector(profile)
        profile = profile(round(size(profile, 1) / 2), :);
    end
    p = profile(:).' / max(profile(:) + eps);
    mainMask = abs(axis_um) <= fwhm_um / 2;
    sideMask = abs(axis_um) >= 1.15 * fwhm_um & abs(axis_um) <= 28;
    mainPeak = max(p(mainMask));
    sidePeak = max(p(sideMask));
    db = 20 * log10(mainPeak / max(sidePeak, eps));
end

function cmap = pam_colormap(n)
    if nargin < 1
        n = 256;
    end
    anchor = [
        0.02 0.02 0.20
        0.10 0.08 0.55
        0.05 0.28 0.90
        0.00 0.78 1.00
        0.95 0.95 0.05
        1.00 1.00 0.92
    ];
    x = linspace(0, 1, size(anchor, 1));
    xi = linspace(0, 1, n);
    cmap = interp1(x, anchor, xi, 'pchip');
    cmap = max(min(cmap, 1), 0);
end

function show_pam_image(ax, axis_um, img, cmap, climVals, labelText)
    imagesc(ax, axis_um, axis_um, img);
    axis(ax, 'image');
    axis(ax, 'off');
    colormap(ax, cmap);
    caxis(ax, climVals);
    text(ax, 0.02, 1.025, labelText, ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'bottom', ...
        'FontName', 'Times New Roman', ...
        'FontWeight', 'bold', ...
        'FontSize', 16, ...
        'Color', 'k', ...
        'BackgroundColor', 'w', ...
        'Margin', 4);
end

function labelText = format_bessel_label(caseName, slrValue)
    if isnan(slrValue)
        labelText = sprintf('%s  |  below threshold', caseName);
    else
        labelText = sprintf('%s  |  SLR %.1f dB', caseName, slrValue);
    end
end

