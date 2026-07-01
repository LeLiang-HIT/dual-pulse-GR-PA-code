%% ========================================================================
%  Manuscript figure generation script
%
%  Project: Dual-pulse Gruneisen-relaxation nonlinear photoacoustics
%  Purpose: Generate manuscript-supporting figures and source data for Fig.2--Fig.4.
%
%  Notes for reviewers/users:
%  1. Place Best_Impulse_Response_Avg.mat in the same folder as this script.
%  2. Run this script from the MATLAB current folder.
%  3. Generated PNG/FIG figures and source-data MAT files are saved to ./outputs/.
%  4. Exploratory/debugging figure blocks were removed; the core calculation flow is preserved.
%
%  Removed internal figure blocks from the original working script:
%  Fig.10, Fig.12, Fig.17, Fig.18, Fig.23, Fig.24, Fig.25, and Fig.31.
% ========================================================================
clear; clc; close all;
clear col

%% Output folder
output_dir = fullfile(pwd, 'outputs');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Keep the output folder clean: this review package does not export PDF files.
% Any PDF files left from previous runs are removed at startup.
existing_pdf_files = dir(fullfile(output_dir, '*.pdf'));
for ii_pdf = 1:numel(existing_pdf_files)
    delete(fullfile(output_dir, existing_pdf_files(ii_pdf).name));
end
clear existing_pdf_files ii_pdf;

%% Global display settings
set(groot,'defaultTextInterpreter','tex'); 
set(groot,'defaultAxesTickLabelInterpreter','tex'); 
set(groot,'defaultLegendInterpreter','tex'); 
set(groot,'defaultAxesFontName','Arial'); 
set(groot,'defaultTextFontName','Arial'); 
set(groot,'defaultAxesFontSize',11); 
set(groot,'defaultLineLineWidth',1.2);
% Line colors
col.P1    = [0 114 178]/255;   % Pulse1
col.P2    = [213 94 0]/255;    % Pulse2
col.DELTA = [123 97 255]/255;  % ΔPA
col.GRAY  = [0.30 0.30 0.30];
% Regime shading colors
cLIN  = [0.83 0.83 0.83];
cWNL  = [1.00 0.95 0.75];      % 1–5%
cNL   = [1.00 0.85 0.70];      % ≥5%
cSAFE = [0.80 1.00 0.80];
%% Baseline parameters
t1          = 5e-6;          % first pulse time [s]
Delta_t0    = 40e-6;         % inter-pulse delay [s]
tauL_ns0    = 10;            % pulse-width FWHM [ns]
F_mJcm2_0   = 5;             % fluence per pulse [mJ/cm^2]
mu_a_cm0    = 240;           % absorption coefficient [1/cm]
w0_um0      = 5;             % characteristic heating dimension [um]
eta_th      = 1;           % photothermal conversion efficiency
Gamma0      = 0.12;          % baseline Gruneisen parameter
dGamma_dT   = 0.01;          % dΓ/dT [1/K]
c0   = 1500;  
rho = 1000;  
Cp = 3600;
alpha_coeff = 0.75;  
alpha_power = 1.5;
alpha_th    = 1.3e-7;        % thermal diffusivity [m^2/s]
% Parameter scan ranges
scanF   = linspace(0,25,10);                       % [mJ/cm^2]
scanMu  = linspace(10,1000,10);                   % [1/cm]
scanTau = linspace(10,40,10);                     % [ns]
scanDt  = linspace(10e-6,200e-6,10);               % [s]
scanW0  = linspace(1,20,10);                       % [um]
MPE_mJcm2 = 20;                                    % safety threshold


if exist('Best_Impulse_Response_Avg.mat', 'file')
    H_data = load('Best_Impulse_Response_Avg.mat');
    h_raw  = double(H_data.h_t_final(:).'); 
else
    error('Best_Impulse_Response_Avg.mat was not found. Place the empirical impulse-response MAT file in the current folder.');
end

fs_exp = 250e6;       
dt_exp = 1/fs_exp;    % 4 ns

dt = 1e-9;            

t_exp_axis = (0:length(h_raw)-1) * dt_exp;
t_sim_axis = 0:dt:t_exp_axis(end);

ref_wave = interp1(t_exp_axis, h_raw, t_sim_axis, 'pchip', 0);

ref_wave = ref_wave / max(abs(ref_wave)); 

[~, kpk] = max(abs(ref_wave));

start_idx = max(1, kpk - round(100e-9/dt)); 
ref_wave = ref_wave(start_idx:end);

if exist('t1','var') && exist('Delta_t0','var') && exist('scanDt','var')
    Delta_t_max = max([Delta_t0, max(scanDt)]);
    T_end = max(t1 + Delta_t_max + 30e-6, 60e-6);
else
    T_end = 80e-6;
end

Nt = ceil(T_end / dt);
t  = (0:Nt-1) * dt;

fprintf('------------------------------------------------------\n');
fprintf('Empirical impulse response loaded successfully.\n');
fprintf('  > Sampling grid: 250 MHz experimental data -> 1 GHz simulation grid.\n');
fprintf('  > Amplitude normalization: enabled; outputs are reported in a.u.\n');
fprintf('  > Laser pulse broadening: 10 ns broadening is neglected in this version.\n');
fprintf('------------------------------------------------------\n');

clear kgrid h_interp H_data;



cm1_to_m1    = @(x) x*100;
mJcm2_to_Jm2 = @(x) x*10;
um_to_m      = @(x) x*1e-6;
pp           = @(x) max(x)-min(x);
gaussI = @(F,t,t0,tauFWHM,dt) ...
    (F * exp(-0.5*((t - t0)/(tauFWHM/(2*sqrt(2*log(2))))).^2) ) / ...
    (sum(exp(-0.5*((t - t0)/(tauFWHM/(2*sqrt(2*log(2))))).^2))*dt);
tau_th_from_w0 = @(w0_m, a_th) w0_m.^2 ./ a_th;



F0  = mJcm2_to_Jm2(F_mJcm2_0);
mu0 = cm1_to_m1(mu_a_cm0);
w00 = um_to_m(w0_um0);
tauL0 = tauL_ns0*1e-9;
tau_th0 = tau_th_from_w0(w00, alpha_th);      % ≈192 µs
[pinc1,pinc2] = build_pincs_prevOnly(t,dt,F0,F0,t1,Delta_t0,tauL0,mu0,eta_th,rho,Cp,Gamma0,dGamma_dT,tau_th0,gaussI);
p1 = conv(pinc1,ref_wave,'full'); 
p2 = conv(pinc2,ref_wave,'full');
t_long = (0:numel(p1)-1)*dt;
fprintf('baseline: tau_th=%.1f us,  PA1=%.2f a.u., PA2=%.2f a.u., ΔPA=%.2f a.u.\n', ...
        tau_th0*1e6, pp(p1), pp(p2), (pp(p2)-pp(p1)));

F_list = mJcm2_to_Jm2(scanF);
PA1_F=zeros(size(F_list)); PA2_F=PA1_F; DP_F=PA1_F;
for i=1:numel(F_list)
    [pinc1i,pinc2i]=build_pincs_prevOnly(t,dt,F_list(i),F_list(i),t1,Delta_t0,tauL0,mu0,eta_th,rho,Cp,Gamma0,dGamma_dT,tau_th0,gaussI);
    p1i=conv(pinc1i,ref_wave,'full'); p2i=conv(pinc2i,ref_wave,'full');
    PA1_F(i)=pp(p1i); PA2_F(i)=pp(p2i); DP_F(i)=PA2_F(i)-PA1_F(i);
end
F_plot = 10; Fp = mJcm2_to_Jm2(F_plot);
[pinc1p,pinc2p]=build_pincs_prevOnly(t,dt,Fp,Fp,t1,Delta_t0,tauL0,mu0,eta_th,rho,Cp,Gamma0,dGamma_dT,tau_th0,gaussI);
p1p=conv(pinc1p,ref_wave,'full'); p2p=conv(pinc2p,ref_wave,'full'); tLp=(0:numel(p1p)-1)*dt;
%ax = nexttile;
figure('Color','w','Name','A1: Example waveform');
ax = axes;
hold(ax,'on');
plot(ax, tLp*1e6, p1p, '--', 'Color', col.P1);
plot(ax, tLp*1e6, p2p, '-',  'Color', col.P2);
grid(ax,'on'); xlim(ax,[0 60]);
xlabel(ax,'Time [\mus]'); ylabel(ax,'Pressure [a.u.]');
%title(ax, sprintf('Synthesized waveforms at F = %.1f mJ/cm^2', F_plot));
pp1_kPa = (max(p1p)-min(p1p));  pp2_kPa = (max(p2p)-min(p2p));
legend(ax, sprintf('PA_1 (pp = %.2f a.u.)', pp1_kPa), sprintf('PA_2 (pp = %.2f a.u.)', pp2_kPa), 'Location','best');
%ax = nexttile;
figure('Color','w','Name','A2: DeltaPA vs F linear');
ax = axes;
plot(ax, scanF, DP_F, 'o-', 'Color', col.DELTA); grid(ax,'on');
xlim(ax,[min(scanF) max(scanF)]); xlabel(ax,'F [mJ/cm^2]'); ylabel(ax,'\DeltaPA [a.u.]');
%title(ax,'\DeltaPA vs F (linear)');
annotate_mpe_noOverlap(ax, MPE_mJcm2, max(scanF));
maskF = DP_F>0; pfit = polyfit(log10(F_list(maskF)), log10(DP_F(maskF)),1);
slopeF = pfit(1);
%ax = nexttile;
figure('Color','w','Name','A3: DeltaPA vs F log-log');
ax = axes;
loglog(ax, F_list(maskF),DP_F(maskF),'o','Color',col.DELTA); hold(ax,'on');
loglog(ax, F_list(maskF), 10.^polyval(pfit,log10(F_list(maskF))),'-','Color',col.GRAY); grid(ax,'on');
xlabel(ax,'F [J/m^2]'); ylabel(ax,'\DeltaPA [a.u.]'); 
%title(ax,sprintf('log–log: slope=%.2f',slopeF));
% ax = nexttile; hold(ax,'on');
% plot(ax, scanF, PA1_F,'-','Color',col.P1,'DisplayName','PA_1');
% plot(ax, scanF, PA2_F,'-','Color',col.P2,'DisplayName','PA_2');
% grid(ax,'on'); xlim(ax,[min(scanF) max(scanF)]);
% xlabel(ax,'F [mJ/cm^2]'); ylabel(ax,'PA [a.u.]');
% title(ax,'PA (pp) vs F with linear/nonlinear zones'); legend(ax,'Location','southeast');
% smart_ylim(ax, [PA1_F; PA2_F], 0.05);
% annotate_linear_regions_F(ax, scanF, mu0, Gamma0, dGamma_dT, rho, Cp, Delta_t0, tau_th0, MPE_mJcm2, cLIN, cWNL, cNL, cSAFE, true);
% A4：PA vs F
%ax = nexttile; 
figure('Color','w','Name','A4: PA vs F with regimes');
ax = axes;
hold(ax,'on');
hP1 = plot(ax, scanF, PA1_F,'-','Color',col.P1,'DisplayName','PA_1');
hP2 = plot(ax, scanF, PA2_F,'-','Color',col.P2,'DisplayName','PA_2');
grid(ax,'on'); xlim(ax,[min(scanF) max(scanF)]);
xlabel(ax,'F [mJ/cm^2]'); ylabel(ax,'PA [a.u.]');
%title(ax,'PA (pp) vs F ');
smart_ylim(ax, [PA1_F; PA2_F], 0.05);
[h1,h5,hMPE] = annotate_linear_regions_F(ax, scanF, mu0, Gamma0, dGamma_dT, rho, Cp, Delta_t0, tau_th0, MPE_mJcm2, cLIN, cWNL, cNL, cSAFE, true);
legend(ax, [hP1 hP2 h1 h5 hMPE], 'Location','southeast');

mu_list_cm = scanMu; mu_list = cm1_to_m1(mu_list_cm);
PA1_mu=zeros(size(mu_list)); PA2_mu=PA1_mu; DP_mu=PA1_mu;
for i=1:numel(mu_list)
    [pinc1i,pinc2i]=build_pincs_prevOnly(t,dt,F0,F0,t1,Delta_t0,tauL0,mu_list(i),eta_th,rho,Cp,Gamma0,dGamma_dT,tau_th0,gaussI);
    p1i=conv(pinc1i,ref_wave,'full'); p2i=conv(pinc2i,ref_wave,'full');
    PA1_mu(i)=pp(p1i); PA2_mu(i)=pp(p2i); DP_mu(i)=PA2_mu(i)-PA1_mu(i);
end
mu_plot_cm = 400; mu_plot = cm1_to_m1(mu_plot_cm);
[pinc1p,pinc2p]=build_pincs_prevOnly(t,dt,F0,F0,t1,Delta_t0,tauL0,mu_plot,eta_th,rho,Cp,Gamma0,dGamma_dT,tau_th0,gaussI);
p1p=conv(pinc1p,ref_wave,'full'); p2p=conv(pinc2p,ref_wave,'full'); tLp=(0:numel(p1p)-1)*dt;
%ax = nexttile;
figure('Color','w','Name','B1: Example waveform');
ax = axes;
hold(ax,'on');
plot(ax, tLp*1e6, p1p, '--', 'Color', col.P1);
plot(ax, tLp*1e6, p2p, '-',  'Color', col.P2);
grid(ax,'on'); xlim(ax,[0 60]);
xlabel(ax,'Time [\mus]'); ylabel(ax,'Pressure [a.u.]');
%title(ax, sprintf('Synthesized waveforms at \\mu_a = %d cm^{-1}', mu_plot_cm));
pp1_kPa = (max(p1p)-min(p1p));  pp2_kPa = (max(p2p)-min(p2p));
legend(ax, sprintf('PA_1 (pp = %.2f a.u.)', pp1_kPa), sprintf('PA_2 (pp = %.2f a.u.)', pp2_kPa), 'Location','best');
% B2：ΔPA vs μa
%ax = nexttile;
figure('Color','w','Name','B2: DeltaPA vs mu_a semilog');
ax = axes;
semilogx(ax, mu_list_cm, DP_mu,'o-','Color',col.DELTA); grid(ax,'on');
xlabel(ax,'\mu_a [1/cm]'); ylabel(ax,'\DeltaPA [a.u.]'); 
%title(ax,'\DeltaPA vs \mu_a');
maskMu = DP_mu>0; pfit = polyfit(log10(mu_list(maskMu)), log10(DP_mu(maskMu)),1);
slopeMu = pfit(1);
%ax = nexttile;
figure('Color','w','Name','B3: DeltaPA vs mu_a log-log');
ax = axes;
loglog(ax, mu_list(maskMu),DP_mu(maskMu),'o','Color',col.DELTA); hold(ax,'on');
loglog(ax, mu_list(maskMu), 10.^polyval(pfit,log10(mu_list(maskMu))),'-','Color',col.GRAY); grid(ax,'on');
xlabel(ax,'\mu_a [1/m]'); ylabel(ax,'\DeltaPA [a.u.]'); 
%title(ax,sprintf('log–log: slope=%.2f',slopeMu));
% plot(ax, mu_list_cm, PA1_mu,'-','Color',col.P1,'DisplayName','PA_1');
% plot(ax, mu_list_cm, PA2_mu,'-','Color',col.P2,'DisplayName','PA_2');
% grid(ax,'on'); xlabel(ax,'\mu_a [1/cm]'); ylabel(ax,'PA [a.u.]');
% title(ax,'PA (pp) vs \mu_a with 3-zone shading'); legend(ax,'Location','southeast');
% smart_ylim(ax, [PA1_mu; PA2_mu], 0.05);
%ax = nexttile; 
figure('Color','w','Name','B4: PA vs mu_a with regimes');
ax = axes;
hold(ax,'on'); xlim(ax,[10 1000]);
hP1 = plot(ax, mu_list_cm, PA1_mu,'-','Color',col.P1,'DisplayName','PA_1');
hP2 = plot(ax, mu_list_cm, PA2_mu,'-','Color',col.P2,'DisplayName','PA_2');
grid(ax,'on'); xlabel(ax,'\mu_a [1/cm]'); ylabel(ax,'PA [a.u.]');
%title(ax,'PA (pp) vs \mu_a ');
smart_ylim(ax, [PA1_mu; PA2_mu], 0.05);
[h1,h5] = annotate_threezone_mu(ax, [10 1000], F0, Gamma0, dGamma_dT, rho, Cp, Delta_t0, tau_th0, cLIN, cWNL, cSAFE, true);
legend(ax, [hP1 hP2 h1 h5], 'Location','southeast');

hold(ax, 'on'); 
cat_colors_bright.uv  = [0.8 0 0.8];
cat_colors_bright.vis = [0.1 0.8 0.1];
cat_colors_bright.nir = [1 0.5 0];
cat_markers.uv  = 'o';
cat_markers.vis = 'p';
cat_markers.nir = 'd';
absorber_data = {
    'DNA/RNA',   80,  'uv';
    'Melanin',   900, 'uv';
    'HbO_2',     250, 'vis';
    'HbR',       350, 'vis';
    'Lipid',     20,  'nir'
};
for i = 1:size(absorber_data, 1)
    name   = absorber_data{i, 1};
    mu_val = absorber_data{i, 2};
    cat    = absorber_data{i, 3};
    
    marker_color = cat_colors_bright.(cat);
    marker_style = cat_markers.(cat);
    
    current_marker_size = 8;
    switch marker_style
        case 'o'
            current_marker_size = 7; 
        case 'p'
            current_marker_size = 11; 
        case 'd'
            current_marker_size = 7; 
    end
    x_limits = get(ax, 'XLim');
    if mu_val >= x_limits(1) && mu_val <= x_limits(2)
        pa2_val = interp1(mu_list_cm, PA2_mu, mu_val);
        
        plot(ax, mu_val, pa2_val, marker_style, ...
            'MarkerFaceColor', marker_color, ... 
            'MarkerEdgeColor', 'none', ...
            'MarkerSize', current_marker_size, ...
            'HandleVisibility', 'off');
        
        text(ax, mu_val + 15, pa2_val, name, ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'Color', 'k', ...
            'BackgroundColor', [1 1 1 0.7], 'EdgeColor', 'none', 'Margin', 2, ...
            'Interpreter', 'tex');
    end
end
hold(ax, 'off');
% =================================================================================
% =================================================================================

tauL_list = scanTau*1e-9;
PA1_tL=zeros(size(tauL_list)); PA2_tL=PA1_tL; DP_tL=PA1_tL;
for i=1:numel(tauL_list)
    [pinc1i,pinc2i]=build_pincs_prevOnly(t,dt,F0,F0,t1,Delta_t0,tauL_list(i),mu0,eta_th,rho,Cp,Gamma0,dGamma_dT,tau_th0,gaussI);
    p1i=conv(pinc1i,ref_wave,'full'); p2i=conv(pinc2i,ref_wave,'full');
    PA1_tL(i)=pp(p1i); PA2_tL(i)=pp(p2i); DP_tL(i)=PA2_tL(i)-PA1_tL(i);
end
Dt_list = scanDt;
PA1_Dt=zeros(size(Dt_list)); PA2_Dt=PA1_Dt; DP_Dt=PA1_Dt;
for i=1:numel(Dt_list)
    [pinc1i,pinc2i]=build_pincs_prevOnly(t,dt,F0,F0,t1,Dt_list(i),tauL0,mu0,eta_th,rho,Cp,Gamma0,dGamma_dT,tau_th0,gaussI);
    p1i=conv(pinc1i,ref_wave,'full'); p2i=conv(pinc2i,ref_wave,'full');
    PA1_Dt(i)=pp(p1i); PA2_Dt(i)=pp(p2i); DP_Dt(i)=PA2_Dt(i)-PA1_Dt(i);
end
w0_list_um = scanW0; w0_list = um_to_m(w0_list_um);
PA1_w0=zeros(size(w0_list)); PA2_w0=PA1_w0; DP_w0=PA1_w0;
for i=1:numel(w0_list)
    tau_th_i = tau_th_from_w0(w0_list(i), alpha_th);
    [pinc1i,pinc2i]=build_pincs_prevOnly(t,dt,F0,F0,t1,Delta_t0,tauL0,mu0,eta_th,rho,Cp,Gamma0,dGamma_dT,tau_th_i,gaussI);
    p1i=conv(pinc1i,ref_wave,'full'); p2i=conv(pinc2i,ref_wave,'full');
    PA1_w0(i)=pp(p1i); PA2_w0(i)=pp(p2i); DP_w0(i)=PA2_w0(i)-PA1_w0(i);
end
figure('Color','w','Name','C: tau_L / Delta_t / d scans');
tiledlayout(3,2,'Padding','compact','TileSpacing','compact');
ax = nexttile; hold(ax,'on');
plot(ax, scanTau, PA1_tL,'--o','Color',col.P1,'MarkerSize',3,'DisplayName','PA_1');
plot(ax, scanTau, PA2_tL,'-^','Color',col.P2,'MarkerSize',3,'DisplayName','PA_2');
grid(ax,'on'); xlabel(ax,'\tau_L [ns]'); ylabel(ax,'PA [a.u.]');
%title(ax,'PA_1 / PA_2 – \tau_L'); legend(ax,'Location','northeast');
smart_ylim(ax, [PA1_tL; PA2_tL], 0.08);
ax = nexttile;
plot(ax, scanTau, DP_tL,'o-','Color',col.DELTA,'MarkerSize',3);
grid(ax,'on'); xlabel(ax,'\tau_L [ns]'); ylabel(ax,'\DeltaPA [a.u.]');
%title(ax,'\DeltaPA – \tau_L');
ax = nexttile; hold(ax,'on');
plot(ax, Dt_list*1e6, PA1_Dt,'--o','Color',col.P1,'MarkerSize',3,'DisplayName','PA_1');
plot(ax, Dt_list*1e6, PA2_Dt,'-^','Color',col.P2,'MarkerSize',3,'DisplayName','PA_2');
grid(ax,'on'); xlabel(ax,'\Delta t [\mus]'); ylabel(ax,'PA [a.u.]');
%title(ax,'PA_1 / PA_2 – \Delta t'); legend(ax,'Location','northeast');
smart_ylim(ax, [PA1_Dt; PA2_Dt], 0.08);
ax = nexttile;
plot(ax, Dt_list*1e6, DP_Dt,'o-','Color',col.DELTA,'MarkerSize',3);
grid(ax,'on'); xlabel(ax,'\Delta t [\mus]'); ylabel(ax,'\DeltaPA [a.u.]');
%title(ax,'\DeltaPA – \Delta t');
ax = nexttile; hold(ax,'on');
plot(ax, w0_list_um, PA1_w0,'--o','Color',col.P1,'MarkerSize',3,'DisplayName','PA_1');
plot(ax, w0_list_um, PA2_w0,'-^','Color',col.P2,'MarkerSize',3,'DisplayName','PA_2');
grid(ax,'on'); xlabel(ax,'w_0 [\mum]'); ylabel(ax,'PA [a.u.]');
%title(ax,'PA_1 / PA_2 – w_0'); legend(ax,'Location','southeast');
smart_ylim(ax, [PA1_w0; PA2_w0], 0.08);
ax = nexttile;
plot(ax, w0_list_um, DP_w0,'o-','Color',col.DELTA,'MarkerSize',3);
grid(ax,'on'); xlabel(ax,'w_0 [\mum]'); ylabel(ax,'\DeltaPA [a.u.]');
%title(ax,'\DeltaPA – w_0');
%% ======================== E1: β-map（Δt × w0） ============================
F_fix_mJ    = 5;
mu_fix_cm   = 240;
tauL_fix_ns = tauL_ns0;
Dt_vec = round(linspace(10e-6, 200e-6, 151) / dt) * dt;
w0_vec = linspace(2, 20, 181);
[W0, DT] = meshgrid(w0_vec, Dt_vec);
DP_map = zeros(size(W0));                   % kPa
for ii = 1:numel(W0)
    DP_map(ii) = evalDP_equalF( ...
        F_fix_mJ, mu_fix_cm, DT(ii), W0(ii), tauL_fix_ns, ...
        t, dt, t1, ref_wave, eta_th, rho, Cp, Gamma0, dGamma_dT, alpha_th);
end
dp_rng = [min(DP_map(:)), max(DP_map(:))];
figure('Color','w','Name','E1a: ΔPA vs (w0, Δt)');
ax1 = axes;
h1 = imagesc(ax1, w0_vec, Dt_vec*1e6, DP_map);
axis(ax1,'xy'); grid(ax1,'on'); caxis(ax1, dp_rng);
set(h1,'AlphaData',~isnan(DP_map));
if isprop(h1,'Interpolation'), set(h1,'Interpolation','bilinear'); end
cb1 = colorbar(ax1); cb1.Label.String = '\DeltaPA [a.u.]';
xlabel(ax1,'w_0 [\mum]'); ylabel(ax1,'\Delta t [\mus]');
%title(ax1, '\DeltaPA vs (w_0, \Delta t)');
beta_levels = [0.03 0.1 0.3 1];
clr = [0 114 178; 230 159 0; 0 158 115; 148 0 211]/255;
hold(ax1,'on');
for k = 1:numel(beta_levels)
    b = beta_levels(k);
    Dt_beta = b * (w0_vec*1e-6).^2 / alpha_th;   % Δt = β τ_th
    plot(ax1, w0_vec, Dt_beta*1e6, '-', 'Color', clr(k,:), 'LineWidth', 2.0, ...
         'DisplayName', sprintf('\\beta = %.2g', b));
end
legend(ax1,'Location','northwest');
hold(ax1, 'off');

F_vec   = linspace(0, 25, 40);     % mJ/cm^2
mu_vec  = linspace(100, 1000, 40); % 1/cm
[MU, FF] = meshgrid(mu_vec, F_vec);
DP_Fmu = zeros(size(FF));          % ΔPA [a.u.]
Delta_t_use = Delta_t0;
w0_use_um   = w0_um0;
for ii = 1:numel(FF)
    DP_Fmu(ii) = evalDP_equalF( ...
        FF(ii), MU(ii), Delta_t_use, w0_use_um, tauL_ns0, ...
        t, dt, t1, ref_wave, eta_th, rho, Cp, Gamma0, dGamma_dT, alpha_th);
end
dp_rng2 = [min(DP_Fmu(:)) max(DP_Fmu(:))];
figure('Color','w','Name','E2: \DeltaPA vs (F,\mu_a)');
ax = axes;
imagesc(ax, mu_vec, F_vec, DP_Fmu); axis(ax,'xy'); grid(ax,'on');
caxis(ax, dp_rng2);
cb = colorbar(ax); cb.Label.String = '\DeltaPA [a.u.]';
xlabel(ax,'\mu_a [1/cm]'); ylabel(ax,'F per pulse [mJ/cm^2]');
%title(ax, ' \DeltaPA vs (F,\mu_a) ');
F1_line = zeros(size(mu_vec));
F5_line = zeros(size(mu_vec));
for i = 1:numel(mu_vec)
    [F1_line(i), F5_line(i)] = F_thresholds( ...
        cm1_to_m1(mu_vec(i)), Gamma0, dGamma_dT, rho, Cp, ...
        Delta_t_use, (w0_use_um*1e-6)^2/alpha_th);
end
hold(ax,'on');
%h1 = plot(ax, mu_vec, F1_line, ':',  'Color',[0.50 0.50 0.50], 'LineWidth',1.6, 'DisplayName','1%'); % ★ 1%
%h2 = plot(ax, mu_vec, F5_line, '--', 'Color',[0.00 0.45 0.85], 'LineWidth',1.8, 'DisplayName','5%');
h1 = plot(ax, mu_vec, F1_line, '-.', 'Color',[0.80 0.00 0.80], 'LineWidth',2.8, 'DisplayName','1%');
h2 = plot(ax, mu_vec, F5_line, '--',  'Color',[0.00 0.60 0.00], 'LineWidth',3.0, 'DisplayName','5%');
legend(ax,[h1 h2],'Location','northwest');
yline(ax, MPE_mJcm2, '--', 'Color',[0.85 0.2 0.2], 'LineWidth',1.6, 'HandleVisibility','off');
yl = ylim(ax);
if yl(2) > MPE_mJcm2
    patch('XData',[mu_vec(1) mu_vec(end) mu_vec(end) mu_vec(1)], ...
          'YData',[MPE_mJcm2 MPE_mJcm2 yl(2) yl(2)], ...
          'FaceColor',[0.5 0.5 0.5], 'FaceAlpha',0.15, ...
          'EdgeColor','none', 'HandleVisibility','off');
    uistack(h1,'top'); uistack(h2,'top');
end
text(ax, mu_vec(end)*0.94, MPE_mJcm2*1.03, 'MPE', ...
     'Color',[0.85 0.2 0.2], 'FontWeight','bold', ...
     'HorizontalAlignment','right','VerticalAlignment','bottom', ...
     'BackgroundColor','w','Margin',2);
legend(ax,[h1 h2],'Location','northwest');

beta_list = [0.03 0.1 0.3 1];
F_line = linspace(0, 20, 20);
fig_E4_final_style = figure('Color','w','Name','E4 Final Style');
ax = axes(fig_E4_final_style);
hold(ax, 'on'); 
box(ax, 'off');
grid(ax, 'on');
plot_handles = gobjects(1, numel(beta_list));
for i = 1:numel(beta_list)
    b = beta_list(i);
    Dt_line = b * ((w0_um0*1e-6)^2 / alpha_th);
    dp = zeros(size(F_line));
    for j = 1:numel(F_line)
        dp(j) = evalDP_equalF(F_line(j), mu_a_cm0, Dt_line, w0_um0, tauL_ns0, ...
                              t, dt, t1, ref_wave, eta_th, rho, Cp, Gamma0, dGamma_dT, alpha_th);
    end
    plot_handles(i) = plot(ax, F_line, dp, 'o-', 'LineWidth', 1.5, ...
        'MarkerFaceColor', 'auto', ...
        'DisplayName', sprintf('\\beta=%.2g', b));
end
xlabel(ax, 'F per pulse [mJ/cm^2]');
ylabel(ax, '\DeltaPA [a.u.]');
legend(plot_handles, 'Location','northwest', 'FontSize', 10);
set(ax, 'FontSize', 11, 'FontName', 'Arial', 'LineWidth', 1);
xlim(ax, [0 20]);
hold(ax, 'off');
%% ===================== E5: Sensitivity（Tornado + SRRC） ==================
dp_base_kPa = evalDP_equalF(F_mJcm2_0, mu_a_cm0, Delta_t0, w0_um0, tauL_ns0, ...
                            t, dt, t1, ref_wave, eta_th, rho, Cp, Gamma0, dGamma_dT, alpha_th);
lofac = 0.8; 
hifac = 1.2;
pars_lo = [mu_a_cm0*lofac, F_mJcm2_0*lofac, w0_um0*lofac, Delta_t0*lofac, tauL_ns0*lofac];
pars_hi = [mu_a_cm0*hifac, F_mJcm2_0*hifac, w0_um0*hifac, Delta_t0*hifac, tauL_ns0*hifac];
dp_lo = zeros(1,5); 
dp_hi = zeros(1,5);
for i = 1:5
    mu_i   = mu_a_cm0; F_i  = F_mJcm2_0; w0_i = w0_um0; Dt_i = Delta_t0; tau_i = tauL_ns0;
    switch i
        case 1, mu_i = pars_lo(1);   % μa
        case 2, F_i  = pars_lo(2);   % F
        case 3, w0_i = pars_lo(3);   % w0
        case 4, Dt_i = pars_lo(4);   % Δt
        case 5, tau_i= pars_lo(5);   % τL
    end
    dp_lo(i) = evalDP_equalF(F_i, mu_i, Dt_i, w0_i, tau_i, ...
                             t, dt, t1, ref_wave, eta_th, rho, Cp, Gamma0, dGamma_dT, alpha_th);
    mu_i   = mu_a_cm0; F_i  = F_mJcm2_0; w0_i = w0_um0; Dt_i = Delta_t0; tau_i = tauL_ns0;
    switch i
        case 1, mu_i = pars_hi(1);
        case 2, F_i  = pars_hi(2);
        case 3, w0_i = pars_hi(3);
        case 4, Dt_i = pars_hi(4);
        case 5, tau_i= pars_hi(5);
    end
    dp_hi(i) = evalDP_equalF(F_i, mu_i, Dt_i, w0_i, tau_i, ...
                             t, dt, t1, ref_wave, eta_th, rho, Cp, Gamma0, dGamma_dT, alpha_th);
end
dp_base = dp_base_kPa;
dp_lo   = dp_lo;
dp_hi   = dp_hi;
rng(1); N = 600;
U = lhsdesign(N,5,'criterion','maximin','iterations',10);  % [0,1] LHS
F_s   = F_mJcm2_0   * (lofac + (hifac-lofac)*U(:,2));
mu_s  = mu_a_cm0    * (lofac + (hifac-lofac)*U(:,1));
w0_s  = w0_um0      * (lofac + (hifac-lofac)*U(:,3));
Dt_s  = Delta_t0    * (lofac + (hifac-lofac)*U(:,4));
tau_s = tauL_ns0    * (lofac + (hifac-lofac)*U(:,5));
Y = zeros(N,1);
for n = 1:N
    Y(n) = evalDP_equalF(F_s(n), mu_s(n), Dt_s(n), w0_s(n), tau_s(n), ...
                         t, dt, t1, ref_wave, eta_th, rho, Cp, Gamma0, dGamma_dT, alpha_th);
end
X = [mu_s F_s w0_s Dt_s tau_s];
Xr = tiedrank(X);   Yr = tiedrank(Y);
Xz = zscore(Xr);    Yz = zscore(Yr);
betaSRRC  = Xz \ Yz;                               % 5×1
[absSRRCs, k] = sort(abs(betaSRRC(:)),'descend');
sgn = sign(betaSRRC(k));
xtags_all = {'\mu_a','F','d','\Deltat','\tau_L'};
names_s   = xtags_all(k);
dp_lo_s = dp_lo(k);
dp_hi_s = dp_hi(k);
figure('Color','w','Name','E5: Sensitivity ranking (Tornado + SRRC)');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
axL = nexttile; hold(axL,'on'); box(axL,'on'); grid(axL,'on');
y  = 1:5;
bw = 16;
offY = 0.22;
cLo = [0.30 0.45 0.85];
cHi = [0.93 0.49 0.18];
for r = 1:5
    xLr = min(dp_lo_s(r), dp_hi_s(r));
    xRr = max(dp_lo_s(r), dp_hi_s(r));
    plot(axL, [xLr xRr], [r r], '-', 'Color',[0.80 0.80 0.80], ...
         'LineWidth', bw, 'HandleVisibility','off');
    plot(axL, dp_lo_s(r), r, 'o', 'MarkerFaceColor', cLo, 'MarkerEdgeColor','w','HandleVisibility','off');
    plot(axL, dp_hi_s(r), r, 'o', 'MarkerFaceColor', cHi, 'MarkerEdgeColor','w','HandleVisibility','off');
    xmid = 0.5*(dp_lo_s(r)+dp_hi_s(r));
    half = 0.5*abs(dp_hi_s(r)-dp_lo_s(r));
    ht = text(axL, xmid, r+offY, sprintf('\\pm %d a.u.', round(half)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'Color',[0.25 0.25 0.25],'FontWeight','bold', ...
        'BackgroundColor','w','Margin',2,'Clipping','off');
    uistack(ht,'top');
end
ylim(axL, [0.5, 5+0.6]);
hBase = xline(axL, dp_base, ':', 'Color',[.35 .35 .35], 'LineWidth',1.2, 'DisplayName','baseline \DeltaPA');
hLo = plot(axL, NaN, NaN, 'o', 'MarkerFaceColor', cLo, 'MarkerEdgeColor','w','DisplayName','lower bound');
hHi = plot(axL, NaN, NaN, 'o', 'MarkerFaceColor', cHi, 'MarkerEdgeColor','w','DisplayName','upper bound');
set(axL,'YTick',y, 'YTickLabel', names_s);
xlabel(axL,'\DeltaPA [a.u.]');
%title(axL,' Local sensitivity (tornado)','Interpreter','tex');
legend(axL, [hLo hHi hBase], 'Location','northwest');
axR = nexttile; box(axR,'on'); grid(axR,'on');
bar(axR, absSRRCs, 'FaceColor',[0.23 0.49 0.77]);
xticks(axR, 1:5); xticklabels(axR, names_s);
ylabel(axR,'|SRRC|'); ylim(axR,[0 max(absSRRCs)*1.15]);
%title(axR, sprintf(' Global sensitivity (|SRRC|, LHS %d)', N));
sgn_tag = repmat({'(-)'}, size(sgn));  sgn_tag(sgn > 0) = {'(+)'};
for i = 1:numel(absSRRCs)
    text(axR, i, absSRRCs(i)*1.02, sprintf('%.2f %s', absSRRCs(i), sgn_tag{i}), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom');
end
%% =======================================================================
% =======================================================================
col_new.P1    = [239 116 95]/255;
col_new.P2    = [74 152 237]/255;
col_new.DELTA = [249 167 16]/255;
figure('Color','w','Name','Figure 3: PA and DeltaPA dual-axis');
tl = tiledlayout(3,1,'Padding','compact','TileSpacing','compact');
h_legend = gobjects(1,3);
ax1 = nexttile; 
yyaxis(ax1, 'left');
hold(ax1, 'on');
h_legend(1) = plot(ax1, Dt_list*1e6, PA1_Dt, '-o', 'Color', col_new.P1, 'MarkerFaceColor', col_new.P1, 'MarkerSize', 4, 'DisplayName', 'PA_1');
h_legend(2) = plot(ax1, Dt_list*1e6, PA2_Dt, '-o', 'Color', col_new.P2, 'MarkerFaceColor', col_new.P2, 'MarkerSize', 4, 'DisplayName', 'PA_2');
ylabel(ax1, 'PA [a.u.]');
ax1.YColor = [0 0 0]; 
grid(ax1, 'on');
smart_ylim(ax1, [PA1_Dt; PA2_Dt], 0.1);
yyaxis(ax1, 'right');
h_legend(3) = plot(ax1, Dt_list*1e6, DP_Dt, '-^', 'Color', col_new.DELTA, 'MarkerFaceColor', col_new.DELTA, 'MarkerSize', 4, 'LineWidth', 1.5, 'DisplayName', '\DeltaPA');
ylabel(ax1, '\DeltaPA [a.u.]');
ax1.YColor = col_new.DELTA;
smart_ylim(ax1, DP_Dt, 0.2);
hold(ax1, 'off');
xlabel(ax1, '\Delta t [\mus]');
box(ax1, 'on');
ax2 = nexttile; 
yyaxis(ax2, 'left');
hold(ax2, 'on');
plot(ax2, w0_list_um, PA1_w0, '-o', 'Color', col_new.P1, 'MarkerFaceColor', col_new.P1, 'MarkerSize', 4, 'HandleVisibility','off');
plot(ax2, w0_list_um, PA2_w0, '-o', 'Color', col_new.P2, 'MarkerFaceColor', col_new.P2, 'MarkerSize', 4, 'HandleVisibility','off');
ylabel(ax2, 'PA [a.u.]');
ax2.YColor = [0 0 0];
grid(ax2, 'on');
smart_ylim(ax2, [PA1_w0; PA2_w0], 0.1);
yyaxis(ax2, 'right');
plot(ax2, w0_list_um, DP_w0, '-^', 'Color', col_new.DELTA, 'MarkerFaceColor', col_new.DELTA, 'MarkerSize', 4, 'LineWidth', 1.5, 'HandleVisibility','off');
ylabel(ax2, '\DeltaPA [a.u.]');
ax2.YColor = col_new.DELTA;
smart_ylim(ax2, DP_w0, 0.2);
hold(ax2, 'off');
xlabel(ax2, 'w_0 [\mum]');
box(ax2, 'on');
ax3 = nexttile; 
yyaxis(ax3, 'left');
hold(ax3, 'on');
plot(ax3, scanTau, PA1_tL, '-o', 'Color', col_new.P1, 'MarkerFaceColor', col_new.P1, 'MarkerSize', 4, 'HandleVisibility','off');
plot(ax3, scanTau, PA2_tL, '-o', 'Color', col_new.P2, 'MarkerFaceColor', col_new.P2, 'MarkerSize', 4, 'HandleVisibility','off');
ylabel(ax3, 'PA [a.u.]');
ax3.YColor = [0 0 0]; 
grid(ax3, 'on');
smart_ylim(ax3, [PA1_tL; PA2_tL], 0.1); 
yyaxis(ax3, 'right');
plot(ax3, scanTau, DP_tL, '-^', 'Color', col_new.DELTA, 'MarkerFaceColor', col_new.DELTA, 'MarkerSize', 4, 'LineWidth', 1.5, 'HandleVisibility','off');
ylabel(ax3, '\DeltaPA [a.u.]');
ax3.YColor = col_new.DELTA; 
smart_ylim(ax3, DP_tL, 0.2); 
hold(ax3, 'off');
xlabel(ax3, '\tau_L [ns]');
box(ax3, 'on');
lg = legend(h_legend);
lg.Layout.Tile = 'south'; 
lg.Orientation = 'horizontal';


scan_f0_MHz = [15 20 25 30 35]; 
nF0 = numel(scan_f0_MHz);

scan_BW_MHz = [5 10 15 20]; 
nBW = numel(scan_BW_MHz);

F_list_sim = mJcm2_to_Jm2(scanF); 
nF_sim = numel(F_list_sim);

DP_virt_all = zeros(nF0, nBW, nF_sim);
Ratio_virt_all = zeros(nF0, nBW, nF_sim);

fprintf('Computing virtual-transducer effects (total combinations: %d)...\n', nF0*nBW);

for iF0 = 1:nF0
    for iBW = 1:nBW
        
        curr_f0 = scan_f0_MHz(iF0);
        curr_bw = scan_BW_MHz(iBW);
        
        h_virtual = apply_bandlimit(ref_wave, dt, curr_f0, curr_bw);
        
        % h_virtual = h_virtual / max(abs(h_virtual)); 
        
        for iF = 1:nF_sim
            Fi = F_list_sim(iF);
            
            [pinc1, pinc2] = build_pincs_prevOnly(t, dt, Fi, Fi, ...
                t1, Delta_t0, tauL0, mu0, eta_th, rho, Cp, ...
                Gamma0, dGamma_dT, tau_th0, gaussI);
            
            p1_v = conv(pinc1, h_virtual, 'full');
            p2_v = conv(pinc2, h_virtual, 'full');
            
            pa1_val = pp(p1_v);
            pa2_val = pp(p2_v);
            
            DP_virt_all(iF0, iBW, iF) = pa2_val - pa1_val;
            
            if pa1_val > 1e-6
                Ratio_virt_all(iF0, iBW, iF) = (pa2_val - pa1_val) / pa1_val;
            else
                Ratio_virt_all(iF0, iBW, iF) = 0;
            end
        end
    end
    fprintf('  > Completed f0 = %d MHz group.\n', curr_f0);
end

%% =======================================================================
% =======================================================================
set(groot,'defaultTextInterpreter','tex');
set(groot,'defaultAxesFontSize',10);
set(groot,'defaultLineLineWidth',1.5);


dt = 1e-9;
T_end = 150e-6; 
t = 0:dt:T_end;

base_F      = 150;           % [J/m^2] (15 mJ/cm^2)
base_mu     = 24000;         % [1/m]   (400 1/cm)
base_d      = 5e-6;          % [m]     (5 um)
base_Dt     = 40e-6;         % [s]     (40 us)
base_Tau    = 10e-9;         % [s]     (10 ns)

rho = 1000; Cp = 3600; 
eta_th = 1; Gamma0 = 0.13; 
dGamma_dT = 0.01;
alpha_th = 1.3e-7;

gaussI = @(F,tv,t0,tau,dtv) (F * exp(-0.5*((tv-t0)/(tau/(2*sqrt(2*log(2))))).^2)) / ...
                            (sum(exp(-0.5*((tv-t0)/(tau/(2*sqrt(2*log(2))))).^2))*dtv);
tau_th_func = @(d) d.^2 / alpha_th;
pp = @(x) max(x) - min(x);


ref_impulse = zeros(size(t));
start_idx = round(5e-6/dt); % 5us start
ref_impulse(start_idx) = 1; 

list_f0_MHz = [15, 20, 25, 30, 35];
fbw_exp     = 0.68; 
Filters_F0  = cell(1, numel(list_f0_MHz));
for k = 1:numel(list_f0_MHz)
    f_c = list_f0_MHz(k);
    bw  = f_c * fbw_exp; 
    Filters_F0{k} = apply_bandlimit_simple(ref_impulse, dt, f_c, bw);
end

fix_f0_MHz  = 25;
list_BW_MHz = [5, 10, 15, 20]; 
Filters_BW  = cell(1, numel(list_BW_MHz));
for k = 1:numel(list_BW_MHz)
    bw = list_BW_MHz(k);
    Filters_BW{k} = apply_bandlimit_simple(ref_impulse, dt, fix_f0_MHz, bw);
end


Tasks(1).name = 'mu_a';
Tasks(1).label = 'Absorption Coefficient \mu_a [1/cm]';
Tasks(1).x_scan = linspace(10, 1000, 20); % 1/cm
Tasks(1).x_unit_scale = 100; % cm^-1 -> m^-1

Tasks(2).name = 'd';
Tasks(2).label = 'Spot Radius d [\mum]';
Tasks(2).x_scan = linspace(2, 20, 20);    % um
Tasks(2).x_unit_scale = 1e-6; % um -> m

Tasks(3).name = 'Delta_t';
Tasks(3).label = 'Pulse Interval \Delta t [\mus]';
Tasks(3).x_scan = linspace(10, 100, 20);  % us
Tasks(3).x_unit_scale = 1e-6; % us -> s

Tasks(4).name = 'tau_L';
Tasks(4).label = 'Pulse Width \tau_L [ns]';
Tasks(4).x_scan = linspace(2, 40, 20);    % ns
Tasks(4).x_unit_scale = 1e-9; % ns -> s


colors = turbo(max(numel(list_f0_MHz), numel(list_BW_MHz)) + 1);

for iTask = 1:4
    T = Tasks(iTask);
    fprintf('Processing panel %d: %s ...\n', iTask, T.name);
    
    nPoints = length(T.x_scan);
    
    Res_F0_DPA   = zeros(nPoints, numel(list_f0_MHz));
    Res_F0_Ratio = zeros(nPoints, numel(list_f0_MHz));
    Res_BW_DPA   = zeros(nPoints, numel(list_BW_MHz));
    Res_BW_Ratio = zeros(nPoints, numel(list_BW_MHz));
    
    for iX = 1:nPoints
        val = T.x_scan(iX) * T.x_unit_scale;
        
        curr_mu  = base_mu;
        curr_d   = base_d;
        curr_Dt  = base_Dt;
        curr_Tau = base_Tau;
        
        switch T.name
            case 'mu_a',    curr_mu  = val;
            case 'd',       curr_d   = val;
            case 'Delta_t', curr_Dt  = val;
            case 'tau_L',   curr_Tau = val;
        end
        
        curr_tau_th = tau_th_func(curr_d);
        
        t1 = 10e-6; 
        I1 = gaussI(base_F, t, t1, curr_Tau, dt);
        I2 = gaussI(base_F, t, t1+curr_Dt, curr_Tau, dt);
        q1 = curr_mu * eta_th * I1;
        q2 = curr_mu * eta_th * I2;
        
        heat_impulse = exp(-t/curr_tau_th);
        T_rise = (1/(rho*Cp)) * conv(q1, heat_impulse, 'full') * dt;
        T_rise = T_rise(1:length(t));
        
        pinc1 = Gamma0 * q1;
        
        pinc2 = (Gamma0 + dGamma_dT * T_rise) .* q2; 
        
        for k = 1:numel(list_f0_MHz)
            h = Filters_F0{k};
            p1 = conv(pinc1, h, 'full'); 
            p2 = conv(pinc2, h, 'full');
            pa1 = pp(p1); pa2 = pp(p2);
            dpa = pa2 - pa1;
            
            Res_F0_DPA(iX, k) = dpa;
            if pa1 > 1e-9, Res_F0_Ratio(iX, k) = dpa / pa1; end
        end
        
        for k = 1:numel(list_BW_MHz)
            h = Filters_BW{k};
            p1 = conv(pinc1, h, 'full'); 
            p2 = conv(pinc2, h, 'full');
            pa1 = pp(p1); pa2 = pp(p2);
            dpa = pa2 - pa1;
            
            Res_BW_DPA(iX, k) = dpa;
            if pa1 > 1e-9, Res_BW_Ratio(iX, k) = dpa / pa1; end
        end
    end
    

    fig = figure('Color','w','Position',[100 100 1100 700], 'Name', ['Analysis: ' T.name]);
    tl = tiledlayout(2,2, 'Padding','compact', 'TileSpacing','compact');
    
    title(tl, ['Impact of Acoustic Parameters on ' T.name ' Sensitivity'], 'FontSize', 14, 'FontWeight','bold');
    
    % Plot 1: f0 Effect on ΔPA
    ax1 = nexttile; hold(ax1,'on'); grid(ax1,'on'); box(ax1,'on');
    for k = 1:numel(list_f0_MHz)
        plot(ax1, T.x_scan, Res_F0_DPA(:,k), 'o-', 'LineWidth',1.5, ...
            'Color', colors(k,:), 'MarkerFaceColor', colors(k,:), ...
            'MarkerSize', 4, ...
            'DisplayName', sprintf('f_0=%d MHz', list_f0_MHz(k)));
    end
    xlabel(ax1, T.label); ylabel(ax1, '\DeltaPA [a.u.]');
    title(ax1, 'Effect of f_0 on \DeltaPA (FBW=68%)');
    legend(ax1,'Location','best');
    
    % Plot 2: f0 Effect on Ratio
    ax2 = nexttile; hold(ax2,'on'); grid(ax2,'on'); box(ax2,'on');
    for k = 1:numel(list_f0_MHz)
        plot(ax2, T.x_scan, Res_F0_Ratio(:,k), 'o-', 'LineWidth',1.5, ...
            'Color', colors(k,:), 'MarkerFaceColor', colors(k,:), ...
            'MarkerSize', 4, 'HandleVisibility','off');
    end
    xlabel(ax2, T.label); ylabel(ax2, 'Ratio \DeltaPA/PA_1');
    title(ax2, 'Effect of f_0 on Ratio (FBW=68%)');
    
    % Plot 3: BW Effect on ΔPA
    ax3 = nexttile; hold(ax3,'on'); grid(ax3,'on'); box(ax3,'on');
    colors_bw = cool(numel(list_BW_MHz));
    for k = 1:numel(list_BW_MHz)
        plot(ax3, T.x_scan, Res_BW_DPA(:,k), 's-', 'LineWidth',1.5, ...
            'Color', colors_bw(k,:), 'MarkerFaceColor', colors_bw(k,:), ...
            'MarkerSize', 5, ...
            'DisplayName', sprintf('BW=%d MHz', list_BW_MHz(k)));
    end
    xlabel(ax3, T.label); ylabel(ax3, '\DeltaPA [a.u.]');
    title(ax3, sprintf('Effect of BW on \\DeltaPA (f_0=%d MHz)', fix_f0_MHz));
    legend(ax3,'Location','best');
    
    % Plot 4: BW Effect on Ratio
    ax4 = nexttile; hold(ax4,'on'); grid(ax4,'on'); box(ax4,'on');
    for k = 1:numel(list_BW_MHz)
        plot(ax4, T.x_scan, Res_BW_Ratio(:,k), 's-', 'LineWidth',1.5, ...
            'Color', colors_bw(k,:), 'MarkerFaceColor', colors_bw(k,:), ...
            'MarkerSize', 5, 'HandleVisibility','off');
    end
    xlabel(ax4, T.label); ylabel(ax4, 'Ratio \DeltaPA/PA_1');
    title(ax4, sprintf('Effect of BW on Ratio (f_0=%d MHz)', fix_f0_MHz));
    
end

fprintf('All plotting tasks completed.\n');


%% ================== Figure 4: hardware influence（retained version） ==================


plot_BW_MHz = [5 10 15];

idx_BW_plot = zeros(size(plot_BW_MHz));
for kk = 1:numel(plot_BW_MHz)
    [~, idx_BW_plot(kk)] = min(abs(scan_BW_MHz - plot_BW_MHz(kk)));
end

[~, idx_BW_fix] = min(abs(scan_BW_MHz - 10));
[~, idx_f0_fix] = min(abs(scan_f0_MHz - 25));


colors_f0 = [
    0.25 0.10 0.35   % f0 = 15 MHz
    0.20 0.70 0.90   % f0 = 20 MHz
    0.60 0.90 0.20   % f0 = 25 MHz
    0.95 0.50 0.10   % f0 = 30 MHz
    0.55 0.00 0.00   % f0 = 35 MHz
];

colors_bw = [
    0.10 0.90 0.95   % BW = 5 MHz
    0.35 0.60 0.95   % BW = 10 MHz
    1.00 0.00 1.00   % BW = 15 MHz
];

fig_fig4 = figure('Color','w', ...
    'Name','Figure 4: Hardware Influence', ...
    'Position',[100 80 900 980]);

tiledlayout(3, 2, ...
    'Padding','compact', ...
    'TileSpacing','compact');

%% ================= (a) spectral overlap =================
ax_a = nexttile([1 2]);
hold(ax_a, 'on');
box(ax_a, 'on');
grid(ax_a, 'off');

% experimental transducer response
[f_h, A_h] = one_sided_amp_spectrum(ref_wave, dt);
A_h_norm = A_h / max(A_h);
f_MHz = f_h / 1e6;

mask_f = f_MHz <= 60;
f_plot = f_MHz(mask_f);

A_h_plot = A_h_norm(mask_f);
A_h_plot = smoothdata(A_h_plot, 'gaussian', 9);
A_h_plot = A_h_plot / max(A_h_plot);

fc = 12.0;
p_shape = 2.0;

A_src_plot = exp(-(f_plot / fc).^p_shape);
A_src_plot = A_src_plot / max(A_src_plot);

h_src = patch(ax_a, ...
    [f_plot, fliplr(f_plot)], ...
    [A_src_plot, zeros(size(A_src_plot))], ...
    [0.55 0.55 0.55], ...
    'FaceAlpha', 0.75, ...
    'EdgeColor', 'none', ...
    'DisplayName', 'Analytical Source Spectrum');

h_tr = patch(ax_a, ...
    [f_plot, fliplr(f_plot)], ...
    [A_h_plot, zeros(size(A_h_plot))], ...
    [0.55 0.77 0.92], ...
    'FaceAlpha', 0.90, ...
    'EdgeColor', 'none', ...
    'DisplayName', 'Experimental Transducer Response');

xlabel(ax_a, 'Frequency [MHz]');
ylabel(ax_a, 'Normalized Amplitude [a.u.]');

xlim(ax_a, [0 60]);
ylim(ax_a, [0 1.02]);

legend(ax_a, [h_src h_tr], ...
    'Location', 'northeast', ...
    'Box', 'on');

text(ax_a, 4.0, 0.88, 'Source Energy', ...
    'Color', [0.20 0.20 0.20], ...
    'FontWeight', 'bold');

text(ax_a, 33.5, 0.40, 'Transducer Window', ...
    'Color', [0.00 0.45 0.80], ...
    'FontWeight', 'bold');

text(ax_a, -0.07, 1.03, '(a)', ...
    'Units', 'normalized', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom');

hold(ax_a, 'off');

%% ================= (b) ΔPA vs F for different f0 =================
ax_b = nexttile;
hold(ax_b, 'on');
box(ax_b, 'on');
grid(ax_b, 'off');

for i = 1:numel(scan_f0_MHz)
    plot(ax_b, scanF, squeeze(DP_virt_all(i, idx_BW_fix, :)), ...
        'o-', ...
        'LineWidth', 1.5, ...
        'Color', colors_f0(i,:), ...
        'MarkerFaceColor', colors_f0(i,:), ...
        'DisplayName', sprintf('f_0 = %d MHz', scan_f0_MHz(i)));
end

xlabel(ax_b, 'F [mJ/cm^2]');
ylabel(ax_b, '\DeltaPA [a.u.]');

legend(ax_b, 'Location', 'northwest', 'Box', 'on');

text(ax_b, -0.16, 1.03, '(b)', ...
    'Units', 'normalized', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom');

hold(ax_b, 'off');

%% ================= (c) nonlinear gain vs F for different f0 =================
ax_c = nexttile;
hold(ax_c, 'on');
box(ax_c, 'on');
grid(ax_c, 'off');

markers_f0 = {'o','s','d','^','v'};

for i = 1:numel(scan_f0_MHz)
    mk = markers_f0{mod(i-1, numel(markers_f0)) + 1};

    plot(ax_c, scanF(2:end), squeeze(Ratio_virt_all(i, idx_BW_fix, 2:end)), ...
        [mk '-'], ...
        'LineWidth', 1.5, ...
        'Color', colors_f0(i,:), ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', 'none', ...
        'DisplayName', sprintf('f_0 = %d MHz', scan_f0_MHz(i)));
end

xlabel(ax_c, 'F [mJ/cm^2]');
ylabel(ax_c, 'Nonlinear Gain \eta (\DeltaPA/PA_1)');
ylim(ax_c, [0 0.2]);

legend(ax_c, 'Location', 'northwest', 'Box', 'on');

text(ax_c, -0.16, 1.03, '(c)', ...
    'Units', 'normalized', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom');

hold(ax_c, 'off');

%% ================= (d) ΔPA vs F for different BW =================
ax_d = nexttile;
hold(ax_d, 'on');
box(ax_d, 'on');
grid(ax_d, 'off');

for jj = 1:numel(plot_BW_MHz)
    j = idx_BW_plot(jj);

    plot(ax_d, scanF, squeeze(DP_virt_all(idx_f0_fix, j, :)), ...
        's-', ...
        'LineWidth', 1.5, ...
        'Color', colors_bw(jj,:), ...
        'MarkerFaceColor', colors_bw(jj,:), ...
        'DisplayName', sprintf('BW = %d MHz', plot_BW_MHz(jj)));
end

xlabel(ax_d, 'F [mJ/cm^2]');
ylabel(ax_d, '\DeltaPA [a.u.]');

legend(ax_d, 'Location', 'northwest', 'Box', 'on');

text(ax_d, -0.16, 1.03, '(d)', ...
    'Units', 'normalized', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom');

hold(ax_d, 'off');

%% ================= (e) nonlinear gain vs F for different BW =================
ax_e = nexttile;
hold(ax_e, 'on');
box(ax_e, 'on');
grid(ax_e, 'off');

markers_bw = {'s','d','^'};

for jj = 1:numel(plot_BW_MHz)
    j = idx_BW_plot(jj);
    mk = markers_bw{jj};

    plot(ax_e, scanF(2:end), squeeze(Ratio_virt_all(idx_f0_fix, j, 2:end)), ...
        [mk '-'], ...
        'LineWidth', 1.5, ...
        'Color', colors_bw(jj,:), ...
        'MarkerFaceColor', 'none', ...
        'MarkerSize', 8, ...
        'DisplayName', sprintf('BW = %d MHz', plot_BW_MHz(jj)));
end

xlabel(ax_e, 'F [mJ/cm^2]');
ylabel(ax_e, 'Nonlinear Gain \eta (\DeltaPA/PA_1)');
ylim(ax_e, [0 0.2]);

legend(ax_e, 'Location', 'northwest', 'Box', 'on');

text(ax_e, -0.16, 1.03, '(e)', ...
    'Units', 'normalized', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom');

hold(ax_e, 'off');

%% ================= global formatting =================
all_ax = findall(fig_fig4, 'Type', 'axes');

set(all_ax, ...
    'FontName', 'Arial', ...
    'FontSize', 11, ...
    'LineWidth', 0.8, ...
    'Box', 'on', ...
    'XGrid', 'off', ...
    'YGrid', 'off', ...
    'TickDir', 'out');

set(findall(fig_fig4, '-property', 'FontName'), 'FontName', 'Arial');
set(findall(fig_fig4, '-property', 'FontSize'), 'FontSize', 11);
%% ================== Figure 4: hardware influence（retained polished version） ==================

scan_BW_MHz = [5 10 15];

idx_BW_fix = find(scan_BW_MHz == 10, 1);
idx_f0_fix = find(scan_f0_MHz == 25, 1);

if isempty(idx_BW_fix)
    error('scan_BW_MHz does not contain 10 MHz. Please check scan_BW_MHz.');
end

if isempty(idx_f0_fix)
    error('scan_f0_MHz does not contain 25 MHz. Please check scan_f0_MHz.');
end

colors_f0 = [
    0.25 0.10 0.35
    0.20 0.70 0.90
    0.60 0.90 0.20
    0.95 0.50 0.10
    0.55 0.00 0.00
];

colors_bw = [
    0.10 0.90 0.95
    0.35 0.60 0.95
    1.00 0.00 1.00
];

fig_fig4 = figure('Color','w', ...
    'Name','Figure 4: Hardware Influence', ...
    'Position', [100 80 900 980]);

tl = tiledlayout(3, 2, ...
    'Padding', 'compact', ...
    'TileSpacing', 'compact');

%% ================= (a) spectral overlap =================
ax_a = nexttile([1 2]);
hold(ax_a, 'on');
box(ax_a, 'on');
grid(ax_a, 'off');

% experimental transducer response
[f_h, A_h] = one_sided_amp_spectrum(ref_wave, dt);
A_h_norm = A_h / max(A_h);
f_MHz = f_h / 1e6;

fc = 12.0;
p_shape = 2.0;

A_src_smooth = exp(-(f_MHz / fc).^p_shape);
A_src_smooth = A_src_smooth / max(A_src_smooth);

mask_f = f_MHz <= 60;
f_plot = f_MHz(mask_f);
A_src_plot = A_src_smooth(mask_f);
A_h_plot = A_h_norm(mask_f);

h_src = patch(ax_a, ...
    [f_plot, fliplr(f_plot)], ...
    [A_src_plot, zeros(size(A_src_plot))], ...
    [0.55 0.55 0.55], ...
    'FaceAlpha', 0.75, ...
    'EdgeColor', 'none', ...
    'DisplayName', 'analytical source spectrum');

h_tr = patch(ax_a, ...
    [f_plot, fliplr(f_plot)], ...
    [A_h_plot, zeros(size(A_h_plot))], ...
    [0.55 0.77 0.92], ...
    'FaceAlpha', 0.90, ...
    'EdgeColor', 'none', ...
    'DisplayName', 'experimental transducer response');

xlabel(ax_a, 'frequency [MHz]');
ylabel(ax_a, 'normalized amplitude [a.u.]');

xlim(ax_a, [0 60]);
ylim(ax_a, [0 1.02]);

legend(ax_a, [h_src h_tr], ...
    'Location', 'northeast', ...
    'Box', 'on');

text(ax_a, 4.2, 0.88, 'source energy', ...
    'Color', [0.20 0.20 0.20], ...
    'FontWeight', 'bold');

text(ax_a, 32.8, 0.50, 'transducer window', ...
    'Color', [0.00 0.45 0.80], ...
    'FontWeight', 'bold');

text(ax_a, -0.06, 1.03, '(a)', ...
    'Units', 'normalized', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom');

hold(ax_a, 'off');

%% ================= (b) ΔPA vs F for different f0 =================
ax_b = nexttile;
hold(ax_b, 'on');
box(ax_b, 'on');
grid(ax_b, 'off');

for i = 1:numel(scan_f0_MHz)
    plot(ax_b, scanF, squeeze(DP_virt_all(i, idx_BW_fix, :)), ...
        'o-', ...
        'LineWidth', 1.5, ...
        'Color', colors_f0(i,:), ...
        'MarkerFaceColor', colors_f0(i,:), ...
        'DisplayName', sprintf('f_0 = %d MHz', scan_f0_MHz(i)));
end

xlabel(ax_b, 'F [mJ/cm^2]');
ylabel(ax_b, '\DeltaPA [a.u.]');

legend(ax_b, 'Location', 'northwest', 'Box', 'on');

text(ax_b, -0.16, 1.03, '(b)', ...
    'Units', 'normalized', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom');

hold(ax_b, 'off');

%% ================= (c) nonlinear gain vs F for different f0 =================
ax_c = nexttile;
hold(ax_c, 'on');
box(ax_c, 'on');
grid(ax_c, 'off');

markers_f0 = {'o', 's', 'd', '^', 'v'};

for i = 1:numel(scan_f0_MHz)
    mk = markers_f0{mod(i-1, numel(markers_f0)) + 1};

    plot(ax_c, scanF(2:end), squeeze(Ratio_virt_all(i, idx_BW_fix, 2:end)), ...
        [mk '-'], ...
        'LineWidth', 1.5, ...
        'Color', colors_f0(i,:), ...
        'MarkerSize', 8, ...
        'MarkerFaceColor', 'none', ...
        'DisplayName', sprintf('f_0 = %d MHz', scan_f0_MHz(i)));
end

xlabel(ax_c, 'F [mJ/cm^2]');
ylabel(ax_c, 'nonlinear gain \eta (\DeltaPA/PA_1)');
ylim(ax_c, [0 0.2]);

legend(ax_c, 'Location', 'northwest', 'Box', 'on');

text(ax_c, -0.16, 1.03, '(c)', ...
    'Units', 'normalized', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom');

hold(ax_c, 'off');

%% ================= (d) ΔPA vs F for different BW =================
ax_d = nexttile;
hold(ax_d, 'on');
box(ax_d, 'on');
grid(ax_d, 'off');

for j = 1:numel(scan_BW_MHz)
    plot(ax_d, scanF, squeeze(DP_virt_all(idx_f0_fix, j, :)), ...
        's-', ...
        'LineWidth', 1.5, ...
        'Color', colors_bw(j,:), ...
        'MarkerFaceColor', colors_bw(j,:), ...
        'DisplayName', sprintf('BW = %d MHz', scan_BW_MHz(j)));
end

xlabel(ax_d, 'F [mJ/cm^2]');
ylabel(ax_d, '\DeltaPA [a.u.]');

legend(ax_d, 'Location', 'northwest', 'Box', 'on');

text(ax_d, -0.16, 1.03, '(d)', ...
    'Units', 'normalized', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom');

hold(ax_d, 'off');

%% ================= (e) nonlinear gain vs F for different BW =================
ax_e = nexttile;
hold(ax_e, 'on');
box(ax_e, 'on');
grid(ax_e, 'off');

markers_bw = {'s', 'd', '^'};

for j = 1:numel(scan_BW_MHz)
    mk = markers_bw{mod(j-1, numel(markers_bw)) + 1};

    plot(ax_e, scanF(2:end), squeeze(Ratio_virt_all(idx_f0_fix, j, 2:end)), ...
        [mk '-'], ...
        'LineWidth', 1.5, ...
        'Color', colors_bw(j,:), ...
        'MarkerFaceColor', 'none', ...
        'MarkerSize', 8, ...
        'DisplayName', sprintf('BW = %d MHz', scan_BW_MHz(j)));
end

xlabel(ax_e, 'F [mJ/cm^2]');
ylabel(ax_e, 'nonlinear gain \eta (\DeltaPA/PA_1)');
ylim(ax_e, [0 0.2]);

legend(ax_e, 'Location', 'northwest', 'Box', 'on');

text(ax_e, -0.16, 1.03, '(e)', ...
    'Units', 'normalized', ...
    'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'bottom');

hold(ax_e, 'off');

%% ================= global formatting =================
all_ax = findall(fig_fig4, 'Type', 'axes');

set(all_ax, ...
    'FontName', 'Arial', ...
    'FontSize', 11, ...
    'LineWidth', 0.8, ...
    'Box', 'on', ...
    'XGrid', 'off', ...
    'YGrid', 'off', ...
    'TickDir', 'out');

set(findall(fig_fig4, '-property', 'FontName'), 'FontName', 'Arial');
set(findall(fig_fig4, '-property', 'FontSize'), 'FontSize', 11);


%% =======================================================================
% =======================================================================
fig_composite = figure('Color', 'w', 'Name', 'Composite figure (a)-(g)', 'Position', [100 100 1100 850]);

tl_main = tiledlayout(fig_composite, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');


tl_left = tiledlayout(tl_main, 2, 1, 'Padding', 'none', 'TileSpacing', 'compact');
tl_left.Layout.Tile = 1;

ax_a = nexttile(tl_left); hold(ax_a,'on');
hP1_a = plot(ax_a, scanF, PA1_F, '-', 'Color', col.P1, 'DisplayName', 'PA_1');
hP2_a = plot(ax_a, scanF, PA2_F, '-', 'Color', col.P2, 'DisplayName', 'PA_2');
grid(ax_a,'on'); xlim(ax_a, [min(scanF) max(scanF)]);
xlabel(ax_a, 'F [mJ/cm^2]'); ylabel(ax_a, 'PA [a.u.]');
smart_ylim(ax_a, [PA1_F; PA2_F], 0.05);
[h1_a, h5_a, hMPE_a] = annotate_linear_regions_F(ax_a, scanF, mu0, Gamma0, dGamma_dT, rho, Cp, Delta_t0, tau_th0, MPE_mJcm2, cLIN, cWNL, cNL, cSAFE, true);
legend(ax_a, [hP1_a, hP2_a, h1_a, h5_a, hMPE_a], 'Location', 'southeast');
title(ax_a, '(a)', 'Units', 'normalized', 'Position', [-0.08, 1.02], 'HorizontalAlignment', 'left', 'FontWeight', 'bold');
box(ax_a, 'on');

ax_b = nexttile(tl_left); hold(ax_b,'on'); xlim(ax_b, [10 1000]);
hP1_b = plot(ax_b, mu_list_cm, PA1_mu, '-', 'Color', col.P1, 'DisplayName', 'PA_1');
hP2_b = plot(ax_b, mu_list_cm, PA2_mu, '-', 'Color', col.P2, 'DisplayName', 'PA_2');
grid(ax_b,'on'); xlabel(ax_b, '\mu_a [1/cm]'); ylabel(ax_b, 'PA [a.u.]');
smart_ylim(ax_b, [PA1_mu; PA2_mu], 0.05);
[h1_b, h5_b] = annotate_threezone_mu(ax_b, [10 1000], F0, Gamma0, dGamma_dT, rho, Cp, Delta_t0, tau_th0, cLIN, cWNL, cSAFE, true);
legend(ax_b, [hP1_b, hP2_b, h1_b, h5_b], 'Location', 'southeast');
title(ax_b, '(b)', 'Units', 'normalized', 'Position', [-0.08, 1.02], 'HorizontalAlignment', 'left', 'FontWeight', 'bold');
box(ax_b, 'on');

for i = 1:size(absorber_data, 1)
    name = absorber_data{i, 1}; mu_val = absorber_data{i, 2}; cat = absorber_data{i, 3};
    marker_color = cat_colors_bright.(cat); marker_style = cat_markers.(cat);
    current_marker_size = 8;
    if strcmp(marker_style, 'o'), current_marker_size = 7;
    elseif strcmp(marker_style, 'p'), current_marker_size = 11;
    elseif strcmp(marker_style, 'd'), current_marker_size = 7; end
    
    if mu_val >= 10 && mu_val <= 1000
        pa2_val = interp1(mu_list_cm, PA2_mu, mu_val);
        plot(ax_b, mu_val, pa2_val, marker_style, 'MarkerFaceColor', marker_color, ...
            'MarkerEdgeColor', 'none', 'MarkerSize', current_marker_size, 'HandleVisibility', 'off');
        text(ax_b, mu_val + 15, pa2_val, name, 'HorizontalAlignment', 'left', ...
            'VerticalAlignment', 'middle', 'Color', 'k', 'Interpreter', 'tex', 'Margin', 2);
    end
end


tl_right = tiledlayout(tl_main, 4, 1, 'Padding', 'none', 'TileSpacing', 'compact');
tl_right.Layout.Tile = 2;

tl_cd = tiledlayout(tl_right, 1, 2, 'Padding', 'none', 'TileSpacing', 'compact');
tl_cd.Layout.Tile = 1;

maskF = DP_F > 0; pfit_F = polyfit(log10(F_list(maskF)), log10(DP_F(maskF)), 1);
maskMu = DP_mu > 0; pfit_mu = polyfit(log10(mu_list(maskMu)), log10(DP_mu(maskMu)), 1);

% >>> (c) Log-log: Delta PA vs F
ax_c = nexttile(tl_cd);
loglog(ax_c, F_list(maskF), DP_F(maskF), 'o', 'Color', col.DELTA); hold(ax_c,'on');
loglog(ax_c, F_list(maskF), 10.^polyval(pfit_F, log10(F_list(maskF))), '-', 'Color', col.GRAY);
grid(ax_c,'on'); xlabel(ax_c, 'F [J/m^2]'); ylabel(ax_c, '\DeltaPA [a.u.]');
title(ax_c, '(c)', 'Units', 'normalized', 'Position', [-0.15, 1.02], 'HorizontalAlignment', 'left', 'FontWeight', 'bold');
box(ax_c, 'on');

% >>> (d) Log-log: Delta PA vs mu_a
ax_d = nexttile(tl_cd);
loglog(ax_d, mu_list(maskMu), DP_mu(maskMu), 'o', 'Color', col.DELTA); hold(ax_d,'on');
loglog(ax_d, mu_list(maskMu), 10.^polyval(pfit_mu, log10(mu_list(maskMu))), '-', 'Color', col.GRAY);
grid(ax_d,'on'); xlabel(ax_d, '\mu_a [1/m]'); ylabel(ax_d, '\DeltaPA [a.u.]');
title(ax_d, '(d)', 'Units', 'normalized', 'Position', [-0.15, 1.02], 'HorizontalAlignment', 'left', 'FontWeight', 'bold');
box(ax_d, 'on');

c_p1 = col_new.P1; c_p2 = col_new.P2; c_dp = col_new.DELTA;
h_leg_bot = gobjects(1,3);

ax_e = nexttile(tl_right);
yyaxis(ax_e, 'left'); hold(ax_e, 'on');
h_leg_bot(1) = plot(ax_e, Dt_list*1e6, PA1_Dt, '-o', 'Color', c_p1, 'MarkerFaceColor', c_p1, 'MarkerSize', 4, 'DisplayName', 'PA_1');
h_leg_bot(2) = plot(ax_e, Dt_list*1e6, PA2_Dt, '-o', 'Color', c_p2, 'MarkerFaceColor', c_p2, 'MarkerSize', 4, 'DisplayName', 'PA_2');
ylabel(ax_e, 'PA [a.u.]'); ax_e.YColor = [0 0 0]; grid(ax_e, 'on');
smart_ylim(ax_e, [PA1_Dt; PA2_Dt], 0.1);
yyaxis(ax_e, 'right');
h_leg_bot(3) = plot(ax_e, Dt_list*1e6, DP_Dt, '-^', 'Color', c_dp, 'MarkerFaceColor', c_dp, 'MarkerSize', 4, 'LineWidth', 1.5, 'DisplayName', '\DeltaPA');
ylabel(ax_e, '\DeltaPA [a.u.]'); ax_e.YColor = c_dp; smart_ylim(ax_e, DP_Dt, 0.2);
xlabel(ax_e, '\Delta t [\mus]'); box(ax_e, 'on');
title(ax_e, '(e)', 'Units', 'normalized', 'Position', [-0.08, 1.02], 'HorizontalAlignment', 'left', 'FontWeight', 'bold');

ax_f = nexttile(tl_right);
yyaxis(ax_f, 'left'); hold(ax_f, 'on');
plot(ax_f, w0_list_um, PA1_w0, '-o', 'Color', c_p1, 'MarkerFaceColor', c_p1, 'MarkerSize', 4);
plot(ax_f, w0_list_um, PA2_w0, '-o', 'Color', c_p2, 'MarkerFaceColor', c_p2, 'MarkerSize', 4);
ylabel(ax_f, 'PA [a.u.]'); ax_f.YColor = [0 0 0]; grid(ax_f, 'on');
smart_ylim(ax_f, [PA1_w0; PA2_w0], 0.1);
yyaxis(ax_f, 'right');
plot(ax_f, w0_list_um, DP_w0, '-^', 'Color', c_dp, 'MarkerFaceColor', c_dp, 'MarkerSize', 4, 'LineWidth', 1.5);
ylabel(ax_f, '\DeltaPA [a.u.]'); ax_f.YColor = c_dp; smart_ylim(ax_f, DP_w0, 0.2);
xlabel(ax_f, 'd [\mum]'); box(ax_f, 'on');  
title(ax_f, '(f)', 'Units', 'normalized', 'Position', [-0.08, 1.02], 'HorizontalAlignment', 'left', 'FontWeight', 'bold');

ax_g = nexttile(tl_right);
yyaxis(ax_g, 'left'); hold(ax_g, 'on');
plot(ax_g, scanTau, PA1_tL, '-o', 'Color', c_p1, 'MarkerFaceColor', c_p1, 'MarkerSize', 4);
plot(ax_g, scanTau, PA2_tL, '-o', 'Color', c_p2, 'MarkerFaceColor', c_p2, 'MarkerSize', 4);
ylabel(ax_g, 'PA [a.u.]'); ax_g.YColor = [0 0 0]; grid(ax_g, 'on');
smart_ylim(ax_g, [PA1_tL; PA2_tL], 0.1);
yyaxis(ax_g, 'right');
plot(ax_g, scanTau, DP_tL, '-^', 'Color', c_dp, 'MarkerFaceColor', c_dp, 'MarkerSize', 4, 'LineWidth', 1.5);
ylabel(ax_g, '\DeltaPA [a.u.]'); ax_g.YColor = c_dp; smart_ylim(ax_g, DP_tL, 0.2);
xlabel(ax_g, '\tau_L [ns]'); box(ax_g, 'on');
title(ax_g, '(g)', 'Units', 'normalized', 'Position', [-0.08, 1.02], 'HorizontalAlignment', 'left', 'FontWeight', 'bold');

lg = legend(ax_g, h_leg_bot, 'Orientation', 'horizontal');
lg.Layout.Tile = 'south';


%% =======================================================================
% =======================================================================

fig_final = figure('Color', 'w', 'Name', 'Comprehensive Multi-parameter Analysis', 'Position', [100 100 1000 900]);
tl = tiledlayout(fig_final, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

set(findall(fig_final, '-property', 'FontName'), 'FontName', 'Arial');
set(findall(fig_final, '-property', 'FontSize'), 'FontSize', 11);


ax1 = nexttile(tl, 1); hold(ax1, 'on'); box(ax1, 'on'); grid(ax1, 'on');

imagesc(ax1, w0_vec, Dt_vec * 1e6, DP_map);
axis(ax1, 'xy');
colormap(ax1, parula);
cb1 = colorbar(ax1); 
cb1.Label.String = '\DeltaPA [a.u.]';
xlabel(ax1, 'd [\mum]'); 
ylabel(ax1, '\Delta t [\mus]');

beta_levels = [0.03, 0.1, 0.3, 1];
clr_beta = [0 114 178; 230 159 0; 0 158 115; 148 0 211]/255;
h_beta = gobjects(numel(beta_levels), 1);

for k = 1:numel(beta_levels)
    Dt_beta = beta_levels(k) * (w0_vec * 1e-6).^2 / alpha_th;
    h_beta(k) = plot(ax1, w0_vec, Dt_beta * 1e6, '-', 'Color', clr_beta(k,:), 'LineWidth', 2, ...
        'DisplayName', sprintf('\\beta = %g', beta_levels(k)));
end

legend(ax1, h_beta, 'Location', 'northwest', 'Box', 'on');
xlim(ax1, [min(w0_vec) max(w0_vec)]);
ylim(ax1, [min(Dt_vec * 1e6) max(Dt_vec * 1e6)]);


ax2 = nexttile(tl, 2); hold(ax2, 'on'); box(ax2, 'on'); grid(ax2, 'on');

imagesc(ax2, mu_vec, F_vec, DP_Fmu);
axis(ax2, 'xy');
colormap(ax2, parula);
cb2 = colorbar(ax2); 
cb2.Label.String = '\DeltaPA [a.u.]';
xlabel(ax2, '\mu_a [1/cm]'); 
ylabel(ax2, 'F per pulse [mJ/cm^2]');

h1 = plot(ax2, mu_vec, F1_line, '--', 'Color', [0.8 0 0.8], 'LineWidth', 2.5, 'DisplayName', '1%');
h5 = plot(ax2, mu_vec, F5_line, '--', 'Color', [0 0.8 0], 'LineWidth', 2.5, 'DisplayName', '5%');

yline(ax2, MPE_mJcm2, '--', 'Color', [0.85 0.2 0.2], 'LineWidth', 1.5, 'HandleVisibility', 'off');
yl = ylim(ax2);
if yl(2) > MPE_mJcm2
    patch(ax2, [mu_vec(1) mu_vec(end) mu_vec(end) mu_vec(1)], ...
        [MPE_mJcm2 MPE_mJcm2 yl(2) yl(2)], [0.5 0.5 0.5], ...
        'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    text(ax2, mu_vec(end)*0.95, MPE_mJcm2*1.05, 'MPE', 'Color', [0.85 0.2 0.2], ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'right', 'BackgroundColor', 'w');
end

legend(ax2, [h1, h5], 'Location', 'northwest', 'Box', 'on');
xlim(ax2, [min(mu_vec) max(mu_vec)]);
ylim(ax2, [min(F_vec) max(F_vec)]);


ax3 = nexttile(tl, 3); hold(ax3, 'on'); box(ax3, 'on');

num_vars = numel(names_s);
y_pos = 1:num_vars;
cLo = [0.30 0.45 0.85];
cHi = [0.93 0.49 0.18];

for r = 1:num_vars
    xL = min(dp_lo_s(r), dp_hi_s(r));
    xR = max(dp_lo_s(r), dp_hi_s(r));
    
    plot(ax3, [xL xR], [r r], '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 18, 'HandleVisibility', 'off');
    
    plot(ax3, dp_lo_s(r), r, 'o', 'MarkerFaceColor', cLo, 'MarkerEdgeColor', 'w', 'MarkerSize', 6, 'HandleVisibility', 'off');
    plot(ax3, dp_hi_s(r), r, 'o', 'MarkerFaceColor', cHi, 'MarkerEdgeColor', 'w', 'MarkerSize', 6, 'HandleVisibility', 'off');
    
    half_diff = round(abs(dp_hi_s(r) - dp_lo_s(r)) / 2);
    xmid = (xL + xR) / 2;
    text(ax3, xmid, r + 0.35, sprintf('\\pm %d a.u.', half_diff), 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontWeight', 'bold', 'BackgroundColor', 'w', 'Margin', 1);
end

hBase = xline(ax3, dp_base, 'k:', 'LineWidth', 1.2, 'DisplayName', 'baseline');
hLo = plot(ax3, NaN, NaN, 'o', 'MarkerFaceColor', cLo, 'MarkerEdgeColor', 'w', 'DisplayName', 'lower bound');
hHi = plot(ax3, NaN, NaN, 'o', 'MarkerFaceColor', cHi, 'MarkerEdgeColor', 'w', 'DisplayName', 'upper bound');

yticks(ax3, y_pos);
yticklabels(ax3, names_s);
xlabel(ax3, '\DeltaPA [a.u.]');
ylim(ax3, [0.5, num_vars + 0.8]);
legend(ax3, [hLo, hHi, hBase], 'Location', 'northeast', 'Box', 'on');


ax4 = nexttile(tl, 4); hold(ax4, 'on'); box(ax4, 'on'); grid(ax4, 'on');

bar(ax4, y_pos, absSRRCs, 0.65, 'FaceColor', [0.23 0.49 0.77], 'EdgeColor', 'k', 'LineWidth', 1);

for r = 1:num_vars
    text(ax4, r, absSRRCs(r) + 0.02, sprintf('%.2f %s', absSRRCs(r), sgn_tag{r}), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontWeight', 'bold');
end

xticks(ax4, y_pos);
xticklabels(ax4, names_s);
ylabel(ax4, '|SRRC|');
ylim(ax4, [0, max(absSRRCs) * 1.2]);

fprintf('Composite 2x2 figure generated.\n');

%% =======================================================================
%  Figure 1(c) final version - updated
%  beta-zone map + eta design map + sensitivity ranking
%
% =======================================================================



fontName = 'Arial';

fsTick   = 10;
fsLabel  = 10;
fsTitle  = 10.5;
fsBeta   = 10.5;
fsEta    = 10.0;
lwAxis   = 0.85;
lwLine   = 1.45;



beta_zone_colors = [
    221  89  76
    238 141  75
    244 204 100
    166 202 112
     93 145 194
] / 255;

eta_anchor = [
     63 110 183
    106 190 170
    255 239 170
    242 156  76
    210  70  64
] / 255;

rank_anchor = [
    210  70  64
    235 116  55
    245 181  67
    157 195  89
     84 132 190
] / 255;

cmap_eta = make_cmap(eta_anchor, 256);

col_th1  = [0.96 0.82 0.96];
col_th5  = [1.00 1.00 1.00];
col_mpe  = [0.52 0.80 0.96];
fill_mpe = [0.82 0.89 0.96];



Dt_vec = round(linspace(10e-6, 200e-6, 151) / dt) * dt;
d_vec  = linspace(2, 20, 181);

[D_um, DT] = meshgrid(d_vec, Dt_vec);

Beta_map = DT ./ ((D_um * 1e-6).^2 / alpha_th);

Beta_zone = nan(size(Beta_map));
Beta_zone(Beta_map < 0.03)                    = 1;
Beta_zone(Beta_map >= 0.03 & Beta_map < 0.1) = 2;
Beta_zone(Beta_map >= 0.1  & Beta_map < 0.3) = 3;
Beta_zone(Beta_map >= 0.3  & Beta_map < 1.0) = 4;
Beta_zone(Beta_map >= 1.0)                   = 5;



F_vec  = linspace(0.5, 25, 80);
mu_vec = linspace(100, 1000, 80);
[MU, FF] = meshgrid(mu_vec, F_vec);

PA1_Fmu = zeros(size(FF));
PA2_Fmu = zeros(size(FF));
DP_Fmu  = zeros(size(FF));

Delta_t_use = Delta_t0;
d_use_um    = w0_um0;
tauL_use    = tauL_ns0 * 1e-9;
tau_th_use  = (d_use_um * 1e-6)^2 / alpha_th;

for ii = 1:numel(FF)

    F_i  = mJcm2_to_Jm2(FF(ii));
    mu_i = cm1_to_m1(MU(ii));

    [pinc1i, pinc2i] = build_pincs_prevOnly( ...
        t, dt, F_i, F_i, t1, Delta_t_use, tauL_use, mu_i, ...
        eta_th, rho, Cp, Gamma0, dGamma_dT, tau_th_use, gaussI);

    p1i = conv(pinc1i, ref_wave, 'full');
    p2i = conv(pinc2i, ref_wave, 'full');

    PA1_Fmu(ii) = pp(p1i);
    PA2_Fmu(ii) = pp(p2i);
    DP_Fmu(ii)  = PA2_Fmu(ii) - PA1_Fmu(ii);
end

Eta_Fmu = DP_Fmu ./ max(PA1_Fmu, eps);

Z2 = Eta_Fmu;
cb2_label = '\eta = \DeltaPA / PA_1';

F1_line = zeros(size(mu_vec));
F5_line = zeros(size(mu_vec));

for i = 1:numel(mu_vec)
    [F1_line(i), F5_line(i)] = F_thresholds( ...
        cm1_to_m1(mu_vec(i)), Gamma0, dGamma_dT, rho, Cp, ...
        Delta_t_use, tau_th_use);
end


%% ===================== 4. E5: sensitivity ranking =====================
lofac = 0.8;
hifac = 1.2;

rng(1);
N = 600;

U = lhsdesign(N, 5, 'criterion', 'maximin', 'iterations', 10);

mu_s  = mu_a_cm0    * (lofac + (hifac - lofac) * U(:,1));
F_s   = F_mJcm2_0   * (lofac + (hifac - lofac) * U(:,2));
d_s   = w0_um0      * (lofac + (hifac - lofac) * U(:,3));
Dt_s  = Delta_t0    * (lofac + (hifac - lofac) * U(:,4));
tau_s = tauL_ns0    * (lofac + (hifac - lofac) * U(:,5));

Y = zeros(N, 1);

for n = 1:N
    Y(n) = evalDP_equalF( ...
        F_s(n), mu_s(n), Dt_s(n), d_s(n), tau_s(n), ...
        t, dt, t1, ref_wave, eta_th, rho, Cp, Gamma0, dGamma_dT, alpha_th);
end

X = [mu_s, F_s, d_s, Dt_s, tau_s];

Xr = tiedrank(X);
Yr = tiedrank(Y);

Xz = zscore(Xr);
Yz = zscore(Yr);

betaSRRC = Xz \ Yz;

[absSRRCs, k_rank] = sort(abs(betaSRRC(:)), 'descend');
sgnSRRC = sign(betaSRRC(k_rank));

name_all = {'\mu_a', 'F', 'd', '\Delta t', '\tau_L'};
names_s  = name_all(k_rank);

rankValue   = absSRRCs(:);
rankName    = names_s(:);
rankSign    = sgnSRRC(:);
rank_colors = make_cmap(rank_anchor, numel(rankValue));



figC = figure('Color','w', ...
    'Name','Figure_1c_final_design_maps_sensitivity_updated', ...
    'Units','centimeters', ...
    'Position',[3 3 25.5 7.8]);

tl = tiledlayout(figC, 1, 3, ...
    'Padding','compact', ...
    'TileSpacing','compact');


%% -----------------------------------------------------------------------
%  Panel 1: beta-zone map
% ------------------------------------------------------------------------
ax1 = nexttile(tl, 1);

imagesc(ax1, Dt_vec * 1e6, d_vec, Beta_zone.');
axis(ax1, 'xy');
colormap(ax1, beta_zone_colors);
clim(ax1, [1 5]);

hold(ax1, 'on');

beta_levels = [0.03, 0.1, 0.3, 1];
label_frac = [0.78, 0.66, 0.56, 0.44];

for kk = 1:numel(beta_levels)

    bval = beta_levels(kk);
    d_beta_um = sqrt(Dt_vec * alpha_th / bval) * 1e6;
    valid = d_beta_um >= min(d_vec) & d_beta_um <= max(d_vec);

    plot(ax1, Dt_vec(valid) * 1e6, d_beta_um(valid), ...
        '--', 'Color', [1 1 1], 'LineWidth', 1.55);

    idx_valid = find(valid);
    if ~isempty(idx_valid)
        idx_lab = idx_valid(max(1, round(label_frac(kk) * numel(idx_valid))));
        text(ax1, Dt_vec(idx_lab) * 1e6, d_beta_um(idx_lab), ...
            sprintf('\\beta = %.2g', bval), ...
            'FontName', fontName, ...
            'FontSize', fsBeta, ...
            'FontWeight', 'bold', ...
            'Color', [0.10 0.10 0.10], ...
            'BackgroundColor', [1 1 1], ...
            'Margin', 1.8, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'Clipping', 'on');
    end
end

text(ax1, 0.04, 0.96, ...
    {'small \beta', 'strong thermal memory'}, ...
    'Units','normalized', ...
    'FontName',fontName, ...
    'FontSize',8.0, ...
    'FontWeight','bold', ...
    'Color',[0.20 0.10 0.08], ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top', ...
    'BackgroundColor',[1 1 1], ...
    'Margin',2);

text(ax1, 0.96, 0.06, ...
    {'large \beta', 'weak thermal memory'}, ...
    'Units','normalized', ...
    'FontName',fontName, ...
    'FontSize',8.0, ...
    'FontWeight','bold', ...
    'Color',[0.08 0.14 0.24], ...
    'HorizontalAlignment','right', ...
    'VerticalAlignment','bottom', ...
    'BackgroundColor',[1 1 1], ...
    'Margin',2);

xlabel(ax1, '\Delta t [\mus]', 'FontName', fontName, 'FontSize', fsLabel);
ylabel(ax1, 'd [\mum]',        'FontName', fontName, 'FontSize', fsLabel);
title(ax1, '\beta-zone design map', ...
    'FontName', fontName, 'FontSize', fsTitle, 'FontWeight', 'normal');

set(ax1, ...
    'FontName', fontName, ...
    'FontSize', fsTick, ...
    'LineWidth', lwAxis, ...
    'Box', 'off', ...
    'TickLength', [0 0], ...
    'Layer', 'top');

grid(ax1, 'off');


%% -----------------------------------------------------------------------
%  Panel 2: eta map in (F, mu_a) space
% ------------------------------------------------------------------------
ax2 = nexttile(tl, 2);

contourf(ax2, F_vec, mu_vec, Z2.', 100, 'LineColor','none');
axis(ax2, 'xy');
colormap(ax2, cmap_eta);

zValid = Z2(isfinite(Z2));
if ~isempty(zValid)
    clim(ax2, prctile(zValid(:), [2 98]));
end

cb2 = colorbar(ax2);
cb2.Label.String = cb2_label;
cb2.Label.FontName = fontName;
cb2.Label.FontSize = fsLabel;
cb2.FontName = fontName;
cb2.FontSize = fsTick;
cb2.Box = 'off';

hold(ax2, 'on');

yl = [min(mu_vec), max(mu_vec)];

if max(F_vec) > MPE_mJcm2
    patch(ax2, ...
        [MPE_mJcm2, max(F_vec), max(F_vec), MPE_mJcm2], ...
        [yl(1), yl(1), yl(2), yl(2)], ...
        fill_mpe, ...
        'FaceAlpha', 0.28, ...
        'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
end

hth1 = plot(ax2, F1_line, mu_vec, '-.', ...
    'Color', col_th1, ...
    'LineWidth', 1.45, ...
    'DisplayName', '1% threshold');

hth5 = plot(ax2, F5_line, mu_vec, '--', ...
    'Color', col_th5, ...
    'LineWidth', 1.75, ...
    'DisplayName', '5% threshold');

idx_eta1 = round(0.58 * numel(mu_vec));
idx_eta5 = round(0.40 * numel(mu_vec));

text(ax2, F1_line(idx_eta1), mu_vec(idx_eta1), ...
    '\eta = 1%', ...
    'FontName', fontName, ...
    'FontSize', fsEta, ...
    'FontWeight', 'bold', ...
    'Color', [0.18 0.10 0.22], ...
    'BackgroundColor', [1 1 1], ...
    'Margin', 1.5, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'Clipping', 'on');

text(ax2, F5_line(idx_eta5), mu_vec(idx_eta5), ...
    '\eta = 5%', ...
    'FontName', fontName, ...
    'FontSize', fsEta, ...
    'FontWeight', 'bold', ...
    'Color', [0.12 0.12 0.12], ...
    'BackgroundColor', [1 1 1], ...
    'Margin', 1.5, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'Clipping', 'on');

xline(ax2, MPE_mJcm2, '--', ...
    'Color', col_mpe, ...
    'LineWidth', 1.55, ...
    'HandleVisibility', 'off');

text(ax2, MPE_mJcm2 * 0.98, yl(2) * 0.96, ...
    'MPE boundary', ...
    'FontName', fontName, ...
    'FontSize', 8.2, ...
    'FontWeight', 'bold', ...
    'Color', [0.08 0.32 0.50], ...
    'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'top', ...
    'BackgroundColor', [1 1 1], ...
    'Margin', 1.5);

xlabel(ax2, 'F [mJ/cm^2]', 'FontName', fontName, 'FontSize', fsLabel);
ylabel(ax2, '\mu_a [cm^{-1}]', 'FontName', fontName, 'FontSize', fsLabel);
title(ax2, '\eta design map', ...
    'FontName', fontName, 'FontSize', fsTitle, 'FontWeight', 'normal');

xlim(ax2, [0.5 25]);

legend(ax2, [hth1, hth5], ...
    'Location', 'northwest', ...
    'FontName', fontName, ...
    'FontSize', 7.5, ...
    'Box', 'off');

set(ax2, ...
    'FontName', fontName, ...
    'FontSize', fsTick, ...
    'LineWidth', lwAxis, ...
    'Box', 'off', ...
    'TickLength', [0 0], ...
    'Layer', 'top');

grid(ax2, 'off');


%% -----------------------------------------------------------------------
%  Panel 3: sensitivity ranking
% ------------------------------------------------------------------------
ax3 = nexttile(tl, 3);

y = 1:numel(rankValue);

b = barh(ax3, y, rankValue, 0.62);
b.FaceColor = 'flat';

for ii = 1:numel(rankValue)
    b.CData(ii,:) = rank_colors(ii,:);
end

set(ax3, ...
    'YTick', y, ...
    'YTickLabel', rankName, ...
    'YDir', 'reverse', ...
    'FontName', fontName, ...
    'FontSize', fsTick, ...
    'LineWidth', lwAxis, ...
    'Box', 'off', ...
    'TickLength', [0 0], ...
    'Layer', 'top');

xlabel(ax3, '|SRRC|', 'FontName', fontName, 'FontSize', fsLabel);
title(ax3, 'Sensitivity ranking', ...
    'FontName', fontName, 'FontSize', fsTitle, 'FontWeight', 'normal');

xlim(ax3, [0, max(rankValue) * 1.20]);

grid(ax3, 'off');

for ii = 1:numel(rankValue)

    if rankSign(ii) >= 0
        signText = '+';
    else
        signText = '-';
    end

    text(ax3, rankValue(ii) + max(rankValue) * 0.025, ii, ...
        sprintf('%.2f (%s)', rankValue(ii), signText), ...
        'FontName', fontName, ...
        'FontSize', 8.0, ...
        'VerticalAlignment', 'middle', ...
        'HorizontalAlignment', 'left', ...
        'Color', [0.20 0.20 0.20]);
end

text(ax3, 0.98, 1.04, ...
    'Higher \rightarrow more sensitive', ...
    'Units', 'normalized', ...
    'FontName', fontName, ...
    'FontSize', 8.2, ...
    'HorizontalAlignment', 'right', ...
    'Color', [0.25 0.25 0.25]);



set(findall(figC, '-property', 'FontName'), 'FontName', fontName);

exportgraphics(figC, fullfile(output_dir, 'Fig1c_final_design_maps_sensitivity_updated.png'), 'Resolution', 600);

fprintf('Updated design-map/sensitivity figure exported to outputs/.\n');


%% ========================================================================
%  Save source data and remaining open figures
% ========================================================================
fprintf('Saving source data and open figures to: %s\n', output_dir);

try
    Fig2 = struct();
    Fig2.scanF_mJcm2 = scanF;
    Fig2.PA1_F = PA1_F;
    Fig2.PA2_F = PA2_F;
    Fig2.DP_F  = DP_F;
    Fig2.scanMu_cm = scanMu;
    Fig2.PA1_mu = PA1_mu;
    Fig2.PA2_mu = PA2_mu;
    Fig2.DP_mu  = DP_mu;
    Fig2.scanDt_us = Dt_list * 1e6;
    Fig2.PA1_Dt = PA1_Dt;
    Fig2.PA2_Dt = PA2_Dt;
    Fig2.DP_Dt  = DP_Dt;
    Fig2.scanD_um = w0_list_um;
    Fig2.PA1_d = PA1_w0;
    Fig2.PA2_d = PA2_w0;
    Fig2.DP_d  = DP_w0;
    Fig2.scanTau_ns = scanTau;
    Fig2.PA1_tauL = PA1_tL;
    Fig2.PA2_tauL = PA2_tL;
    Fig2.DP_tauL  = DP_tL;
    Fig2.baseline = struct('F_mJcm2', F_mJcm2_0, 'mu_a_cm', mu_a_cm0, ...
        'd_um', w0_um0, 'Delta_t_us', Delta_t0*1e6, 'tauL_ns', tauL_ns0, ...
        'Gamma0', Gamma0, 'dGamma_dT', dGamma_dT, 'rho', rho, 'Cp', Cp, ...
        'alpha_th', alpha_th, 'MPE_mJcm2', MPE_mJcm2);
    save(fullfile(output_dir, 'Fig2_source_data.mat'), 'Fig2');
catch ME
    warning('Could not save Fig2 source data: %s', ME.message);
end

try
    Fig3 = struct();
    Fig3.w0_vec_um = w0_vec;
    Fig3.Dt_vec_us = Dt_vec * 1e6;
    Fig3.DP_map = DP_map;
    Fig3.beta_levels = beta_levels;
    Fig3.F_vec_mJcm2 = F_vec;
    Fig3.mu_vec_cm = mu_vec;
    Fig3.DP_Fmu = DP_Fmu;
    Fig3.F1_threshold_mJcm2 = F1_line;
    Fig3.F5_threshold_mJcm2 = F5_line;
    Fig3.tornado = struct('baseline', dp_base, 'lower', dp_lo, 'upper', dp_hi, ...
        'names', {xtags_all});
    Fig3.SRRC = struct('beta', betaSRRC, 'abs_sorted', absSRRCs, ...
        'names_sorted', {names_s}, 'sign_sorted', sgn);
    save(fullfile(output_dir, 'Fig3_source_data.mat'), 'Fig3');
catch ME
    warning('Could not save Fig3 source data: %s', ME.message);
end

try
    Fig4 = struct();
    Fig4.scanF_mJcm2 = scanF;
    Fig4.scan_f0_MHz = scan_f0_MHz;
    Fig4.scan_BW_MHz = scan_BW_MHz;
    Fig4.DP_virt_all = DP_virt_all;
    Fig4.Ratio_virt_all = Ratio_virt_all;
    Fig4.ref_wave = ref_wave;
    Fig4.dt = dt;
    save(fullfile(output_dir, 'Fig4_source_data.mat'), 'Fig4');
catch ME
    warning('Could not save Fig4 source data: %s', ME.message);
end

try
    figs = findall(0, 'Type', 'figure');
    [~, order] = sort([figs.Number]);
    figs = figs(order);
    for kk = 1:numel(figs)
        fig = figs(kk);
        figName = get(fig, 'Name');
        if isempty(figName)
            figName = sprintf('Figure_%02d', fig.Number);
        end
        safeName = regexprep(figName, '[^A-Za-z0-9_\-]+', '_');
        safeName = regexprep(safeName, '_+', '_');
        safeName = regexprep(safeName, '^_|_$', '');
        if isempty(safeName)
            safeName = sprintf('Figure_%02d', fig.Number);
        end
        baseFile = sprintf('%02d_%s', kk, safeName);
        exportgraphics(fig, fullfile(output_dir, [baseFile '.png']), 'Resolution', 600);
        savefig(fig, fullfile(output_dir, [baseFile '.fig']));
    end
catch ME
    warning('Could not export all open figures: %s', ME.message);
end

fprintf('Done. Main outputs are available in ./outputs/.\n');

%% ========================================================================
%  Local functions
% ========================================================================

function pretty_axis(ax, fontName, fontSize_axis, fontSize_label)
    box(ax, 'on');
    grid(ax, 'on');

    ax.FontName = fontName;
    ax.FontSize = fontSize_axis;
    ax.LineWidth = 0.95;
    ax.TickDir = 'out';
    ax.TickLength = [0.012 0.012];
    ax.Layer = 'top';

    ax.GridAlpha = 0.12;
    ax.MinorGridAlpha = 0.06;
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';

    ax.XColor = [0.18 0.18 0.18];

    if numel(ax.YAxis) == 1
        ax.YColor = [0.18 0.18 0.18];
    end

    ax.XLabel.FontSize = fontSize_label;
    ax.YLabel.FontSize = fontSize_label;
end


function lim = padded_ylim(y, frac)
    y = y(:);
    y = y(isfinite(y));

    if isempty(y)
        lim = [0 1];
        return;
    end

    ymin = min(y);
    ymax = max(y);

    if abs(ymax - ymin) < eps
        pad = max(abs(ymax), 1) * frac;
    else
        pad = (ymax - ymin) * frac;
    end

    lim = [ymin - pad, ymax + pad];

    if lim(1) == lim(2)
        lim = lim + [-1 1];
    end
end


function x0 = find_crossing_level(x, y, level)
    x = x(:);
    y = y(:);

    valid = isfinite(x) & isfinite(y);
    x = x(valid);
    y = y(valid);

    if numel(x) < 2
        x0 = NaN;
        return;
    end

    [x, id] = sort(x);
    y = y(id);

    if all(y < level)
        x0 = NaN;
        return;
    end

    if y(1) >= level
        x0 = x(1);
        return;
    end

    idx = find(y(1:end-1) < level & y(2:end) >= level, 1, 'first');

    if isempty(idx)
        x0 = NaN;
        return;
    end

    if abs(y(idx+1) - y(idx)) < eps
        x0 = x(idx);
    else
        x0 = interp1(y(idx:idx+1), x(idx:idx+1), level, 'linear');
    end
end


function draw_x_zone(ax, x1, x2, colorRGB, alphaVal)
    if isnan(x1) || isnan(x2)
        return;
    end

    xl = xlim(ax);
    yl = ylim(ax);

    x1 = max(x1, xl(1));
    x2 = min(x2, xl(2));

    if x2 <= x1
        return;
    end

    p = patch(ax, ...
        [x1 x2 x2 x1], ...
        [yl(1) yl(1) yl(2) yl(2)], ...
        colorRGB, ...
        'FaceAlpha', alphaVal, ...
        'EdgeColor', 'none', ...
        'HandleVisibility', 'off');

    try
        uistack(p, 'bottom');
    catch
    end
end


function draw_vline_label(ax, x0, labelStr, colorRGB, lineStyle)
    if isnan(x0)
        return;
    end

    xl = xlim(ax);
    yl = ylim(ax);

    if x0 < xl(1) || x0 > xl(2)
        return;
    end

    plot(ax, [x0 x0], yl, lineStyle, ...
        'Color', colorRGB, ...
        'LineWidth', 1.45, ...
        'HandleVisibility', 'off');

    text(ax, x0, yl(2) - 0.045 * range(yl), labelStr, ...
        'Rotation', 90, ...
        'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', ...
        'Color', colorRGB, ...
        'FontWeight', 'bold', ...
        'FontSize', 10, ...
        'BackgroundColor', 'w', ...
        'Margin', 1.0, ...
        'Clipping', 'on');
end


function add_panel_label(ax, str, fs)
    text(ax, -0.080, 1.045, str, ...
        'Units', 'normalized', ...
        'FontSize', fs, ...
        'FontWeight', 'bold', ...
        'Color', [0 0 0], ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'bottom', ...
        'Clipping', 'off');
end


function plot_dual_y(ax, x, pa1, pa2, dpa, xlabelStr, ylabelLeft, ylabelRight, C, lw, ms)
    x   = x(:).';
    pa1 = pa1(:).';
    pa2 = pa2(:).';
    dpa = dpa(:).';

    yyaxis(ax, 'left');
    hold(ax, 'on');

    plot(ax, x, pa1, '-o', ...
        'Color', C.pa1, ...
        'MarkerFaceColor', C.pa1, ...
        'MarkerEdgeColor', 'w', ...
        'MarkerSize', ms, ...
        'LineWidth', lw);

    plot(ax, x, pa2, '-o', ...
        'Color', C.pa2, ...
        'MarkerFaceColor', C.pa2, ...
        'MarkerEdgeColor', 'w', ...
        'MarkerSize', ms, ...
        'LineWidth', lw);

    ylabel(ax, ylabelLeft);
    ylim(ax, padded_ylim([pa1 pa2], 0.11));
    ax.YColor = [0.18 0.18 0.18];

    yyaxis(ax, 'right');
    hold(ax, 'on');

    plot(ax, x, dpa, '-^', ...
        'Color', C.dpa, ...
        'MarkerFaceColor', C.dpa, ...
        'MarkerEdgeColor', 'w', ...
        'MarkerSize', ms, ...
        'LineWidth', lw);

    ylabel(ax, ylabelRight);
    ylim(ax, padded_ylim(dpa, 0.16));
    ax.YColor = C.dpa;

    xlabel(ax, xlabelStr);

    xlim(ax, [min(x), max(x)]);
    box(ax, 'on');
    grid(ax, 'on');
end

function h_filt = apply_bandlimit_simple(h, dt, f0_MHz, BW_MHz)
    h = h(:)';
    N = numel(h);
    H = fft(h);
    f = (0:N-1)/(N*dt); % Hz
    f0 = f0_MHz * 1e6;
    bw = BW_MHz * 1e6;
    mask = (abs(f - f0) <= bw/2) | (abs(f - (1/dt) + f0) <= bw/2);
    H_new = zeros(size(H));
    H_new(mask) = H(mask);
    h_filt = real(ifft(H_new));
end



function [p_inc1, p_inc2] = build_pincs_prevOnly( ...
        t, dt, F1, F2, t1, Dt, tauL, mu, eta, rho, Cp, G0, dGdT, tau_th, gaussI)
    I1 = gaussI(F1,t,t1,      tauL,dt);
    I2 = gaussI(F2,t,t1 + Dt, tauL,dt);
    q1 = mu*eta*I1;  q2 = mu*eta*I2;          % J/m^3/s
    T1 = 0;  T1_hist = zeros(size(t));
    for n = 1:numel(t)
        T1_hist(n) = T1;
        T1 = T1*exp(-dt/tau_th) + q1(n)*dt/(rho*Cp);
    end
    p_inc1 = G0               * q1 * dt;     % Pa
    p_inc2 = (G0 + dGdT*T1_hist) .* q2 * dt; % Pa
    %p_inc1 = G0 * q1;
    %p_inc2 = (G0 + dGdT*T1_hist) .* q2;
end
% function annotate_mpe_noOverlap(ax, x_thresh_mJcm2, x_max)
%     yl = ylim(ax);
%     p = patch(ax, [x_thresh_mJcm2 x_max x_max x_thresh_mJcm2], ...
%                  [yl(1) yl(1) yl(2) yl(2)], [1 0.8 0.8], ...
%                  'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
%     uistack(p,'bottom');
%     % xline(ax, x_thresh_mJcm2, '--', 'Color',[0.85 0.2 0.2], 'LineWidth',1.2, 'HandleVisibility','off');
%     % text(ax, x_thresh_mJcm2*1.01, yl(2)*0.95, 'MPE 20 mJ/cm^2', ...
%     %     'Color',[0.85 0.2 0.2], 'FontWeight','bold', ...
%     %     'HorizontalAlignment','left','VerticalAlignment','top', ...
%     %     'BackgroundColor','none');
%     hMPE = xline(ax, x_thresh_mJcm2, '--', 'Color',[0.85 0.2 0.2], ...
%                  'LineWidth',1.2, 'DisplayName','MPE ');
% end
function hMPE = annotate_mpe_noOverlap(ax, x_thresh_mJcm2, x_max)
    yl = ylim(ax);
    p = patch(ax, [x_thresh_mJcm2 x_max x_max x_thresh_mJcm2], ...
                 [yl(1) yl(1) yl(2) yl(2)], [1 0.8 0.8], ...
                 'FaceAlpha',0.18,'EdgeColor','none','HandleVisibility','off');
    uistack(p,'bottom');
    hMPE = xline(ax, x_thresh_mJcm2, '--', 'Color',[0.85 0.2 0.2], ...
                 'LineWidth',1.2, 'DisplayName','MPE 20 mJ/cm^2');
end
% function annotate_linear_regions_F(ax, scanF_mJcm2, mu_m1, G0, dGdT, rho, Cp, Dt, tau_th, MPE_mJcm2, cLIN, cWNL, cNL, cSAFE, drawBlue)
%     xL = [min(scanF_mJcm2), max(scanF_mJcm2)];
%     [F1, F5] = F_thresholds(mu_m1, G0, dGdT, rho, Cp, Dt, tau_th); % [mJ/cm^2]
%     ph = gobjects(0);
%     ph(end+1) = fill_band(ax, xL(1), min(F1,xL(2)), cLIN);          % <1%
%     ph(end+1) = fill_band(ax, max(F1,xL(1)), min(F5,xL(2)), cWNL);  % 1–5%
%     ph(end+1) = fill_band(ax, max(F5,xL(1)), xL(2), cNL);           % ≥5%
%     FsafeL = max(F5, xL(1)); FsafeR = min(MPE_mJcm2, xL(2));
%     if FsafeR > FsafeL
%         ph(end+1) = fill_band(ax, FsafeL, FsafeR, cSAFE, 0.18);
%         % yl = ylim(ax);
%         % text(ax, mean([FsafeL FsafeR]), yl(2)*0.92, 'safe NL', ...
%         %     'HorizontalAlignment','center','VerticalAlignment','top','Color',[0 0.45 0],'FontWeight','bold');
%     end
%     stretch_patches_to_ylim(ax, ph);
%     % yl = ylim(ax);
%     % if F1>xL(1) && F1<xL(2)
%     %     xline(ax, F1, ':', 'Color',[0.55 0.55 0.55], 'LineWidth',1.2, 'HandleVisibility','off');
%     %     text(ax, F1, yl(2)*0.90, '1%', 'Color',[0.35 0.35 0.35], ...
%     %          'HorizontalAlignment','center','VerticalAlignment','bottom');
%     % end
%     % if drawBlue && F5>xL(1) && F5<xL(2)
%     %     xline(ax, F5, '--', 'Color',[0 0.45 0.85], 'LineWidth',1.3, 'HandleVisibility','off');
%     %     text(ax, F5, yl(2)*0.84, '5%', 'Color',[0 0.45 0.85], ...
%     %          'HorizontalAlignment','center','VerticalAlignment','bottom');
%     % end
%     % annotate_mpe_noOverlap(ax, MPE_mJcm2, xL(2));
% 
%     h1 = []; h5 = [];
%     if F1>xL(1) && F1<xL(2), h1 = xline(ax, F1, '-.', 'Color',[0.00 0.45 0.74], 'LineWidth',1.5, 'DisplayName','1%'); end
%     if F5>xL(1) && F5<xL(2), h5 = xline(ax, F5, '--', 'Color',[0 0.45 0.85],   'LineWidth',1.5, 'DisplayName','5%');  end
% 
%     hMPE = annotate_mpe_noOverlap(ax, MPE_mJcm2, xL(2));
% end
function [h1,h5,hMPE] = annotate_linear_regions_F(ax, scanF_mJcm2, mu_m1, G0, dGdT, rho, Cp, Dt, tau_th, MPE_mJcm2, cLIN, cWNL, cNL, cSAFE, drawBlue)
    xL = [min(scanF_mJcm2), max(scanF_mJcm2)];
    [F1, F5] = F_thresholds(mu_m1, G0, dGdT, rho, Cp, Dt, tau_th);
    ph = gobjects(0);
    ph(end+1) = fill_band(ax, xL(1), min(F1,xL(2)), cLIN, 0.25);
    ph(end+1) = fill_band(ax, max(F1,xL(1)), min(F5,xL(2)), cWNL, 0.18);
    ph(end+1) = fill_band(ax, max(F5,xL(1)), xL(2),            cNL,  0.18);
    FsafeL = max(F5, xL(1)); FsafeR = min(MPE_mJcm2, xL(2));
    if FsafeR > FsafeL
        ph(end+1) = fill_band(ax, FsafeL, FsafeR, cSAFE, 0.18);
    end
    stretch_patches_to_ylim(ax, ph);
    h1 = gobjects(1); h5 = gobjects(1);
    if drawBlue && F1>xL(1) && F1<xL(2)
        h1 = xline(ax, F1, '-.', 'Color',[0.00 0.45 0.74], 'LineWidth',1.5, 'DisplayName','1%');
    else, h1 = gobjects(0); end
    if drawBlue && F5>xL(1) && F5<xL(2)
        h5 = xline(ax, F5, '--',  'Color',[0 0.45 0.85],   'LineWidth',1.5, 'DisplayName','5%');
    else, h5 = gobjects(0); end
    % MPE
    hMPE = annotate_mpe_noOverlap(ax, MPE_mJcm2, xL(2));
end
% function annotate_threezone_mu(ax, xL_cm, F_Jm2, G0, dGdT, rho, Cp, Dt, tau_th, cLIN, cWNL, cSAFE, drawBlue)
%     [mu1_cm, mu5_cm] = mu_thresholds(F_Jm2, G0, dGdT, rho, Cp, Dt, tau_th); % [1/cm]
%     ph = gobjects(0);
%     ph(end+1) = fill_band(ax, xL_cm(1), min(mu1_cm, xL_cm(2)), cLIN);      % <1%
%     ph(end+1) = fill_band(ax, max(mu1_cm, xL_cm(1)), min(mu5_cm, xL_cm(2)), cWNL); % 1–5%
%     ph(end+1) = fill_band(ax, max(mu5_cm, xL_cm(1)), xL_cm(2), cSAFE, 0.18);       % ≥5%
%     % yl = ylim(ax);
%     % if mu1_cm>xL_cm(1) && mu1_cm<xL_cm(2)
%     %     xline(ax, mu1_cm, ':', 'Color',[0.55 0.55 0.55], 'LineWidth',1.2, 'HandleVisibility','off');
%     %     text(ax, mu1_cm, yl(2)*0.90, '1%', 'Color',[0.35 0.35 0.35], ...
%     %          'HorizontalAlignment','center','VerticalAlignment','bottom');
%     % end
%     % if drawBlue && mu5_cm>xL_cm(1) && mu5_cm<xL_cm(2)
%     %     xline(ax, mu5_cm, '--', 'Color',[0 0.45 0.85], 'LineWidth',1.3, 'HandleVisibility','off');
%     %     text(ax, mu5_cm, yl(2)*0.84, '5%', 'Color',[0 0.45 0.85], ...
%     %          'HorizontalAlignment','center','VerticalAlignment','bottom');
%     % end
%      h1 = []; h5 = [];
%     if mu1_cm>xL_cm(1) && mu1_cm<xL_cm(2)
%         h1 = xline(ax, mu1_cm, '-.', 'Color',[0.00 0.45 0.74], 'LineWidth',1.5, 'DisplayName','1%');
%     end
%     if mu5_cm>xL_cm(1) && mu5_cm<xL_cm(2)
%         h5 = xline(ax, mu5_cm, '--',  'Color',[0 0.45 0.85],   'LineWidth',1.5, 'DisplayName','5%');
%     end
% end
function [h1,h5] = annotate_threezone_mu(ax, xL_cm, F_Jm2, G0, dGdT, rho, Cp, Dt, tau_th, cLIN, cWNL, cSAFE, drawBlue)
    [mu1_cm, mu5_cm] = mu_thresholds(F_Jm2, G0, dGdT, rho, Cp, Dt, tau_th);
    ph = gobjects(0);
    ph(end+1) = fill_band(ax, xL_cm(1), min(mu1_cm, xL_cm(2)), cLIN, 0.25);
    ph(end+1) = fill_band(ax, max(mu1_cm, xL_cm(1)), min(mu5_cm, xL_cm(2)), cWNL, 0.18);
    ph(end+1) = fill_band(ax, max(mu5_cm, xL_cm(1)), xL_cm(2), cSAFE, 0.18);
    stretch_patches_to_ylim(ax, ph);
    h1 = gobjects(0); h5 = gobjects(0);
    if drawBlue && mu1_cm>xL_cm(1) && mu1_cm<xL_cm(2)
        h1 = xline(ax, mu1_cm, '-.', 'Color',[0.00 0.45 0.74], 'LineWidth',1.5, 'DisplayName','1%');
    end
    if drawBlue && mu5_cm>xL_cm(1) && mu5_cm<xL_cm(2)
        h5 = xline(ax, mu5_cm, '--',  'Color',[0 0.45 0.85],   'LineWidth',1.5, 'DisplayName','5%');
    end
end
function [F1_mJcm2, F5_mJcm2] = F_thresholds(mu_m1, G0, dGdT, rho, Cp, Dt, tau_th)
    eta1 = 0.01; eta5 = 0.05;
    fac  = (G0*rho*Cp)/(dGdT*mu_m1) * exp(Dt/tau_th); % [J/m^2] per unit η
    F1_mJcm2 = (eta1*fac)/10;   % J/m^2 -> mJ/cm^2
    F5_mJcm2 = (eta5*fac)/10;
end
function [mu1_cm, mu5_cm] = mu_thresholds(F_Jm2, G0, dGdT, rho, Cp, Dt, tau_th)
    eta1 = 0.01; eta5 = 0.05;
    fac  = (G0*rho*Cp)/(dGdT*F_Jm2) * exp(Dt/tau_th); % [1/m] per unit η
    mu1_cm = (eta1*fac)/100;   % 1/m -> 1/cm
    mu5_cm = (eta5*fac)/100;
end
function p = fill_band(ax, xL, xR, rgb, alpha)
    if nargin<5, alpha=0.28; end
    if xR<=xL, p = gobjects(1); return; end
    yl = ylim(ax);
    p = patch(ax, [xL xR xR xL], [yl(1) yl(1) yl(2) yl(2)], rgb, ...
          'FaceAlpha',alpha,'EdgeColor','none','HandleVisibility','off');
    uistack(p,'bottom');
end
function stretch_patches_to_ylim(ax, ph)
    yl = ylim(ax);
    for k = 1:numel(ph)
        if isgraphics(ph(k),'patch')
            ph(k).YData = [yl(1) yl(1) yl(2) yl(2)];
        end
    end
end
function [f, A1] = one_sided_amp_spectrum(x, dt)
    N = numel(x); NFFT = 2^nextpow2(N);
    X = fft(x, NFFT); A = abs(X)/NFFT;
    A1 = 2*A(1:NFFT/2); f = (0:NFFT/2-1) / (NFFT*dt);
end
function [fpk,fL,fR,thr] = peak_and_bw_db(A, f, dB)
    if nargin<3, dB=6; end
    A1=A; A1(1)=0; [Amax, idx]=max(A1); fpk=f(idx); thr=Amax*10^(-dB/20);
    fL=NaN; fR=NaN;
    for i=idx-1:-1:2, if A1(i)<=thr, fL=interp1([A1(i) A1(i+1)],[f(i) f(i+1)],thr); break; end, end
    for i=idx+1:numel(A1)-1, if A1(i)<=thr, fR=interp1([A1(i-1) A1(i)],[f(i-1) f(i)],thr); break; end, end
end
function fancy_spectrum(ax, f, A, col, label, dB)
    if nargin<6, dB=6; end
    [fpk,fL,fR,thr] = peak_and_bw_db(A,f,dB);
    plot(ax, f/1e6, A, 'Color', col, 'LineWidth',1.5,'DisplayName',label); hold(ax,'on');
    if ~isnan(fL) && ~isnan(fR)
        idx = f>=fL & f<=fR;
        fill(ax, [f(idx)/1e6, fliplr(f(idx)/1e6)], [A(idx), zeros(1,sum(idx))], ...
             col, 'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
    end
    [~,ipk]=max(A); ypk=A(ipk);
    plot(ax, fpk/1e6, ypk, 'o','MarkerFaceColor',col,'MarkerEdgeColor','w','HandleVisibility','off');
    text(ax, fpk/1e6, ypk*1.02, sprintf(' %.2f MHz', fpk/1e6), ...
         'Color',col,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','bottom');
    if ~isnan(fL), xline(ax, fL/1e6, ':', 'Color', col, 'HandleVisibility','off'); end
    if ~isnan(fR), xline(ax, fR/1e6, ':', 'Color', col, 'HandleVisibility','off'); end
    if ~isnan(fL) && ~isnan(fR)
        bw=(fR-fL)/1e6; text(ax,(fL+fR)/2/1e6, thr*0.95, sprintf('BW_{-6dB}=%.2f MHz',bw), ...
            'Color',col,'HorizontalAlignment','center','VerticalAlignment','top');
    end
    box(ax,'on'); grid(ax,'on');
end
function smart_ylim(ax, y, margin)
    if nargin<3, margin=0.08; end
    ymin = min(y,[],'all'); ymax = max(y,[],'all');
    if ~isfinite(ymin) || ~isfinite(ymax) || ymax<=ymin, return; end
    pad = (ymax - ymin) * margin;
    ylim(ax, [ymin - pad, ymax + pad]);
end
function dp_kPa = evalDP_equalF(F_mJ, mu_cm, Dt, w0_um, tauL_ns, ...
                                t, dt, t1, ref_wave, eta_th, rho, Cp, G0, dGdT, alpha_th)
    FJ   = F_mJ * 10;            % mJ/cm^2 -> J/m^2
    mu_m = mu_cm * 100;          % 1/cm -> 1/m
    tauL = tauL_ns * 1e-9;
    tau_th = (w0_um*1e-6)^2 / alpha_th;
    [pinc1,pinc2] = build_pincs_prevOnly(t,dt,FJ,FJ,t1,Dt,tauL,mu_m,eta_th, ...
                                         rho,Cp,G0,dGdT,tau_th, ...
                                         @(F,tv,t0,tau,dtv) (F*exp(-0.5*((tv-t0)/(tau/(2*sqrt(2*log(2))))).^2)) / ...
                                             (sum(exp(-0.5*((tv-t0)/(tau/(2*sqrt(2*log(2))))).^2))*dtv));
    p1 = conv(pinc1, ref_wave, 'full');
    p2 = conv(pinc2, ref_wave, 'full');
    dp_kPa = ((max(p2)-min(p2)) - (max(p1)-min(p1))) ;
end
function [P1, P2, DP, t_long] = synth_equalF(F_mJ, mu_cm, Dt, w0_um, tauL_ns, ...
                                             t, dt, t1, ref_wave, eta_th, rho, Cp, G0, dGdT, alpha_th)
    FJ     = F_mJ * 10;          % mJ/cm^2 -> J/m^2
    mu_m   = mu_cm * 100;        % 1/cm    -> 1/m
    tauL   = tauL_ns * 1e-9;     % ns      -> s
    tau_th = (w0_um*1e-6)^2 / alpha_th;
    gaussI = @(F,tv,t0,tauFWHM,dtv) ...
        (F * exp(-0.5*((tv - t0)/(tauFWHM/(2*sqrt(2*log(2))))).^2) ) / ...
        (sum(exp(-0.5*((tv - t0)/(tauFWHM/(2*sqrt(2*log(2))))).^2))*dtv);
    [pinc1,pinc2] = build_pincs_prevOnly( ...
        t, dt, FJ, FJ, t1, Dt, tauL, mu_m, eta_th, rho, Cp, G0, dGdT, tau_th, gaussI);
    P1 = conv(pinc1, ref_wave, 'full');
    P2 = conv(pinc2, ref_wave, 'full');
    DP     = P2 - P1;
    t_long = (0:numel(P1)-1) * dt;
end
function forceTimes(fig)
    objs = findall(fig,'-property','FontName');
    for k = 1:numel(objs)
        try
            set(objs(k), 'FontName', 'Times New Roman');
        catch
        end
    end
end
function h_bp = apply_bandlimit(h, dt, f0_MHz, BW_MHz)
% apply_bandlimit:
%   EN: apply a simple rectangular band-pass filter to emulate a transducer
%       with center frequency f0 and bandwidth BW.
    h = h(:)';
    N = numel(h);
    H = fft(h);
    f = (0:N-1)/(N*dt);         % Hz
    f0 = f0_MHz * 1e6;
    BW = BW_MHz * 1e6;
    mask_pos = abs(f - f0) <= BW/2;
    mask_neg = abs((f - 1/dt) + f0) <= BW/2;
    mask = mask_pos | mask_neg;
    Hf = zeros(size(H));
    Hf(mask) = H(mask);
    h_bp = real(ifft(Hf));
end
function y = bandpass_fft(x, dt, f1, f2)
    x = x(:).';
    N = numel(x);
    X = fft(x);
    f = (0:N-1)/(N*dt);
    mask_pos = (f >= f1) & (f <= f2);
    mask_neg = (f >= (1/dt - f2)) & (f <= (1/dt - f1));
    mask = mask_pos | mask_neg;
    Xf = zeros(size(X));
    Xf(mask) = X(mask);
    y = real(ifft(Xf));
end
function ppv = pp_gate_centered(p, tlong, tcenter, gate_us, base_us)
    gate = gate_us * 1e-6;
    base = base_us * 1e-6;
    idxB = (tlong >= (tcenter - gate - base)) & (tlong <= (tcenter - gate));
    if any(idxB)
        p0 = mean(p(idxB));
    else
        p0 = mean(p(1:min(50,end)));
    end
    p2 = p - p0;
    idxG = (tlong >= (tcenter - gate)) & (tlong <= (tcenter + gate));
    if ~any(idxG)
        ppv = max(p2) - min(p2);
        return;
    end
    seg = p2(idxG);
    ppv = max(seg) - min(seg);
end

function cmap = make_cmap(anchor, n)
    if nargin < 2
        n = 256;
    end
    x  = linspace(0, 1, size(anchor, 1));
    xi = linspace(0, 1, n);
    cmap = interp1(x, anchor, xi, 'pchip');
    cmap = max(min(cmap, 1), 0);
end
