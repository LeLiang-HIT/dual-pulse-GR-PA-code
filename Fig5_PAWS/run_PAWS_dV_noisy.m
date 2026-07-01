%% ========================================================================
% PAWS nonlinear photoacoustic wavefront shaping demo
%
% This script implements a two-stage genetic-algorithm PAWS workflow:
%   1. linear pre-focusing by maximizing C1;
%   2. nonlinear optimization by maximizing noisy DeltaV feedback.

clear; clc; close all;
rng(8);                         % Fixed random seed for reproducibility

fprintf('PAWS dV noisy-feedback simulation started.\n');

%% ========================================================================
%  Section 1: Physical parameters and simulation setup
%  ========================================================================

% --- Acoustic medium parameters ---
c0 = 1500;              
rho = 1000;             
Cp = 3600;              
alpha_coeff = 0.75;     
alpha_power = 1.5;      

% --- Gruneisen nonlinearity parameters ---
mu_a_cm = 240;          
mu_a = mu_a_cm*100;     
Gamma0 = 0.12;          
dGamma_dT = 0.01;       

% --- Thermal relaxation and dual-pulse timing ---
alpha_th = 1.3e-7;              
d_c = 5e-6;                     
tau_th = d_c^2/alpha_th;        
t1_us = 10;                     
Delta_t_us = 40;                
t1 = t1_us*1e-6; 
Delta_t = Delta_t_us*1e-6;
memory = exp(-Delta_t/tau_th);  

% --- Fluence setting for waveform demonstration ---
F_plot_mJcm2 = 5;  
F_plot = F_plot_mJcm2 * 10;     
F_list_mJcm2 = F_plot_mJcm2; 
F_list = F_plot;

% --- k-Wave grid settings ---
Nx = 128; Ny = 128; 
dx = 50e-6; dy = dx; 
CFL = 0.1; 
PML = 20;
T_end = 70e-6;

kgrid = kWaveGrid(Nx, dx, Ny, dy);
dt = CFL * dx / c0; 
Nt = ceil(T_end/dt); 
kgrid.setTime(Nt, dt);

medium.sound_speed = c0;
medium.density     = rho;
medium.alpha_coeff = alpha_coeff;
medium.alpha_power = alpha_power;

% --- Absorber definition ---
[xc, yc] = deal(round(Nx/2), round(Ny/2));
abs_r_pix = 10;
source_mask = makeDisc(Nx, Ny, xc, yc, abs_r_pix);

% --- Sensor definition ---
sensor_row = PML + 6; 
aperture_half = 8;
sensor.mask = zeros(Nx, Ny);
sensor.mask(sensor_row, max(1, yc-aperture_half):min(Ny, yc+aperture_half)) = 1;
sensor.record = {'p'};

%% ========================================================================
%  Section 2: k-Wave calibration using a single reference forward simulation
%  ========================================================================

fprintf('[1/5] Running k-Wave reference calibration...\n');
source = struct(); source.p0 = single(source_mask);
input_args = {'PMLSize', PML, 'PlotPML', false, 'PMLAlpha', 2, 'Smooth', true, 'PlotSim', false, 'DataCast', 'single'};
[~, sensor_data] = evalc('kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});');

if isstruct(sensor_data) && isfield(sensor_data,'p')
    P = double(sensor_data.p);
else
    P = double(sensor_data);
end

if isvector(P)
    ref_wave = P(:)';
else
    ref_wave = mean(P, 1);
end
pp0 = max(ref_wave)-min(ref_wave);    

% --- Acoustic sensitivity weighting map A ---
r6_pix = 10;
sigma_pix = r6_pix/sqrt(2*log(2));
[X, Y] = meshgrid(1:Ny, 1:Nx);
A = exp(-(((X-yc).^2 + (Y-xc).^2)/(2*sigma_pix^2))); 
A = A/max(A(:));

sumA_ref = sum(A(source_mask>0));
fprintf('[1/5] k-Wave calibration completed.\n');

%% ========================================================================
%  Section 3: Optical system and speckle model initialization
%  ========================================================================

fprintf('[2/5] Initializing optical speckle model...\n');
SLMsize   = 64;
grainSize = 2;
n_phase   = 64;
padCoord  = 1:grainSize:(SLMsize*grainSize);
base_phase = exp(1i*2*pi*rand(SLMsize));

% --- Multi-screen scattering phase model ---
n_screens = 3; corr_list = [2, 4, 8];
scatter_phase_static = ones(Nx,Ny);
for k = 1:n_screens
    phi = smooth_random_phase(Nx,Ny,corr_list(k));
    scatter_phase_static = scatter_phase_static .* exp(1i*2*pi*normalize01(phi));
end

abs_r_pix_opt = 40;
absorber = makeDisc(Nx, Ny, xc, yc, abs_r_pix_opt);
mask_abs = absorber > 0;

Delta_init = randi([0, n_phase-1], SLMsize, SLMsize);
I_init = build_intensity_from_Delta(Delta_init, base_phase, padCoord, grainSize, n_phase, Nx, Ny, scatter_phase_static);
I_init = normalize_intensity_in_mask(I_init, mask_abs);
fprintf('[2/5] Optical speckle model initialized.\n');

% --- Coefficient mapping functions ---
coef_C1 = @(I) sum(A(mask_abs).*I(mask_abs))      / sumA_ref;
coef_C2 = @(I) sum(A(mask_abs).*(I(mask_abs).^2)) / sumA_ref;

K1_of = @(C1) pp0*(Gamma0*mu_a)*C1;                                  
K2_of = @(C2) pp0*(dGamma_dT*mu_a^2/(rho*Cp))*memory*C2;             

C1_init = coef_C1(I_init); C2_init = coef_C2(I_init);
K1_init = K1_of(C1_init);  K2_init = K2_of(C2_init);

%% ========================================================================
%  Section 4: Ideal enhancement-factor estimate
%  ========================================================================
N_dof = SLMsize^2;
L_ap  = SLMsize*grainSize;           
d_speckle = Nx / L_ap;               
S_ac   = pi * (r6_pix^2);            
S_sp   = d_speckle^2;                
M_est  = S_ac / max(S_sp, eps);      

Glin_ideal    = (pi/4)*((N_dof-1)/max(M_est,eps)) + 1;
Gnonlin_ideal = M_est;               
Gtotal_ideal  = Glin_ideal * Gnonlin_ideal;


%% ========================================================================
%  Section 5: Stage-1 linear GA optimization
%  Objective: maximize C1
%  ========================================================================

fprintf('[3/5] Running stage-1 linear GA and demo nonlinear GA...\n');
popsize = 40;
n_gen_linear    = 100;
n_gen_nonlinear = 100;
pc = 0.7;
pm = 0.02;

fitness_linear = @(Delta_flat) fitness_C1(Delta_flat, base_phase, padCoord, grainSize, n_phase, Nx, Ny, scatter_phase_static, mask_abs, A);
[Delta_lin_flat, track_lin] = ga_optimize_discrete(fitness_linear, SLMsize^2, popsize, n_gen_linear, pc, pm, n_phase, Delta_init(:));

Delta_lin = reshape(Delta_lin_flat, SLMsize, SLMsize);
I_lin  = build_intensity_from_Delta(Delta_lin,  base_phase, padCoord, grainSize, n_phase, Nx, Ny, scatter_phase_static);
I_lin  = normalize_intensity_in_mask(I_lin,  mask_abs);
C1_lin = coef_C1(I_lin);  C2_lin = coef_C2(I_lin);
K1_lin = K1_of(C1_lin);   K2_lin = K2_of(C2_lin);

% --- Single nonlinear optimization used only for waveform demonstration ---
fitness_nonlinear_demo = @(Delta_flat) fitness_C2_pure(Delta_flat, base_phase, padCoord, grainSize, n_phase, Nx, Ny, scatter_phase_static, mask_abs, A);
[Delta_non_flat, track_non] = ga_optimize_discrete(fitness_nonlinear_demo, SLMsize^2, popsize, n_gen_nonlinear, pc, pm, n_phase, Delta_lin(:));
Delta_non = reshape(Delta_non_flat, SLMsize, SLMsize);
I_non  = build_intensity_from_Delta(Delta_non, base_phase, padCoord, grainSize, n_phase, Nx, Ny, scatter_phase_static);
I_non  = normalize_intensity_in_mask(I_non, mask_abs);
C1_non = coef_C1(I_non);  C2_non = coef_C2(I_non);
K1_non = K1_of(C1_non);   K2_non = K2_of(C2_non);
fprintf('[3/5] Stage-1/demo GA completed.\n');

%% ========================================================================
%  Section 6: Dual-pulse waveform generation based on the demo optimization
%  ========================================================================
V1_init_list = K1_init*F_list; dV_init_list = K2_init*(F_list.^2);
V1_lin_list  = K1_lin *F_list; dV_lin_list  = K2_lin *(F_list.^2);
V1_non_list  = K1_non*F_list;  dV_non_list  = K2_non*(F_list.^2);

A1_init = K1_init*F_plot;  A2_init = A1_init + K2_init*(F_plot^2);
A1_lin  = K1_lin *F_plot;  A2_lin  = A1_lin  + K2_lin *(F_plot^2);
A1_non  = K1_non*F_plot;   A2_non  = A1_non  + K2_non*(F_plot^2);

nshift1  = round(t1/dt);
nshift12 = round((t1+Delta_t)/dt);
Nw = numel(ref_wave); Nlong = nshift12 + Nw;
mk_pulse = @(A, nshift) [zeros(1, nshift), A*ref_wave];

p1_i = mk_pulse(A1_init, nshift1);  p2_i = mk_pulse(A2_init, nshift12);
p1_l = mk_pulse(A1_lin, nshift1);   p2_l = mk_pulse(A2_lin, nshift12);
p1_n = mk_pulse(A1_non, nshift1);   p2_n = mk_pulse(A2_non, nshift12);

pad = @(x) [x, zeros(1, Nlong-numel(x))];
p1_i = pad(p1_i); p2_i = pad(p2_i);
p1_l = pad(p1_l); p2_l = pad(p2_l);
p1_n = pad(p1_n); p2_n = pad(p2_n);

t_long_us = (0:Nlong-1)*dt*1e6;
pp = @(x) max(x)-min(x);

V1_init = pp(p1_i); V2_init = pp(p2_i);
V1_lin  = pp(p1_l); V2_lin  = pp(p2_l);
V1_non  = pp(p1_n); V2_non  = pp(p2_n);
DeltaV_init = V2_init - V1_init;
DeltaV_lin  = V2_lin  - V1_lin;
DeltaV_non  = V2_non  - V1_non;

%% ========================================================================
%  Section 7: Stage-2 nonlinear GA parameter scan
%  Core setting: maximize noisy DeltaV feedback
%  ========================================================================

fprintf('[4/5] Running per-case nonlinear GA scan...\n');
objective = 'dV';   % objective: maximize DeltaV

beta_list = [0.03 0.10 0.30 1.00];
MPE_mJcm2 = 20;
zones = {'eta < 1%', 'eta: 1-5%', 'eta > 5% (<=MPE)'}; 

F1_ref = 1.46; 
F5_ref = 7.29; 
F_col1 = 0.80 * F1_ref;              
F_col2 = 0.50 * (F1_ref + F5_ref);   
F_col3 = min(1.25 * F5_ref, 0.90 * MPE_mJcm2); 
F_values_per_zone = [F_col1, F_col2, F_col3];


mu_base_cm = 400; mu_m = mu_base_cm*100;
SLMsize_case = SLMsize;
popsize_case = popsize; 
pc_case = pc; 
pm_case = pm;
n_gen_case   = n_gen_nonlinear;

nB = numel(beta_list);
I_non_cases = cell(nB,3);
V1_non_mat  = zeros(nB,3);
dV_non_mat  = zeros(nB,3);
ratio_non_mat = zeros(nB,3);
track_cases = cell(nB, 3);


for ib = 1:nB
    beta = beta_list(ib);
    mem_i = exp(-beta);
    Dt_us = beta * (tau_th*1e6);
    
    for zj = 1:3
        F_mJ   = F_values_per_zone(zj);
        F_Jm2  = F_mJ * 10;
        
        % --- Use noisy DeltaV as the GA fitness function ---
        fitness_case = @(Delta_flat) fitness_dV_noisy( ...
            Delta_flat, objective, F_Jm2, mem_i, mu_m, ...
            base_phase, padCoord, grainSize, n_phase, Nx, Ny, ...
            scatter_phase_static, mask_abs, A, ...
            pp0, Gamma0, dGamma_dT, rho, Cp);

        % Use the linear optimum as the seed for each nonlinear GA run
        [Delta_non_case_flat, track_case] = ga_optimize_discrete( ...
            fitness_case, SLMsize_case^2, popsize_case, n_gen_case, ...
            pc_case, pm_case, n_phase, Delta_lin(:));

        track_cases{ib, zj} = track_case;

        Delta_non_case = reshape(Delta_non_case_flat, SLMsize_case, SLMsize_case);
        I_non_case = build_intensity_from_Delta(Delta_non_case, ...
                         base_phase, padCoord, grainSize, n_phase, Nx, Ny, scatter_phase_static);
        I_non_case = normalize_intensity_in_mask(I_non_case, mask_abs);
        I_non_cases{ib, zj} = I_non_case;
        
        C1 = sum(A(mask_abs).*I_non_case(mask_abs));
        C2 = sum(A(mask_abs).*(I_non_case(mask_abs).^2));
        K1 = pp0*(Gamma0*mu_m)*C1;
        K2 = pp0*(dGamma_dT*(mu_m^2)/(rho*Cp))*mem_i*C2;
        V1 = K1*F_Jm2;
        dV = K2*(F_Jm2^2);
        
        V1_non_mat(ib,zj)   = V1;
        dV_non_mat(ib,zj)   = dV;
        ratio_non_mat(ib,zj)= dV/max(V1,eps);
        fprintf('  beta=%.2f, F=%.2f mJ/cm2: DeltaV=%.3g kPa, DeltaV/V1=%.3f, peak I=%.1f\n', ...
            beta, F_mJ, dV/1e3, dV/max(V1,eps), max(I_non_case(:)));
        
    end
end
fprintf('[4/5] Per-case nonlinear GA scan completed.\n');

%% ========================================================================
%  Section 8: Basic plotting and analysis
%  ========================================================================

fprintf('[5/5] Generating figures...\n');

% --- Figure 1: 4x3 intensity maps ---
figI = figure('Color','w','Name','Figure 1: I_non per (beta,F) with Independent Colormaps');
set(figI, 'Position', [50, 50, 900, 900]); 
tI   = tiledlayout(figI, nB, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tI, 'Final nonlinear-optimized intensity maps (I_{non})', 'FontSize', 14);

for ib = 1:nB
    for zj = 1:3
        ax = nexttile;
        imagesc(ax, I_non_cases{ib,zj}); 
        axis(ax,'image'); axis(ax,'off');
        title_str = sprintf('\\beta = %.2f\n%s', beta_list(ib), zones{zj});
        title(ax, title_str, 'FontSize', 10, 'FontWeight', 'bold');
        colorbar(ax);
    end
end
colormap(figI, 'parula');

% --- Figure 2: Bar-plot comparison ---
fig = figure('Color','w','Name','Nonlinear effect (per-case GA)');
tl  = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl,1); hold(ax1,'on'); bar(ax1, beta_list, dV_non_mat/1e3, 'grouped');
grid(ax1,'on'); xlabel(ax1,'\beta=\Delta t/\tau_{th}'); ylabel(ax1,'\Delta V [kPa]');
title(ax1,'DeltaV after per-case nonlinear GA'); legend(ax1, zones,'Location','northwest');

ax2 = nexttile(tl,2); hold(ax2,'on'); bar(ax2, beta_list, ratio_non_mat, 'grouped');
grid(ax2,'on'); xlabel(ax2,'\beta'); ylabel(ax2,'\Delta V / V_1 [-]');
title(ax2,'Relative nonlinearity after per-case GA'); legend(ax2, zones,'Location','northeast');


%% ========================================================================
%  Section 9: Energy compensation and normalization
%  ========================================================================
R_init = 1.00; R_lin = 1.00; R_non = 1.00;   

V1_init_comp = V1_init/max(R_init,eps);
V1_lin_comp  = V1_lin /max(R_lin, eps);
V1_non_comp  = V1_non /max(R_non, eps);

V2_init_comp = V2_init/max(R_init,eps);
V2_lin_comp  = V2_lin /max(R_lin, eps);
V2_non_comp  = V2_non /max(R_non, eps);

DeltaV_init_comp = DeltaV_init/max(R_init^2,eps);
DeltaV_lin_comp  = DeltaV_lin /max(R_lin^2, eps);
DeltaV_non_comp  = DeltaV_non /max(R_non^2, eps);

V1_lin_norm  = V1_lin_comp /max(V1_init_comp,  eps);
V1_non_norm  = V1_non_comp /max(V1_init_comp,  eps);
V2_lin_norm  = V2_lin_comp /max(V2_init_comp,  eps);
V2_non_norm  = V2_non_comp /max(V2_init_comp,  eps);
DeltaV_lin_norm = DeltaV_lin_comp /max(DeltaV_init_comp, eps);
DeltaV_non_norm = DeltaV_non_comp /max(DeltaV_init_comp, eps);


p1_i_norm = (p1_i/max(R_init,eps)) / max(V1_init_comp,eps);
p1_l_norm = (p1_l/max(R_lin, eps)) / max(V1_init_comp,eps);
p1_n_norm = (p1_n/max(R_non, eps)) / max(V1_init_comp,eps);
p2_i_norm = (p2_i/max(R_init,eps)) / max(V2_init_comp,eps);
p2_l_norm = (p2_l/max(R_lin, eps)) / max(V2_init_comp,eps);
p2_n_norm = (p2_n/max(R_non, eps)) / max(V2_init_comp,eps);


%% ========================================================================
%  Section 10: Additional visualizations

% Three-panel waveform plot
fig = figure('Color','w','Name','Two-Pulse Waveforms @ F_{plot} (normalized)');
tl = tiledlayout(fig,3,1,'TileSpacing','compact','Padding','compact');
col_init=[0.1 0.1 0.1]; col_lin=[0 0.45 0.74]; col_non=[0.85 0.33 0.10];

ax1 = nexttile; hold(ax1,'on');
plot(ax1, t_long_us, p1_i_norm, '-', 'Color', col_init, 'LineWidth',1.2);
plot(ax1, t_long_us, p2_i_norm, '-', 'Color', col_init, 'LineWidth',1.2);
grid(ax1,'on'); xlim(ax1,[0 60]); ylabel(ax1,'Normalized PA [-]');
title(ax1, {sprintf('Initial: F=%.1f mJ/cm^2', F_plot_mJcm2), 'Normalized peak-to-peak: PA_1 = 1.00, PA_2 = 1.00'});

ax2 = nexttile; hold(ax2,'on');
plot(ax2, t_long_us, p1_l_norm, '-', 'Color', col_lin, 'LineWidth',1.2);
plot(ax2, t_long_us, p2_l_norm, '-', 'Color', col_lin, 'LineWidth',1.2);
grid(ax2,'on'); xlim(ax2,[0 60]); ylabel(ax2,'Normalized PA [-]');
title(ax2, {'After linear optimization', sprintf('Normalized peak-to-peak: PA_1 = %.2f, PA_2 = %.2f', V1_lin_norm, V2_lin_norm)});

ax3 = nexttile; hold(ax3,'on');
plot(ax3, t_long_us, p1_n_norm, '-', 'Color', col_non, 'LineWidth',1.2);
plot(ax3, t_long_us, p2_n_norm, '-', 'Color', col_non, 'LineWidth',1.2);
grid(ax3,'on'); xlim(ax3,[0 60]); ylabel(ax3,'Normalized PA [-]'); xlabel(ax3,'Time [\mus]');
title(ax3, {'After nonlinear optimization', sprintf('Normalized peak-to-peak: PA_1 = %.2f, PA_2 = %.2f', V1_non_norm, V2_non_norm)});

figure('Color','w','Name','Intensity maps (normalized in absorber)');
subplot(2,2,1); imagesc(I_init); axis image off; colorbar; title('I_{init}');
subplot(2,2,2); imagesc(I_lin);  axis image off; colorbar; title('I_{lin}');
subplot(2,2,3); imagesc(I_non);  axis image off; colorbar; title('I_{non}');
subplot(2,2,4); imagesc(A);      axis image off; colorbar; title('Acoustic weight A');

figure('Color','w','Name','V1/V2 Overlay (normalized)');
plot(t_long_us, p1_i_norm, 'k--','LineWidth',1.8); hold on;
plot(t_long_us, p2_i_norm, 'r--','LineWidth',1.8);
plot(t_long_us, p1_n_norm, 'b-','LineWidth',1.5);
plot(t_long_us, p2_n_norm, 'g-','LineWidth',1.5);
grid on; xlim([0,60]);
xlabel('Time [\mus]'); ylabel('Normalized PA [-]');
legend({'Initial V1','Initial V2','Optimized V1','Optimized V2'},'Location','Best');
title('V1/V2 waveform comparison after energy normalization');

if ~exist('fit_lin_init','var'), fit_lin_init = fitness_linear(Delta_init(:)'); end
if ~exist('fit_non_init','var'), fit_non_init = fitness_nonlinear_demo(Delta_lin(:)'); end

G1 = 0:numel(track_lin)-1;          
G2 = 0:numel(track_non)-1;          

EF_lin = track_lin ./ max(eps, fit_lin_init);   
EF_non = track_non ./ max(eps, fit_non_init);   

col_lin = [0 0.45 0.74];      
col_non = [0.85 0.33 0.10];   

figEF = figure('Color','w','Name','PAWS enhancement factors');
tEF   = tiledlayout(figEF,1,2,'TileSpacing','compact','Padding','compact');

axL = nexttile(tEF,1); hold(axL,'on');
plot(axL, G1, EF_lin, '-o', 'LineWidth',1.4, 'MarkerSize',4, 'Color',col_lin);
grid(axL,'on'); xlim(axL,[0 G1(end)]);
xlabel(axL,'GA generation');
ylabel(axL,'Linear enhancement G_{lin} = V_1/V_{1,0}');
title(axL,'Stage 1: linear PAWS');

axR = nexttile(tEF,2); hold(axR,'on');
plot(axR, G2, EF_non, '-o', 'LineWidth',1.4, 'MarkerSize',4, 'Color',col_non);
grid(axR,'on'); xlim(axR,[0 G2(end)]);
xlabel(axR,'GA generation');
ylabel(axR,'Nonlinear enhancement G_{non} = \Delta V/\Delta V_0');
title(axR,'Stage 2: nonlinear PAWS');

%% ========================================================================
%  Section 11: Task 1 - parameter selection
%  =======================================================================
final_V2 = V1_non_mat + dV_non_mat;
V2_flat = final_V2(:);
[sorted_V2, sorted_idx] = sort(V2_flat);

idx_worst_lin  = sorted_idx(1);
idx_best_lin   = sorted_idx(end);
median_pos = 6;
idx_med_lin = sorted_idx(median_pos);

[i_worst, j_worst] = ind2sub(size(final_V2), idx_worst_lin);
[i_med,   j_med]   = ind2sub(size(final_V2), idx_med_lin);
[i_best,  j_best]  = ind2sub(size(final_V2), idx_best_lin);

params_to_plot = struct(...
    'best',  struct('i', i_best,  'j', j_best,  'tag', 'Best case'), ...
    'med',   struct('i', i_med,   'j', j_med,   'tag', 'Median case'), ...
    'worst', struct('i', i_worst, 'j', j_worst, 'tag', 'Worst case') ...
);


%% ========================================================================
%  Section 12: Task 2 - 3x3 comparison and heatmap
%  =======================================================================
fig_3x3_independent = figure('Color','w','Name','Optimization comparison with independent colorbars');
t_3x3_independent = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
set(fig_3x3_independent, 'Position', [50, 50, 1000, 800]);

plot_order = {'best', 'med', 'worst'};
row_titles = { ...
    sprintf('Best case (\\beta=%.2f, %s)', beta_list(i_best), zones{j_best}), ...
    sprintf('Median case (\\beta=%.2f, %s)', beta_list(i_med), zones{j_med}), ...
    sprintf('Worst case (\\beta=%.2f, %s)', beta_list(i_worst), zones{j_worst}) ...
};

for i = 1:3
    case_key = plot_order{i};
    case_indices = params_to_plot.(case_key);
    
    ax1 = nexttile(t_3x3_independent);
    imagesc(ax1, I_init); axis image off; colorbar(ax1); 
    if i == 1, title('Initial (I_{init})', 'FontSize', 12); end
    ylabel(ax1, row_titles{i}, 'FontWeight', 'bold', 'FontSize', 11);
    
    ax2 = nexttile(t_3x3_independent);
    imagesc(ax2, I_lin); axis image off; colorbar(ax2); 
    if i == 1, title('Linear optimization (I_{lin})', 'FontSize', 12); end

    ax3 = nexttile(t_3x3_independent);
    I_non_current = I_non_cases{case_indices.i, case_indices.j};
    imagesc(ax3, I_non_current); axis image off; colorbar(ax3); 
    if i == 1, title('Nonlinear optimization (I_{non})', 'FontSize', 12); end
end
colormap(fig_3x3_independent, 'parula');

beta_labels_heatmap = cellstr(num2str(beta_list(:), 'beta = %.2f'));
fluence_labels_heatmap = zones;

figure('Color', 'w', 'Name', 'Parameter optimization heatmap');
h = heatmap(fluence_labels_heatmap, beta_labels_heatmap, dV_non_mat / 1e3); 
h.Title = 'Nonlinear signal gain (\DeltaV) [kPa]';
h.XLabel = 'Fluence regime';
h.YLabel = 'Nonlinear parameter beta';
h.Colormap = hot; h.CellLabelFormat = '%.3g'; h.FontSize = 11;

%% ========================================================================
%  Section 13: Task 3 - convergence curves
%  =======================================================================
track_best = track_cases{i_best, j_best};
track_med  = track_cases{i_med, j_med};
track_worst= track_cases{i_worst, j_worst};

row_titles = { ...
    sprintf('Best case (\\beta=%.2f, %s)', beta_list(i_best), zones{j_best}), ...
    sprintf('Median case (\\beta=%.2f, %s)', beta_list(i_med), zones{j_med}), ...
    sprintf('Worst case (\\beta=%.2f, %s)', beta_list(i_worst), zones{j_worst}) ...
};

fig_conv_log = figure('Color', 'w', 'Name', 'Convergence comparison: log-scale y-axis');
hold on; grid on; box on;
semilogy(track_best, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', row_titles{1});
semilogy(track_med, 'g-s', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', row_titles{2});
semilogy(track_worst,'b-^', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', row_titles{3});
xlabel('GA generation'); ylabel('Fitness [log scale]');
title('Convergence comparison with log-scale y-axis'); legend('show', 'Location', 'southeast');
set(gca, 'FontSize', 11);

fig_conv_norm = figure('Color', 'w', 'Name', 'Convergence comparison: normalized fitness');
hold on; grid on; box on;
track_best_norm = track_best / max(track_best(1), eps);
track_med_norm  = track_med  / max(track_med(1), eps);
track_worst_norm= track_worst/ max(track_worst(1), eps);
plot(track_best_norm, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', row_titles{1});
plot(track_med_norm, 'g-s', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', row_titles{2});
plot(track_worst_norm,'b-^', 'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', row_titles{3});
xlabel('GA generation'); ylabel('Normalized fitness');
title('Convergence comparison by enhancement factor'); legend('show', 'Location', 'northwest');
set(gca, 'FontSize', 11);
ylim([0.95, max(track_best_norm)*1.05]); 

fig_conv_subplots = figure('Color', 'w', 'Name', 'Convergence comparison: separate panels');
tl = tiledlayout(3, 1, 'TileSpacing','compact'); 
ax1 = nexttile; plot(ax1, track_best, 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4);
grid on; ylabel(ax1, 'Fitness'); title(ax1, row_titles{1}); set(ax1, 'FontSize', 10);
ax2 = nexttile; plot(ax2, track_med, 'g-s', 'LineWidth', 1.5, 'MarkerSize', 4);
grid on; ylabel(ax2, 'Fitness'); title(ax2, row_titles{2}); set(ax2, 'FontSize', 10);
ax3 = nexttile; plot(ax3, track_worst, 'b-^', 'LineWidth', 1.5, 'MarkerSize', 4);
grid on; ylabel(ax3, 'Fitness'); title(ax3, row_titles{3}); xlabel(ax3, 'GA generation'); set(ax3, 'FontSize', 10);

%% ========================================================================
%  Section 15: Task 5 - line-profile comparison
%  =======================================================================
fig_profile = figure('Color', 'w', 'Name', 'Line profiles of selected focal spots');
hold on; grid on; box on;

I_best  = I_non_cases{i_best, j_best};
I_med   = I_non_cases{i_med, j_med};
I_worst = I_non_cases{i_worst, j_worst};

profile_x_axis = ((1:Ny) - xc) * dx * 1e3; 
center_row = xc;

plot(profile_x_axis, I_best(center_row, :), 'r-', 'LineWidth', 1.5, 'DisplayName', row_titles{1});
plot(profile_x_axis, I_med(center_row, :), 'g-', 'LineWidth', 1.5, 'DisplayName', row_titles{2});
plot(profile_x_axis, I_worst(center_row,:), 'b-', 'LineWidth', 1.5, 'DisplayName', row_titles{3});

xlabel('Lateral position (mm)'); ylabel('Normalized intensity');
title('Central line profiles after nonlinear optimization'); legend('show', 'Location', 'northeast');
xlim([min(profile_x_axis), max(profile_x_axis)]);
set(gca, 'FontSize', 11);

%% ========================================================================
%  Section 16: Additional analysis for three representative regimes
%  =======================================================================
[~, idx_best_lin] = max(final_V2(:));
[i_best, j_best] = ind2sub(size(final_V2), idx_best_lin);

[~, j_high_beta] = max(final_V2(nB, :));
i_high_beta = nB;

[~, i_low_F] = max(final_V2(:, 1));
j_low_F = 1;

cases_to_plot = {
    struct('i', i_best, 'j', j_best, 'label', sprintf('Global best (\\beta=%.2f, %s)', beta_list(i_best), zones{j_best})), ...
    struct('i', i_high_beta, 'j', j_high_beta, 'label', sprintf('Sub-optimal: thermal mismatch (\\beta=%.2f, %s)', beta_list(i_high_beta), zones{j_high_beta})), ...
    struct('i', i_low_F, 'j', j_low_F, 'label', sprintf('Sub-optimal: weak nonlinearity (\\beta=%.2f, %s)', beta_list(i_low_F), zones{j_low_F}))
};

fig_compare_specific = figure('Color','w','Name','Representative-regime optimization comparison');
set(fig_compare_specific, 'Position', [100, 100, 1000, 800]);
t_compare = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:3
    current_case = cases_to_plot{i};
    ax1 = nexttile; imagesc(ax1, I_init); axis image off; colorbar(ax1);
    if i == 1, title('Initial (I_{init})', 'FontSize', 12); end
    ylabel(ax1, current_case.label, 'FontWeight', 'bold', 'FontSize', 11);
    
    ax2 = nexttile; imagesc(ax2, I_lin); axis image off; colorbar(ax2);
    if i == 1, title('Linear optimization (I_{lin})', 'FontSize', 12); end
    
    ax3 = nexttile; I_non_current = I_non_cases{current_case.i, current_case.j};
    imagesc(ax3, I_non_current); axis image off; colorbar(ax3);
    if i == 1, title('Nonlinear optimization (I_{non})', 'FontSize', 12); end
end
colormap(fig_compare_specific, 'parula');

% Figures 14--16 were supplementary/redundant display panels and were removed in this GitHub version:
% - normalized convergence curve for selected English cases
% - selected-case profile comparison
% - 2x2 independent colorbar panel

fprintf('[5/5] Figures generated. Simulation finished.\n');

%% ========================================================================
%  Section 17: Local helper functions including noisy-fitness model
%  =======================================================================

function fit = fitness_dV_noisy(Delta_flat, objective, F_Jm2, mem, mu_m, ...
        base_phase, padCoord, grainSize, n_phase, Nx, Ny, scatter_phase, mask, A, ...
        pp0, Gamma0, dGamma_dT, rho, Cp)
    % fitness_dV_noisy (Based on GR_PAWS_V2.m fitness_non_betaF_noisy)
    % Objective: calculate noisy DeltaV as the GA feedback signal.
    
    % 1. Compute the noise-free physical amplitudes V1_true and dV_true.
    SLMsize = sqrt(numel(Delta_flat));
    Delta   = reshape(Delta_flat, SLMsize, SLMsize);

    I = build_intensity_from_Delta(Delta, base_phase, padCoord, grainSize, n_phase, Nx, Ny, scatter_phase);
    I = normalize_intensity_in_mask(I, mask);

    C1 = sum(A(mask) .* I(mask));
    C2 = sum(A(mask) .* (I(mask).^2));

    K1 = pp0 * (Gamma0    * mu_m)                 * C1;
    K2 = pp0 * (dGamma_dT * (mu_m^2)/(rho*Cp))    * mem * C2;

    V1_true = K1 * F_Jm2;
    dV_true = K2 * (F_Jm2^2);

    % 2. Simulate measurement noise and finite averaging.
    % Noise parameters; tune as needed for the experiment.
    N_avg      = 8;      % number of averages
    sigma_elec = 2e3;    % electronic noise [Pa]
    sigma_rel  = 0.03;   % relative noise (3%)

    dV_meas_sum = 0;
    V1_meas_sum = 0;
    
    % helper noise function
    add_noise = @(val) val + sqrt(sigma_elec^2 + (sigma_rel * max(abs(val),1)).^2) * randn();

    for n = 1:N_avg
        % simulate separate measurements of V1 and V2, where V2 = V1 + dV
        V1_meas_n = add_noise(V1_true);
        V2_meas_n = add_noise(V1_true + dV_true);
        dV_meas_n = V2_meas_n - V1_meas_n;

        dV_meas_sum = dV_meas_sum + dV_meas_n;
        V1_meas_sum = V1_meas_sum + V1_meas_n;
    end

    dV_meas_avg = dV_meas_sum / N_avg;
    V1_meas_avg = V1_meas_sum / N_avg;

    % 3. Return fitness.
    switch objective
        case 'ratio'
            fit = dV_meas_avg / max(V1_meas_avg, eps);
        otherwise
            fit = dV_meas_avg; % default: maximize DeltaV
    end
end

function I = build_intensity_from_Delta(Delta, base_phase, padCoord, grainSize, n_phase, Nx, Ny, scatter_phase)
    SLMsize = size(Delta,1);
    pad = zeros(Nx, Ny);
    phase_mask = base_phase .* exp(1i*2*pi*Delta/n_phase);
    for x = 1:SLMsize
        for y = 1:SLMsize
            pad(padCoord(x):padCoord(x)+grainSize-1, padCoord(y):padCoord(y)+grainSize-1) = phase_mask(x,y);
        end
    end
    I = generate_scattered_intensity(pad, scatter_phase);
end

function I_norm = normalize_intensity_in_mask(I, mask)
    I = max(I,0);
    m = mean(I(mask));
    if ~isfinite(m) || m<=0, I_norm = zeros(size(I)); else, I_norm = I/m; end
end

function fit = fitness_C1(Delta_flat, base_phase, padCoord, grainSize, n_phase, Nx, Ny, scatter_phase, mask, A)
    SLMsize = sqrt(numel(Delta_flat));
    Delta = reshape(Delta_flat, SLMsize, SLMsize);
    I = build_intensity_from_Delta(Delta, base_phase, padCoord, grainSize, n_phase, Nx, Ny, scatter_phase);
    I = normalize_intensity_in_mask(I, mask);
    fit = sum(A(mask) .* I(mask));                        
end

function fit = fitness_C2_pure(Delta_flat, base_phase, padCoord, grainSize, n_phase, Nx, Ny, scatter_phase, mask, A)
    SLMsize = sqrt(numel(Delta_flat));
    Delta = reshape(Delta_flat, SLMsize, SLMsize);
    I = build_intensity_from_Delta(Delta, base_phase, padCoord, grainSize, n_phase, Nx, Ny, scatter_phase);
    I = normalize_intensity_in_mask(I, mask);
    fit = sum(A(mask) .* (I(mask).^2));                   
end

function [best_flat, track] = ga_optimize_discrete(fitnessFcn, chromlength, popsize, n_gen, pc, pm, n_phase, seed_flat)
    if nargin < 8 || isempty(seed_flat), seed_flat = randi([0, n_phase-1], chromlength, 1); end
    pop = randi([0, n_phase-1], popsize, chromlength);
    pop(1,:) = seed_flat(:)';         
    fit = zeros(popsize,1);
    best_flat = pop(1,:)'; best_fit = -inf;
    track = zeros(n_gen+1,1);
    for i=1:popsize
        fit(i) = fitnessFcn(pop(i,:));
        if fit(i)>best_fit, best_fit=fit(i); best_flat=pop(i,:)'; end
    end
    track(1) = best_fit;
    for gen=1:n_gen
        newpop = tournament_selection(pop, fit, 3);
        newpop = single_point_crossover(newpop, pc);
        newpop = random_mutation(newpop, pm, n_phase);
        newfit = zeros(popsize,1);
        for i=1:popsize
            newfit(i) = fitnessFcn(newpop(i,:));
            if newfit(i)>best_fit, best_fit=newfit(i); best_flat=newpop(i,:)'; end
        end
        [~, worst_idx] = min(newfit);
        newpop(worst_idx,:) = best_flat(:)'; newfit(worst_idx)=best_fit;
        pop=newpop; fit=newfit; track(gen+1)=best_fit;
    end
end

function newpop = tournament_selection(pop, fit, k)
    [popsize, ~] = size(pop); newpop = zeros(size(pop));
    for i=1:popsize
        idx = randi(popsize,[k,1]); [~,ib] = max(fit(idx)); newpop(i,:) = pop(idx(ib),:);
    end
end

function newpop = single_point_crossover(pop, pc)
    [popsize, L] = size(pop); newpop = pop;
    for i=1:2:popsize-1
        if rand < pc
            c = randi([1,L-1]); a=newpop(i,:); b=newpop(i+1,:);
            newpop(i,:)=[a(1:c), b(c+1:end)]; newpop(i+1,:)=[b(1:c), a(c+1:end)];
        end
    end
end

function newpop = random_mutation(pop, pm, n_phase)
    [popsize, L] = size(pop); newpop = pop;
    mask = rand(popsize,L) < pm;
    randVals = randi([0,n_phase-1], popsize, L);
    newpop(mask) = randVals(mask);
end

function I_scat = generate_scattered_intensity(pad, scatter_phase)
    pad = double(pad); scatter_phase = double(scatter_phase);
    U = fft2(fftshift(pad .* scatter_phase));
    I_raw = abs(fftshift(U)).^2;
    m = max(I_raw(:));
    if ~isfinite(m) || m<=0, I_scat = zeros(size(I_raw)); else, I_scat = I_raw/m; end
end

function phi = smooth_random_phase(Nx,Ny,sigma_pix)
    z = randn(Nx,Ny);
    gsz = max(3, 2*ceil(3*sigma_pix)+1);
    g = fspecial('gaussian', [gsz gsz], sigma_pix);
    z = conv2(z, g, 'same');
    phi = z / max(eps, std(z(:)));
end

function y = normalize01(x)
    x = x - min(x(:)); mx = max(x(:));
    y = x / max(mx, eps);
end
