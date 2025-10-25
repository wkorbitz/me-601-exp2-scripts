%% ECSA Calculation Script for Pt and Pd Electrodes
% Using CV data in acidic electrolyte
% Reference: Ag/AgCl (3M KCl)
% Includes baseline correction and conversion to RHE

clear; clc; close all

%% ======== USER INPUTS =========
filename = ['exp2-jc\10OCT25_PtBlackCatalyst_CV_OER_50mVs.csv'];        % 2-column file: Potential (V vs Ag/AgCl), Current (A)
metal = 'Pt';                    % Choose: 'Pt' or 'Pd'
A_geo = 0.196;                   % Geometric electrode area [cm2]
scan_rate = 50e-3;               % Scan rate [V/s] (e.g., 50 mV/s)
potential_range = [0.8, 1.3];  % Integration window for H-desorption region [V vs RHE]
pH = 0;                          % pH of electrolyte (0 for 0.5 M H2SO4, etc.)

%% ======== CONSTANTS ===========
E_AgAgCl_to_RHE = 0.210 + 0.0591 * pH;  % Conversion [V]
% -> E(RHE) = E(Ag/AgCl) + 0.197 + 0.0591*pH

%% ======== LOAD DATA ===========
data = readmatrix(filename);
E_AgAgCl = data(:,1);
E_AgAgCl = E_AgAgCl(2:numel(E_AgAgCl)-2);

I = data(:,2) * 1e-6;  % Current in A
I = I(2:numel(I)-2);

% Convert potential to RHE
E = E_AgAgCl + E_AgAgCl_to_RHE;

%% ======== DEFINE INTEGRATION RANGE ==========
idx = (E >= potential_range(1)) & (E <= potential_range(2));

% ======== BASELINE CORRECTION (Fixed) ==========
% Identify numeric indices of integration window
idx_all = find(idx);
i_start = idx_all(1);
i_end = idx_all(end);

% Define baseline endpoints
E_base = [E(i_start), E(i_end)];
I_base = [I(i_start), I(i_end)];

% Interpolate linear baseline across integration range
baseline = interp1(E_base, I_base, E(idx), 'linear', 'extrap');

% Subtract baseline
I_corr = I(idx) - baseline;

%% ======== INTEGRATE CHARGE ==========
Q_H = trapz(E(idx), I_corr) / scan_rate;  % Charge in Coulombs

%% ======== REFERENCE CHARGE DENSITY ==========
switch metal
    case 'Pt'
        Q_ref = 2.10e-4; % C/cm2 for H adsorption monolayer on Pt
    case 'Pd'
        Q_ref = 4.20e-4; % C/cm2 for H adsorption monolayer on Pd
    otherwise
        error('Invalid metal type. Choose ''Pt'' or ''Pd''.')
end

%% ======== CALCULATE ECSA AND ROUGHNESS ==========
ECSA = Q_H / Q_ref;        % [cm2]
Roughness = ECSA / A_geo;  % Dimensionless

%% ======== DISPLAY RESULTS ==========
fprintf('\n=== ECSA Analysis for %s Electrode ===\n', metal)
fprintf('Reference: Ag/AgCl (converted to RHE)\n')
fprintf('Integrated H-desorption charge (Q_H): %.3e C\n', Q_H)
fprintf('Reference charge density (Q_ref): %.2e C/cm²\n', Q_ref)
fprintf('ECSA: %.2f cm²\n', ECSA)
fprintf('ECSA / geometric area: %.2f (dimensionless)\n', Roughness)

%% ======== PLOTS ==========
figure('Position',[100 100 700 500]);

% Plot raw CV
subplot(2,1,1)
plot(E, I*1e3, 'k', 'LineWidth', 1.5)
xlabel('Potential (V vs RHE)')
ylabel('Current (mA)')
title(sprintf('Raw CV (%s Electrode)', metal))
grid on

% Plot corrected region
subplot(2,1,2)
plot(E(idx), I(idx)*1e3, 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2); hold on
plot(E(idx), baseline*1e3, 'r--', 'LineWidth', 1)
area(E(idx), I_corr*1e3, 'FaceColor', [0.8 0.8 1], 'EdgeColor', 'none')
xlabel('Potential (V vs RHE)')
ylabel('Current (mA)')
title('Baseline-Corrected Hydrogen Region Integration')
legend('Raw region', 'Baseline', 'Integrated (corrected)', 'Location', 'best')
grid on
