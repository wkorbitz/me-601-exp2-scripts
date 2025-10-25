%% LSV Tafel Analysis Script
% Analyze raw LSV data for ferri/ferrocyanide redox couple
% Goal: Calculate experimental Tafel slope and compare to theoretical value

clear; clc; close all;

%% User Inputs
filename = 'exp2\exp2-jc\10OCT25_PtBlackCatalyst_LSV_HER_5mVs.csv'; % Your raw data file: [Voltage (V), Current (A)]
electrode_area = 0.196e-4; % cm^2 (example: 5 mm diameter disk)
n = 4;                     % number of electrons (ferri/ferrocyanide = 1)
T = 298.15;                % temperature (K)
alpha = 0.5;               % charge transfer coefficient (assume 0.5)
R = 8.314;                 % J/mol/K
F = 96485;                 % C/mol

%% Load Data
data = readmatrix(filename);
E = data(:,1); % V
E = E(2:numel(E)-2);
I = data(:,2); % A
I = I(2:numel(I)-2);

% Convert to current density (A/cm^2)
j = I ./ electrode_area;

%% Calculate Overpotential
% Determine equilibrium potential (approx. by zero current potential)
E_eq = interp1(j(16:17), E(16:17), 0, 'linear', 'extrap'); 
eta = E - E_eq; % overpotential

%% Filter Data for Tafel Region
% Keep only data where current is positive (anodic branch)
idx = j > 0;
eta = eta(idx);
j = j(idx);

% Log current density (semi-log for Tafel)
log_j = log10(j);

% Select linear region manually or automatically
% Option 1: Manual selection (recommended)
figure;
plot(eta, log_j, 'o'); grid on;
xlabel('Overpotential (V)'); ylabel('log_{10}(j) [A/cm^2]');
title('Select Tafel Region: Click on lower and upper bounds');
[xsel, ~] = ginput(2);
mask = (eta >= min(xsel)) & (eta <= max(xsel));

eta_fit = eta(mask);
log_j_fit = log_j(mask);

%% Linear Fit for Tafel Slope
p = polyfit(eta_fit, log_j_fit, 1);
b_exp = 1/p(1); % experimental Tafel slope (V/dec)
j0_exp = 10^(p(2)); % exchange current density

%% Theoretical Tafel Slope
b_theo = 2.303 * R * T / (alpha * n * F); % V/dec

%% Plot Results
figure;
plot(eta, log_j, 'ko','MarkerSize',4); hold on;
plot(eta_fit, polyval(p, eta_fit), 'r-', 'LineWidth', 2);
xlabel('Overpotential, \eta (V)');
ylabel('log_{10}(j) [A/cm^2]');
title('Tafel Plot');
legend('Data','Linear Fit','Location','best');
grid on;

% Annotate results
text(min(eta_fit), max(log_j_fit), sprintf('Experimental Tafel Slope: %.3f V/dec', b_exp), ...
    'VerticalAlignment','top','FontSize',10,'Color','r');
text(min(eta_fit), max(log_j_fit)-0.5, sprintf('Theoretical Tafel Slope: %.3f V/dec', b_theo), ...
    'VerticalAlignment','top','FontSize',10,'Color','b');

%% Display results in Command Window
fprintf('=== Tafel Analysis Results ===\n');
fprintf('Equilibrium Potential: %.4f V\n', E_eq);
fprintf('Experimental Tafel slope: %.3f V/dec\n', b_exp);
fprintf('Theoretical Tafel slope: %.3f V/dec\n', b_theo);
fprintf('Experimental Exchange Current Density: %.3e A/cm^2\n', j0_exp);
