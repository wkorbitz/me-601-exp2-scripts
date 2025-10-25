%% Cyclic Voltammetry Analysis with Baseline Correction, Peak Detection, and Theoretical Comparison
clear; clc; close all;

%% === USER INPUTS ===
filename = 'exp2-lm/CV_100mV_Pt_noMem.csv';
skipRows = 0;
%delimiter = ',';
scanRate = 0.1;  % V/s
n = 4; %number of electrons
d = 0.3; %electrode diameter [cm]
A = pi()*d^2/4; %electrode area [cm2]
C = 10e-6; %species bulk concentration [mol/cm^3] (1 mM)
D = 7.6e-6; % diffusion coefficient [cm^2/s] ferricyanide

%% === READ RAW DATA ===
raw = readmatrix(filename,'NumHeaderLines',skipRows);
E = raw(:,1); I_raw = raw(:,2);

%% === ROBUST BASELINE FIT ===
absI = abs(I_raw);
threshold = prctile(absI,20);
mask = absI < threshold;
p = polyfit(E(mask),I_raw(mask),1);
I_baseline = polyval(p,E);
I_corr = I_raw - I_baseline;

%% === PEAK DETECTION ===
[Ip_a, idx_a] = max(I_corr);
[Ip_c, idx_c] = min(I_corr);
Ep_a = E(idx_a); Ep_c = E(idx_c);
E_half = (Ep_a + Ep_c)/2;
dEp = Ep_a - Ep_c;
Ip_ratio = abs(Ip_a/Ip_c);

%% === THEORETICAL PEAK CURRENT (Randles–Ševčík) ===
Ip_theory = 2.69e5 * n^(1.5) * A * C * sqrt(D) * sqrt(scanRate); % A
Ip_theory_uA = Ip_theory * 1e6;

%% === DOUBLE-LAYER CAPACITANCE ===
Icap_mean = mean(I_baseline(mask)); % mean baseline current in µA
Cdl = abs(Icap_mean * 1e-6) / scanRate; % ensure positive
Cdl_uF = Cdl * 1e6;

percents = [0.01 0.05 0.10];
Icap_targets = percents * Ip_theory_uA;
Cdl_targets_uF = Icap_targets / scanRate; % since Cdl(µF)=10*Icap(µA) at v=0.1

%% === PLOTS ===
figure('Position',[100 100 900 600]);

% --- Plot raw CV and baseline ---
%subplot(2,1,1);
plot(E,I_raw,'k','DisplayName','Raw CV'); hold on;
plot(E,I_baseline,'r--','DisplayName','Baseline fit');
xlabel('Potential [V]'); ylabel('Current [µA]');
title('Raw CV with Baseline');
legend('Location','best'); grid on;

% Annotate baseline level
text(mean(E),mean(I_baseline),sprintf('Baseline ≈ %.2f µA',Icap_mean),...
    'VerticalAlignment','bottom','HorizontalAlignment','center','Color','r');

% --- Plot corrected CV with experimental + theoretical peaks ---
%subplot(2,1,2);
%plot(E,I_corr,'b','DisplayName','Baseline-corrected'); hold on;

% Experimental peaks
plot(E(idx_a),Ip_a,'ro','MarkerFaceColor','r','DisplayName','Exp anodic peak');
plot(E(idx_c),Ip_c,'go','MarkerFaceColor','g','DisplayName','Exp cathodic peak');
xline(E_half,'k--','E_{1/2}','LabelVerticalAlignment','bottom');

% Theoretical peak markers (centered on experimental Ep positions)
%plot(Ep_a,Ip_theory_uA,'r^','MarkerFaceColor','none','LineWidth',1.5,...
%    'DisplayName','Theory Ipa');
%plot(Ep_c,-Ip_theory_uA,'g^','MarkerFaceColor','none','LineWidth',1.5,...
%    'DisplayName','Theory Ipc');

% Annotate peaks (using TeX-compatible \newline)
text(Ep_a,Ip_a,sprintf('Epa=%.3f V\\newline Ipa=%.2f µA',Ep_a,Ip_a), ...
    'VerticalAlignment','bottom','HorizontalAlignment','left','Color','r');

text(Ep_c,Ip_c,sprintf('Epc=%.3f V\\newline Ipc=%.2f µA',Ep_c,Ip_c), ...
    'VerticalAlignment','top','HorizontalAlignment','left','Color','g');

xlabel('Potential [V]'); ylabel('Current [µA]');
title('Baseline-Corrected CV: Experimental vs. Theoretical Peaks');
legend('Location','best'); grid on;

%% === RESULTS OUTPUT ===
fprintf('\n=== Experimental Results ===\n');
fprintf(' Ep_a = %.4f V, Ip_a = %.3f µA\n',Ep_a,Ip_a);
fprintf(' Ep_c = %.4f V, Ip_c = %.3f µA\n',Ep_c,Ip_c);
fprintf(' ΔEp = %.4f V (theory: %.4f V)\n',dEp,0.0566/n);
fprintf(' E1/2 = %.4f V\n',E_half);
fprintf(' |Ip_a/Ip_c| = %.3f (ideal=1)\n',Ip_ratio);
fprintf('\n=== Theoretical Peak Current ===\n');
fprintf(' Ip_theory = %.3f µA (Randles-Sevcik)\n',Ip_theory_uA);
fprintf('\n=== Double-Layer Capacitance ===\n');
fprintf(' Estimated Cdl = %.3f µF\n',Cdl_uF);
for i=1:length(percents)
    fprintf(' Cdl for %.0f%% of Ip_theory: %.3f µF\n',percents(i)*100,Cdl_targets_uF(i));
end
