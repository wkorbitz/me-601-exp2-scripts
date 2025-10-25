clear all; close all; clc;

filename = 'exp2-lm/CV_100mV_Pt.csv';
skipRows = 0;
%delimiter = ',';
scanRate = 0.1;  % V/s
n = 1; %number of electrons
d = 0.3; %electrode diameter [cm]
A = pi()*d^2/4; %electrode area [cm2]
C = 10e-6; %species bulk concentration [mol/cm^3] (1 mM)
D = 7.6e-6; % diffusion coefficient [cm^2/s] ferricyanide

raw = readmatrix(filename,'NumHeaderLines',skipRows);
E = raw(:,19); I_raw = raw(:,20);

plot(E, I_raw);