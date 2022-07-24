% last update: 2022/7/24
% test the performance of NMOP-CFAR-2D

clc; clear; close all;

rng(5)
MC = 5;

% Define Scenario
Nx = 64; % Length of Sinusoid
My = 32; % width of Sinusoid
K = 16;
sigma_n = 1;              % noise variance sigma^2, instead of sigma

S_snap = 1;
SNR = 12;

% algorithm parameters
K_max = 20;
P_oe = 0.01;
CFAR_method = 'CA';
gamma_oversamping = 4;
N_r = 60;





