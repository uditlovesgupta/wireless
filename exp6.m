%Okumura Model

clc; clear; close all;

fc = 900; % MHz
hte = 100; % m
hre = 10; % m
d = 10:10:50; % km
lambda = 3e8 / (fc * 1e6);
EIRP = 60; % dBm
GA = 9;

LF = 10 * log10((4 * pi * d * 1e3 ./ lambda).^2);
Ghte = 20 * log10(hte / 200);
Ghre = 20 * log10(hre / 3);

A_mu_semiopen = [30, 32, 36, 42, 45];
A_mu_open = [25, 27, 30, 35, 38];

L50_semiopen = LF + A_mu_semiopen - Ghte - Ghre - GA;
L50_open = LF + A_mu_open - Ghte - Ghre - GA;

Pr_semiopen = EIRP - L50_semiopen;
Pr_open = EIRP - L50_open;

figure;
plot(d, L50_semiopen, '-o', d, L50_open, '-s', 'LineWidth', 2);
xlabel('Distance (km)');
ylabel('Median Path Loss (dB)');
title('Okumura Path Loss'); legend('Semiopen', 'Open'); grid on;

figure;
plot(d, Pr_semiopen, '-o', d, Pr_open, '-s', 'LineWidth', 2);
xlabel('Distance (km)');
ylabel('Received Power (dBm)');
title('Okumura Received Power'); legend('Semiopen', 'Open'); grid on;

fprintf('Results at %d km:\n', d(end));
fprintf('Semi-open Path Loss: %.2f dB\n', L50_semiopen(end));
fprintf('Open Path Loss: %.2f dB\n', L50_open(end));
fprintf('Semi-open Power: %.2f dBm\n', Pr_semiopen(end));
fprintf('Open Power: %.2f dBm\n', Pr_open(end));