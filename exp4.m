%% PART A - Log-Normal Path Loss vs Distance
clc;
clear all;
close all;

d0 = 0.5;                     % Reference distance in km
d = 1:1:20;                   % Distance in km
n = 3;                        % Path loss exponent
f = 900;                      % Frequency in MHz
a = randn(1, length(d));      % Gaussian noise for log-normal shadowing

% Free-space path loss at reference distance (in dB)
fspl = 32.45 + 20*log10(d0) + 20*log10(f);

% Log-normal path loss (in dB)
PL = fspl + 10 * n * log10(d / d0) + a;

% Plotting path loss vs. distance
figure(1);
plot(d, PL, 'b-o');
title('Log-Normal Path Loss vs Distance');
xlabel('Distance (km)');
ylabel('Path Loss (dB)');
grid on;

%% PART B - Log-Normal Path Loss vs Path Loss Exponent
clear all;
close all;

d0 = 0.5;                       % Reference distance in km
d = 5;                          % Distance in km
n = 2:0.2:4;                    % Path loss exponent values
f = 900;                        % Frequency in MHz
a = randn(1, length(n));        % Gaussian noise for log-normal shadowing

% Free-space path loss at reference distance (in dB)
fspl = 32.45 + 20*log10(d0) + 20*log10(f);

% Log-normal path loss (in dB)
PL = fspl + 10 * n .* log10(d / d0) + a;

% Plotting path loss vs. path loss exponent
figure(2);
plot(n, PL, 'r-*');
title('Log-Normal Path Loss vs Path Loss Exponent');
xlabel('Path Loss Exponent (n)');
ylabel('Path Loss (dB)');
grid on;

%% PART C - Received Power vs Path Loss Exponent
d0 = 0.5;               % Reference distance in km
f = 900;                % Frequency in MHz
d = 5;                  % Distance in km
n = 2:0.2:4;            % Path loss exponent values

pt_dbm = 20;            % Transmitted power in dBm
Pt = pt_dbm - 30;       % Convert to dB (Watts)

% Free-space path loss at reference distance (in dB)
PLFS = 20 * log10(d0) + 20 * log10(f) - 147.55;

% Generate random Gaussian variable for shadowing
X0 = randn(size(n));

% Log-normal shadowing model path loss
PL = PLFS + 10 .* n .* log10(d / d0) + X0;

% Received power = Transmitted power - Path loss
Pr = Pt - PL;

% Plotting received power vs. path loss exponent
figure(3);
plot(n, Pr, 'g-s');
title('Received Power vs Path Loss Exponent');
xlabel('Path Loss Exponent (n)');
ylabel('Received Power (dB)');
grid on;
