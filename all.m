%Exp1
% Question (a)
clc
clear all
close all

B = 30 * 10^6;
C = 50 * 10^3;
S = B / C;

T_Area = input("Enter the graphical area to be covered: ");
A_Cell = input("Enter the area of single cell: ");
total_num_cell_N = T_Area / A_Cell;
cluster_size = input("Enter the number of clusters: ");
k = S / cluster_size;
num_cluster_M = total_num_cell_N / cluster_size;
cap_c = num_cluster_M * S;

% Question (b)
c = 20:20:200;
b = input("Enter the total bandwidth (in MHz): ");
b = b * 10^6;
n = [3,4,7,9,12];
a = input("Enter the total area: ");
ca = input("Enter area of cell: ");
figure
for i = 1:length(n)
    K = (b ./ (c * 1000)) * n(i);
    no_of_cluster = a / ca;
    capacity = K * no_of_cluster;
    plot(c, capacity, 'DisplayName', ['N=' num2str(n(i))]);
    hold on;
end
xlabel("Channel Bandwidth (kHz)");
ylabel("Capacity");
title("Channel Bandwidth vs Capacity");
legend;
grid on;

% Question (c)
r = 100:100:1000;
Total_BW = 1800 * 10^6;
channel_BW = 200 * 10^3;
N = [3,4,7,9,12];
Total_Area = 2100;
figure
for i = 1:length(N)
    Cluster_Area = 2.6 * r.^2 / 10^6;
    M = Total_Area ./ Cluster_Area;
    S = N(i) * (Total_BW / channel_BW);
    C = M * S;
    plot(r, C, 'DisplayName', ['N=' num2str(N(i))]);
    hold on;
end
xlabel("Cell Radius (m)");
ylabel("Capacity");
title("Cell Radius vs Capacity");
legend;
grid on;


%Exp2
% Question (A) - Handover Probability vs Mobile Speed
clc;
clear all;
close all;
r = 500;
avg = 3 * 60;
v = 5:5:50;
u = 1 / avg;
PI = 3.14;
n = (2 .* v) / (PI * r);
Hp = n ./ (n + u);
disp(Hp)
figure(1)
plot(v, Hp)
xlabel("Velocity");
ylabel("Hand Over Probability");
title("Handover Probability vs Mobile Speed");
grid on

% Question (B) - Handover Probability vs Cell Radius
r = 100:100:1000;
avg = 3 * 60;
v = 60 * (5 / 18);
u = 1 / avg;
PI = 3.14;
n = (2 .* v) ./ (PI * r);
Hp = n ./ (n + u);
figure(2)
plot(r, Hp)
xlabel("Cell Radius");
ylabel("Handover Probability");
title("Cell Radius vs Handoff Probability");
grid on

% Question (C) - Distance Between Co-Channel Cells vs SIR
d = 500:100:2000;
r = 5 * 10^3;
n = [2, 4];
q = 6;
capacity_matrix = zeros(length(d), length(n));
for i = 1:length(d)
    for j = 1:length(n)
        s = ((d(i) / r)^n(j)) / q;
        capacity_matrix(i, j) = s;
    end
end
figure(3)
hold on
for i = 1:length(n)
    plot(d, capacity_matrix(:, i), 'DisplayName', ['n=' num2str(n(i))])
end
xlabel("Distance Between Co-Channel Cells (m)");
ylabel("SIR");
title("Distance Between Co-Channel Cells vs SIR");
legend;
grid on

% Question (D) - SIR Table for Different Sectorization Scenarios
n = 3;
N1 = 5;
N2 = 7;
N3 = 19;
i1 = 6;
i2 = 1;
i3 = 2;

SIR11 = 10 * log10(((3 * N1)^(n / 2)) / i1);
SIR12 = 10 * log10(((3 * N1)^(n / 2)) / i2);
SIR13 = 10 * log10(((3 * N1)^(n / 2)) / i3);
SIR21 = 10 * log10(((3 * N2)^(n / 2)) / i1);
SIR22 = 10 * log10(((3 * N2)^(n / 2)) / i2);
SIR23 = 10 * log10(((3 * N2)^(n / 2)) / i3);
SIR31 = 10 * log10(((3 * N3)^(n / 2)) / i1);
SIR32 = 10 * log10(((3 * N3)^(n / 2)) / i2);
SIR33 = 10 * log10(((3 * N3)^(n / 2)) / i3);

Sectorization = ["No sectorization"; "60° sectorization"; "120° sectorization";
                 "No sectorization"; "60° sectorization"; "120° sectorization";
                 "No sectorization"; "60° sectorization"; "120° sectorization"];

SIRDb = [SIR11; SIR12; SIR13; SIR21; SIR22; SIR23; SIR31; SIR32; SIR33];
N = [N1; N1; N1; N2; N2; N2; N3; N3; N3];

t1 = table(Sectorization, N, SIRDb)

%Exp3
clc;
clear all;

%% PART - A -------------------------
% Free Space Propagation Model
freq = 900e6;         % Frequency in Hz
Pt = 10;              % Transmit power in Watts
Gt = 5;               % Transmitter gain
Gr = 3;               % Receiver gain
loss = 1;             % System loss
lambda = 3e8 / freq;  % Wavelength (c/f)

d = 100:10:2000;      % Distance range in meters

% Received Power Calculation
Pr = (Pt * Gt * Gr * (lambda^2)) ./ ((d.^2) * loss * (4 * pi)^2);
figure(1)
plot(d, Pr, 'b');
xlabel("Distance (m)");
ylabel("Received Power (W)");
title("Free Space: Distance vs Received Power");
grid on

% Path Loss Calculation
path_loss = -10 * log10(Pr / Pt);
figure(2)
plot(d, path_loss, 'r');
xlabel("Distance (m)");
ylabel("Path Loss (dB)");
title("Free Space: Distance vs Path Loss");
grid on

%% PART - B -------------------------
% Two-Ray Ground Reflection Model
ht = 40;  % Transmitter height in meters
hr = 3;   % Receiver height in meters

% Received Power Calculation
Pr2ray = Pt * Gt * Gr * (ht^2) * (hr^2) ./ (d.^4);
figure(3)
plot(d, Pr2ray, 'g');
xlabel("Distance (m)");
ylabel("Received Power (W)");
title("Two-Ray: Distance vs Received Power");
grid on

% Path Loss Calculation
path_loss2ray = -10 * log10(Pr2ray / Pt);
figure(4)
plot(d, path_loss2ray, 'm');
xlabel("Distance (m)");
ylabel("Path Loss (dB)");
title("Two-Ray: Distance vs Path Loss");
grid on

%% PART - C -------------------------
% Comparison between Free Space and Two-Ray for one distance range

freq_op = 900e6;
transmit_power = 10;
tx_gain = 10^(0.5);   % in linear scale
rx_gain = 10^(0.3);   % in linear scale
loss_factor = 1;
tx_height = 40;
rx_height = 3;
lambda = 3e8 / freq_op;

% Free Space Received Power
const = ((4 * pi).^2) .* (d.^2);
PrFS = (transmit_power * tx_gain * rx_gain * lambda^2) ./ (loss_factor .* const);
path_loss_fs = 10 * log10(transmit_power) - 10 * log10(PrFS);

% Two-Ray Path Loss
path_loss_2ray = 40 * log10(d) - ...
    (10 * log10(tx_gain) + 10 * log10(rx_gain) + 20 * log10(tx_height) + 20 * log10(rx_height));

% Plotting Comparison
figure(5);
plot(d, path_loss_fs, 'b', 'LineWidth', 1.5);
hold on;
plot(d, path_loss_2ray, 'r--', 'LineWidth', 1.5);
ylabel("Path Loss (dB)");
xlabel("Distance (m)");
title("Path Loss: Free Space vs Two-Ray");
legend('Free Space', 'Two-Ray Model');
grid on



%Exp4
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


%Exp5
%Part A: Hata Model – Distance vs Path Loss

clc; clear; close all;

f = [700, 900]; % MHz
hte = 30; % m
hre = 1.5; % m
d = 1:2:20; % km

a_hre_urban = (1.1 * log10(f) - 0.7) * hre - (1.56 * log10(f) - 0.8);
a_hre_suburban = a_hre_urban - 2 * (log10(f / 28)).^2 - 5.4;
a_hre_rural = a_hre_urban - 4.78 * (log10(f)).^2 + 18.33 * log10(f) - 40.94;

PL_urban = zeros(length(f), length(d));
PL_large_city = zeros(length(f), length(d));
PL_suburban = zeros(length(f), length(d));
PL_rural = zeros(length(f), length(d));

for i = 1:length(f)
    PL_urban(i, :) = 69.55 + 26.16 * log10(f(i)) - 13.82 * log10(hte) - a_hre_urban(i) + (44.9 - 6.55 * log10(hte)) * log10(d);
    PL_large_city(i, :) = PL_urban(i, :) + 3;
    PL_suburban(i, :) = PL_urban(i, :) - 2 * (log10(f(i)/28)).^2 - 5.4;
    PL_rural(i, :) = PL_urban(i, :) - 4.78 * (log10(f(i))).^2 + 18.33 * log10(f(i)) - 40.94;
end

figure;
for i = 1:length(f)
    subplot(2,1,i);
    plot(d, PL_urban(i,:), 'b-o', d, PL_large_city(i,:), 'r-s', d, PL_suburban(i,:), 'c-^', d, PL_rural(i,:), 'm-d', 'LineWidth', 1.5);
    grid on;
    title(['Hata Path Loss at ', num2str(f(i)), ' MHz']);
    xlabel('Distance (km)');
    ylabel('Path Loss (dB)');
    legend('Urban', 'Large City', 'Suburban', 'Rural', 'Location', 'best');
end


%Part B: Hata Model – hte vs Path Loss

clc; clear; close all;

fc = 700; % MHz
hre = 1.5; % m
distances = [5, 10, 20]; % km
hte_values = 10:10:100; % m

alpha_hre = (1.1 * log10(fc) - 0.7) * hre - (1.56 * log10(fc) - 0.8);

figure; hold on;
for d = distances
    L50_urban = 69.55 + 26.16 * log10(fc) - 13.82 * log10(hte_values) - alpha_hre + (44.9 - 6.55 * log10(hte_values)) .* log10(d);
    plot(hte_values, L50_urban, '-o', 'DisplayName', ['d = ', num2str(d), ' km']);
end
xlabel('Transmitting Antenna Height (m)');
ylabel('Path Loss (dB)');
title('Path Loss vs hte (Urban)');
grid on; legend; hold off;

%Part C: Hata vs Two-Ray Model

clc; clear; close all;

f = 900; % MHz
hte = 30; % m
hre = 1.5; % m
d = 1:1:10; % km
c = 3e8;
lambda = c / (f * 1e6); % m

a_hre = (1.1 * log10(f) - 0.7) * hre - (1.56 * log10(f) - 0.8);
PL_Hata = 69.55 + 26.16 * log10(f) - 13.82 * log10(hte) - a_hre + (44.9 - 6.55 * log10(hte)) * log10(d);

d_m = d * 1e3;
PL_TwoRay = 40 * log10(d_m) - 20 * log10(hte * hre);

figure;
plot(d, PL_Hata, 'b-o', d, PL_TwoRay, 'r-s', 'LineWidth', 1.5);
xlabel('Distance (km)');
ylabel('Path Loss (dB)');
title('Hata vs Two-Ray Path Loss');
legend('Hata', 'Two-Ray'); grid on;


%Exp6
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
