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
