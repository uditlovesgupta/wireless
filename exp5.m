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

