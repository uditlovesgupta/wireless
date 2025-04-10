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
