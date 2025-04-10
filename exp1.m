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
