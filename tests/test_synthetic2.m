% synthetic matrix test for AltMin

clear 
clc

addpath(genpath("../src"))

mList = [3e3, 5e3, 1e4, 1e4, 2e4, 3e4];
nList = [3e3, 2e3, 1e3, 1e4, 5e3, 3e3];
r = 50;
spar = [0.1, 0.2, 0.4];
rng(42);

params1 = struct("Lmod", "fullsvd", "Smod", "fullsort", "tol", 1e-6, "disp", false);
params2 = struct("Lmod", "fullsvd", "Smod", "partialsort", "tol", 1e-6, "disp", false);

time_tot_ful = zeros(length(mList), length(spar));
time_L_ful = zeros(length(mList), length(spar));
time_S_ful = zeros(length(mList), length(spar));
time_sort_ful = zeros(length(mList), length(spar));
time_tot_par = zeros(length(mList), length(spar));
time_L_par = zeros(length(mList), length(spar));
time_S_par = zeros(length(mList), length(spar));
time_sort_par = zeros(length(mList), length(spar));

for i = 1 : length(mList)
    m = mList(i);
    n = nList(i);
    for j = 1 : length(spar)
        prob = spar(j);
        S = (rand(m, n) < prob) .* (2 * (rand(m, n) < 0.5) - 1) * 0.5;

        U = randn(m, r) / sqrt(m);
        V = randn(n, r) / sqrt(n);
        L = U * V';

        noise = 0.001 * randn(m, n);
        D = S + L + noise;

        lambda = 1 / sqrt(m);
        mu = sqrt(n / 2);

        [~, ~, res1] = altmin(D, lambda, mu, params1);
        time_tot_ful(i, j) = res1.time;
        time_L_ful(i, j) = res1.time_detail.L;
        time_S_ful(i, j) = res1.time_detail.S;
        time_sort_ful(i, j) = res1.time_detail.sort;

        [~, ~, res2] = altmin(D, lambda, mu, params2);
        time_tot_par(i, j) = res2.time;
        time_L_par(i, j) = res2.time_detail.L;
        time_S_par(i, j) = res2.time_detail.S;
        time_sort_par(i, j) = res2.time_detail.sort;
    end
end

fprintf("m      n    sparsity    mode       time       time L     time S     time sort\n\n")
for i = 1 : length(mList)
    for j = 1 : length(spar)
        fprintf("%d   %d   %.1f   full sort      %.6f   %.6f   %.6f   %.6f\n", mList(i), nList(i), spar(j), time_tot_ful(i, j), time_L_ful(i, j), time_S_ful(i, j), time_sort_ful(i, j));
        fprintf("%d   %d   %.1f   partial sort   %.6f   %.6f   %.6f   %.6f\n\n", mList(i), nList(i), spar(j), time_tot_par(i, j), time_L_par(i, j), time_S_par(i, j), time_sort_par(i, j));
    end
end


respath = "../results";
if ~exist(respath, "dir")
    mkdir(respath);
end
filename = fullfile(respath, "synthetic_unbalance.txt");
fileID = fopen(filename, "w");

fprintf(fileID, "Results for Synthetic Unbalanced Matrix Test:\n");
fprintf(fileID, "Date: %s\n\n", datetime);
fprintf(fileID, "m      n    sparsity    mode       time       time L     time S     time sort\n\n");
for i = 1 : length(mList)
    for j = 1 : length(spar)
        fprintf(fileID, "%d   %d   %.1f   full sort      %.6f   %.6f   %.6f   %.6f\n", mList(i), nList(i), spar(j), time_tot_ful(i, j), time_L_ful(i, j), time_S_ful(i, j), time_sort_ful(i, j));
        fprintf(fileID, "%d   %d   %.1f   partial sort   %.6f   %.6f   %.6f   %.6f\n\n", mList(i), nList(i), spar(j), time_tot_par(i, j), time_L_par(i, j), time_S_par(i, j), time_sort_par(i, j));
    end
end

fprintf(fileID, "======================================\n\n");
fclose(fileID);
