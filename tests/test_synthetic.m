% expansion test for robustness
% test expansion needs of 0 to 6 expansions

clear 
clc

addpath(genpath("../src"))

nList = [1e5, 1e6, 1e7, 1e8, 5e8];
prob = 0.35;

t_hist = zeros(length(nList), 8);

for i = 1 : length(nList)
    n = nList(i);
    S = (rand(n, 1) < 0.5) * 2 - 1.0;
    S = S .* (rand(n, 1) < prob);
    S = S + 0.05 * randn(n, 1);
    fprintf("Case n = %g\n\n", n);

    rho = 1.0 / sqrt(n / 2);

    s_cpy = abs(S);
    timer1 = tic;
    [k, ~, ~, ~] = updateS_full(s_cpy, rho);
    t = toc(timer1);
    fprintf("Fullsort: Time -- %.6f,  k -- %d\n", t, k);
    t_hist(i, 1) = t;

    minWidth = 100;
    maxWidth = 5e4;
    alpha = 3;

    s_cpy = abs(S);
    kleft = max([k - 5e4, min([k - 100, ceil(k / 1.1)])]);
    kright = min([k + 5e4, max([k + 100, floor(k / 0.9)])]);
    k0 = randi([kleft + 1, kright - 1]);
    timer2 = tic;
    [k1, ~, ~, ~] = updateS_acc(s_cpy, rho, k0);
    t = toc(timer2);
    fprintf("Partialsort(one-hit): Time -- %.6f, k -- %d\n", t, k1);
    t_hist(i, 2) = t;
    
    for j = 1 : 6
        s_cpy = abs(S);
        kright = kleft;
        kleft = max([k - 5e4 - 1e5 * alpha^(j-1), min([k - 100 - 200 * alpha^(j-1), floor(k / (1.1 + 0.2 * alpha^(j-1)))])]);
        k0 = randi([kleft + 1, kright - 1]);
        k0 = max([k0, 1]);
        timer3 = tic;
        [k1, ~, ~, ~] = updateS_acc(s_cpy, rho, k0);
        t = toc(timer3);
        fprintf("Partialsort(resort%d): Time -- %.6f, k -- %d\n", j, t, k1);
        t_hist(i, j + 2) = t;
    end

    fprintf("\n\n");
end

respath = "../results";
if ~exist(respath, "dir")
    mkdir(respath);
end
filename = fullfile(respath, "synthetic_expand.txt");
fileID = fopen(filename, "w");

fprintf(fileID, "Results for Synthetic Robustness Test:\n");
fprintf(fileID, "Date: %s\n\n", datetime);
fprintf(fileID, "length   mode     expand     time\n\n");

for i = 1 : length(nList)
    n = nList(i);
    fprintf(fileID, "%d   full sort      0   %.6f\n", n, t_hist(i, 1));
    fprintf(fileID, "%d   partial sort   0   %.6f\n", n, t_hist(i, 2));
    for j = 1 : 6
        fprintf(fileID, "%d   partial sort   %d   %.6f\n", n, j, t_hist(i, j+2));
    end
    fprintf(fileID, "\n");
end

fprintf(fileID, "=========================");
fclose(fileID);
