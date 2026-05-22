% video data test
% compare algorithms on 8 selected datasets

clear
clc

data_path = "/data/ljq/VBM4D_rawRGB";
addpath(genpath("../src"))

datasets = {
    'store',             'M0016';
    'tea set',           'M0017';
    'basketball',        'M0001';
    'street',            'M0010';
    'windmill',          'M0005';
    'table tennis',      'M0008';
    'billiards',         'M0009';
    'flag',              'M0020'
};
len = size(datasets, 1);

time_tab = zeros(len, 8);
sparse_hist = zeros(len, 1);

for idx = 1 : len

video_name = datasets{idx, 1};
folder = datasets{idx, 2};
fprintf("Processing video: %s (%s)...\n", video_name, folder);

im_folder = fullfile(data_path, folder);
files = dir(fullfile(im_folder, "*.png"));
num_frames = length(files);

info = imfinfo(fullfile(im_folder, files(1).name));
m = info.Width;
n = info.Height;
num_pixels = m * n;

Img = zeros(num_pixels, num_frames);

for i = 1 : num_frames
    im_rgb = imread(fullfile(im_folder, files(i).name));
    im_gray = im2double(rgb2gray(im_rgb));
    Img(:, i) = reshape(im_gray, [], 1);
end

lambda = 1.0 / sqrt(num_pixels);
mu = sqrt(num_frames / 2.0);
params1 = struct("Lmod", "fullsvd", "Smod", "fullsort", "disp", false);
params2 = struct("Lmod", "fullsvd", "Smod", "partialsort", "disp", false);

[~, ~, res1] = altmin(Img, lambda, mu, params1);

fprintf("Full sort:\n")
fprintf("Total running time: %.6f\n", res1.time);
fprintf("Time of L-subproblem: %.6f\n", res1.time_detail.L);
fprintf("Time of S-subproblem: %.6f\n", res1.time_detail.S);
fprintf("Time of sorting: %.6f\n", res1.time_detail.sort);

[~, S, res2] = altmin(Img, lambda, mu, params2);

fprintf("Partial sort:\n")
fprintf("Total running time: %.6f\n", res2.time);
fprintf("Time of L-subproblem: %.6f\n", res2.time_detail.L);
fprintf("Time of S-subproblem: %.6f\n", res2.time_detail.S);
fprintf("Time of sorting: %.6f\n\n", res2.time_detail.sort);

time_tab(idx, 1) = res1.time;
time_tab(idx, 2) = res1.time_detail.L;
time_tab(idx, 3) = res1.time_detail.S;
time_tab(idx, 4) = res1.time_detail.sort;
time_tab(idx, 5) = res2.time;
time_tab(idx, 6) = res2.time_detail.L;
time_tab(idx, 7) = res2.time_detail.S;
time_tab(idx, 8) = res2.time_detail.sort;

sparse_hist(idx) = nnz(S);

end


respath = "../results";
if ~exist(respath, "dir")
    mkdir(respath);
end
filename = fullfile(respath, "video.txt");
fileID = fopen(filename, "w");

fprintf(fileID, "Results for Video Test:\n");
fprintf(fileID, "Date: %s\n\n", datetime);
fprintf(fileID, "name    mode    time total    time L    time S    time sort    nnz\n\n");
for idx = 1 : len
    fprintf(fileID, "%s   full sort      %.6f   %.6f   %.6f   %.6f   %d\n", datasets{idx, 1}, time_tab(idx, 1), time_tab(idx, 2), time_tab(idx, 3), time_tab(idx, 4), sparse_hist(idx));
    fprintf(fileID, "%s   partial sort   %.6f   %.6f   %.6f   %.6f   %d\n\n", datasets{idx, 1}, time_tab(idx, 5), time_tab(idx, 6), time_tab(idx, 7), time_tab(idx, 8), sparse_hist(idx));
end
fprintf("====================================\n");

fclose(fileID);
