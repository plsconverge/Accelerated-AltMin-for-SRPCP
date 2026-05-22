% perform AltMin on video data
% save the resulting matrices as images for visual check

clear
clc

data_path = "/data/ljq/VBM4D_rawRGB";
addpath(genpath("../src"))

datasets = {
    'billiards',         'M0009';
    'store',             'M0016';
    'tea set',           'M0017';
    'table tennis',      'M0008';
    'basketball',        'M0001';
    'windmill',          'M0005';
    'street',            'M0010';
    'flag',              'M0020'
};
idx = 5;
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
% params = struct("Lmod", "fullsvd", "Smod", "fullsort");
params = struct("Lmod", "fullsvd", "Smod", "partialsort");

[L, S, res] = altmin(Img, lambda, mu, params);

fprintf("Total running time: %.6f\n", res.time);
fprintf("Time of L-subproblem: %.6f\n", res.time_detail.L);
fprintf("Time of S-subproblem: %.6f\n", res.time_detail.S);

respath = "../results/pic";
if ~exist(respath, "dir")
    mkdir(respath);
end

for frame = [30, 60, 90]
    L_wri = uint8(255 * rescale(reshape(L(:, frame), n, m)));
    filename_L = fullfile(respath, sprintf("%s_frame%d_low_rank.png", video_name, frame));
    imwrite(L_wri, filename_L);

    S_wri = uint8(255 * rescale(reshape(S(:, frame), n, m)));
    filename_S = fullfile(respath, sprintf("%s_frame%d_sparse.png", video_name, frame));
    imwrite(S_wri, filename_S);

    recov = L + S;
    LS_wri = uint8(255 * rescale(reshape(recov(:, frame), n, m)));
    filename_LS = fullfile(respath, sprintf("%s_frame%d_recovery.png", video_name, frame));
    imwrite(LS_wri, filename_LS);
end
