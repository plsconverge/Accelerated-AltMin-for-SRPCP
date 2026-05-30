% test AltMin with accelerated UpdateS on video data
% experiment on a single video for an early check

clear
clc

profile on

data_path = "../data/VBM4D_rawRGB";
addpath(genpath("../src"))

datasets = {
    'billiards',         'M0009';
    'store',             'M0016';
    'tea set',           'M0017';
    'table tennis',      'M0008';
    'basketball player', 'M0001';
    'windmill',          'M0005';
    'street',            'M0010';
    'flag',              'M0020'
};
idx = 3;
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

profile off
