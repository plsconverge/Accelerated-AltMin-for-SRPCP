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
len = length(datasets);

size_list = zeros(len, 1);
num_iters = 50;
sparse_hist = zeros(len, num_iters);

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
params = struct("Lmod", "fullsvd", "Smod", "partialsort", "disp", false, "maxiter", num_iters);

[~, ~, res] = altmin(Img, lambda, mu, params);

size_list(idx) = num_pixels * num_frames;
sparse_hist(idx, :) = res.shist;

end


lineStyles = {'-', '--', ':', '-.', '-', '--', ':', '-.'};
lineColors = lines(8);
lineWidths = [1.5, 1.5, 1.5, 1.5, 1, 1, 1, 1];

fig = figure("Visible", "off", "Renderer", "painters");
hold on;

for i = 1 : len
    rnode = sum(sparse_hist(i, :) > 0);
    plot(1:rnode, sparse_hist(i, 1:rnode) / size_list(i), "LineStyle", lineStyles{i}, ...
        "Color", lineColors(i, :), "LineWidth", lineWidths(i));
end

lgd = legend(datasets(:, 1));
set(lgd, "Location", "northeast");

title("Change History of Sparsity");
xlabel("Iterations");
ylabel("Degree of Sparsity");

ax1 = gca;
ax1.Box = "off";
ax1.TickDir = "in";
ax1.XAxisLocation = "bottom";
ax1.YAxisLocation = "left";

xlim([0, 50]);
ylim([0, 0.6]);
xticks(0:10:50);
yticks(0:0.1:0.6);
ytickformat("%.1f")

ax2 = axes('Position', ax1.Position, ...
           'XAxisLocation', 'top', ...
           'YAxisLocation', 'right', ...
           'Color', 'none', ...
           'XLim', ax1.XLim, ...
           'YLim', ax1.YLim, ...
           'XTick', [], ...
           'YTick', [], ...
           'XColor', 'k', 'YColor', 'k', ...
           'Box', 'on');

axes(ax1);
ax1.LooseInset = ax1.TightInset;

set(gcf, "Color", "w");
set(gcf, "WindowState", "maximize");

respath = "../results";
if ~exist(respath, "dir")
    mkdir(respath);
end
drawnow;
% print(gcf, fullfile(respath, "fig_sparsity.pdf"), "-dpdf", "-bestfit");
exportgraphics(gcf, fullfile(respath, "fig_sparsity.png"), "Resolution", 300);

close(fig);
