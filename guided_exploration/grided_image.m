clc; clear all; close all;

% X  = imread('26.jpg');


addpath('/home/akshay/Desktop/CVC/Multi_target/Code/Multi-level_Sampling');
addpath('/home/akshay/Desktop/Data/Biopsy');
addpath('/home/akshay/Desktop/Data');
X = imread('stain_thumbnail.png');


[rows, columns, numberOfColorChannels] = size(X);

hold on;
imshow(X)

%% Grids
hold on;

%lineSpacing = 20; % Whatever you want.
stepSize = 50;

for row = 1 : stepSize : rows
    line([1, columns], [row, row], 'Color', 'b', 'LineWidth', 1);
end
for col = 1 : stepSize : columns
    line([col, col], [1, rows], 'Color', 'b', 'LineWidth', 1);
end

%% Rectangular ROI
rectangle('Position', [900 750 50 50], 'EdgeColor', 'r','LineWidth',3); %// draw rectangle on image
