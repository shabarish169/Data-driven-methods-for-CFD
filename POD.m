clc;
clear;
close all;

% Set size manually
n = 192;
m = 168;

% Load your image
myImage = imread("C:\Users\shaba\Desktop\jkk.jpg");
% Convert to grayscale if needed
if size(myImage,3) == 3
    myImage = rgb2gray(myImage);
end

% Resize
myImage = imresize(myImage, [n, m]);

% Flatten
testFace = double(myImage(:));

% Display
figure;
imshow(uint8(reshape(testFace, n, m)));
title('My Face');

%% Load the dataset
addpath("C:\Users\shaba\Downloads\DATA\DATA")
load allFaces.mat
  % make sure faces, nfaces, n, m are loaded

%% Prepare training data
trainingFaces = faces(:,1:sum(nfaces(1:36)));
avgFace = mean(trainingFaces,2);
X = trainingFaces - avgFace * ones(1, size(trainingFaces,2));

% Compute SVD
[U, S, V] = svd(X,'econ');

%% Choose a test image
% Pick a face from the training set or your own!
 % for example, pick the 10th face


%% Subtract mean
testFaceCentered = testFace - avgFace;

%% Project onto the basis faces
coefficients = U' * testFaceCentered;  % projections
disp(['Maximum eigenvalues = ' num2str(size(S,1))]);

%% Reconstruct the face using increasing number of basis vectors
reconstructedFace = zeros(size(testFace));
figure;
for k = [1, 5, 10, 20, 30, 40, 50, 100, 200,400, 500, 750,1000,1500,2000]  % number of eigenfaces
    reconstruction = avgFace + U(:,1:k) * coefficients(1:k);
    
    subplot(3,5,find([1, 5, 10, 20, 30, 40, 50, 100, 200,400, 500, 750,1000,1500,2000]==k))
    imshow(uint8(reshape(reconstruction,n,m)));
    title([num2str(k)]);
    
    drawnow;
end

