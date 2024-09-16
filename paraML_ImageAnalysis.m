% Shido Nakajima
% Machine Learning Image Analysis using an adjusted version of previously built network for image classification.

clear;clc;close all;

%% import excel data, get index of movie and ID
% same as paraSort.m
paraData = readtable("data-SPZ-in-skin-to-analyze.xlsx");
paraData = sortrows(paraData,"movie");

% list index of where 'movie' value changes
movieIndex = ischange(paraData.movie);
movieIndex = find(movieIndex);

% list index of where 'PARASITEID' value changes
parasiteidIndex = zeros(length(movieIndex),1);
for i = 1:length(paraData.PARASITEID)
    parasiteidIndex(i) = str2double(extract(paraData.PARASITEID(i), digitsPattern(1,2)));
end
parasiteidIndex = ischange(parasiteidIndex);
parasiteidIndex = cat(1,1,find(parasiteidIndex));
parasiteidIndex = cat(1,parasiteidIndex,length(paraData.PARASITEID)+1);

%% plot each parasite as line graph and export as image
% run only if images dont exist (save time)
if ~exist('parasiteImages\','dir')
    % make folders for NINV and INV
    mkdir parasiteImages\INV\;
    mkdir parasiteImages\NINV\;
    
    for i = 1:length(parasiteidIndex)-1
        xy = [paraData.x_micron_(parasiteidIndex(i):parasiteidIndex(i+1)-1), ...
            paraData.y_micron_(parasiteidIndex(i) : parasiteidIndex(i+1)-1)];
        % plot parasite path, with red circle at start and blue at end
        plot(xy(:,1),xy(:,2),'-o','MarkerIndices',1:length(xy)-1:length(xy));
        axis off;
        hold on;
        plot(xy(1,1),xy(1,2),'o','Color','red');
        ax = gca;
        hold off;
    
        % save by NINV or INV
        id = char(paraData.PARASITEID(parasiteidIndex(i)));
        if (id(1) == 'N')
            exportgraphics(ax,'parasiteImages/NINV/'+string(paraData.PARASITEID(parasiteidIndex(i)))+'.jpg');
        else
            exportgraphics(ax,'parasiteImages/INV/'+string(paraData.PARASITEID(parasiteidIndex(i)))+'.jpg');
        end
    end
end

%% Image analysis using Deep Network Designer
% read and make image dataset of parasite locomotion
datasetPath = fullfile('parasiteImages/');
imds = imageDatastore(datasetPath, ...
    'IncludeSubfolders',true,'LabelSource','foldernames');

% split data into 80% training images and 20% validation images
[imdsTrain,imdsValidation] = splitEachLabel(imds,0.8,"randomized");

% make image dataset to 200x200
augTrain = augmentedImageDatastore([400, 400], imdsTrain);
augValidation = augmentedImageDatastore([400, 400], imdsValidation);

%% Define and train neural network for image classification
% define layers. 3 total training block conducted.
layers = [
    imageInputLayer([400 400 3])
	
    % 2-D convolutional filter that scans image (filter size, filter amount)
    convolution2dLayer(3,64)
    % normalize data to reduce sensitivity to network initialization
    batchNormalizationLayer
    % rectified linear unit (ReLU). Cuts any value below zero
    reluLayer
	
    % downsamples by taking maximum of each pooling region
    maxPooling2dLayer(2,Stride=2)
	
    convolution2dLayer(3,32)
    batchNormalizationLayer
    reluLayer
	
    maxPooling2dLayer(2,Stride=2)
	
    convolution2dLayer(3,64)
    batchNormalizationLayer
    reluLayer
	
    % applies weight and bias vector
    fullyConnectedLayer(2)
    % applies softmax function, making the matrix into probability distribution
    softmaxLayer];

% training options
options = trainingOptions("sgdm", ...
    MaxEpochs=16, ...
    ValidationData=augValidation, ...
    ValidationFrequency=2, ...
    Metrics="accuracy", ...
    Plots="training-progress");    

% train neural network defined above
net = trainnet(augTrain,layers,'crossentropy',options);
