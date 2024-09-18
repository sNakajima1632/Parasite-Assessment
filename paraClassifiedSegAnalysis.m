% Shido Nakajima
% Analysis of classified, segmented data. Goals are to do classification ML
% on the data, and determine which class each parasiteID falls under.

clear;clc;close all;

%% import excel data, get index for each parasiteID
% import original datatable
paraData = readtable("data-SPZ-in-skin-to-analyze.xlsx");
paraData = sortrows(paraData,"movie");

% import segmented classified csv file exported by paraTrajSegAnalysis.m
paraClassSeg = readtable("evaluationExport\analysisData.csv");

% index of ID from paraClassSeg
segIndex = zeros(length(paraClassSeg.ID),1);
for i = 1:length(paraClassSeg.ID)
    segIndex(i) = str2double(extract(paraClassSeg.ID(i),digitsPattern(1,2)));
end
segIndex = ischange(segIndex,'Threshold',0.01);
segIndex = cat(1,1,find(segIndex));
segIndex = cat(1,segIndex,length(paraClassSeg.ID)+1);

%% Move images exported by paraTrajSegAnalysis.m into parasiteImages\segmented30\ into parasiteImages\classifiedSeg30\1, \2, \3
mkdir parasiteImages\classifiedSeg30\1\;
mkdir parasiteImages\classifiedSeg30\2\;
mkdir parasiteImages\classifiedSeg30\3\;

% move images to class folder if there are images in the folder.
segImageCount = dir(fullfile('parasiteImages\segmented30\', '*.jpg'));
if size(segImageCount,1) ~= 0
    for i = 1:length(paraClassSeg.ID)
        currfile = "parasiteImages\segmented30\"+paraClassSeg.ID(i)+"seg"+string(paraClassSeg.SegmentNumber(i))+".jpg";
    
        % move file into each classification folder
        switch paraClassSeg.ClassNum(i)
            case 1
                destfolder = "parasiteImages\classifiedSeg30\1\";
                movefile(currfile, destfolder);
            case 2
                destfolder = "parasiteImages\classifiedSeg30\2\";
                movefile(currfile, destfolder);
            case 3
                destfolder = "parasiteImages\classifiedSeg30\3\";
                movefile(currfile, destfolder);
        end
    end
elseif exist('parasiteImages\segmented30\','dir')
    % delete empty folder
    rmdir parasiteImages\segmented30\;
end

%% determine which class(1-3) each parasite falls under
% initialization of values to put in new table
ID = unique(paraClassSeg.ID,'stable');
c1Percentage = [];
c2Percentage = [];
c3Percentage = [];

for i=1:length(segIndex)-1
    % set indices and number of segments
    currParaIndex = segIndex(i);
    nextParaIndex = segIndex(i+1)-1;
    segNumber = nextParaIndex+1-currParaIndex;

    % count how many segments were of each class
    c1Count = sum(paraClassSeg.ClassNum(currParaIndex:nextParaIndex) == 1);
    c2Count = sum(paraClassSeg.ClassNum(currParaIndex:nextParaIndex) == 2);
    c3Count = sum(paraClassSeg.ClassNum(currParaIndex:nextParaIndex) == 3);

    % precentage of segments in each class for current parasiteID
    c1Percentage = cat(1,c1Percentage,c1Count/segNumber);
    c2Percentage = cat(1,c2Percentage,c2Count/segNumber);
    c3Percentage = cat(1,c3Percentage,c3Count/segNumber);
end

% form a table of percentage of parasite's segments being class 1-3
SegmentClassPercentage = table(ID,c1Percentage,c2Percentage,c3Percentage);

%% apply image analysis to see difference from paraInitialImageAnalysis.m
%{
% read and make image dataset of parasite locomotion
datasetPath = fullfile('parasiteImages\classifiedSeg30\');
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
	
    % training block 1
    % 2-D convolutional filter that scans image (filter size, filter amount)
    convolution2dLayer(3,16)
    % normalize data to reduce sensitivity to network initialization
    batchNormalizationLayer
    % rectified linear unit (ReLU). Cuts any value below zero
    reluLayer
	
    % downsamples by taking maximum of each pooling region
    maxPooling2dLayer(2,Stride=2)
	
    % training block 2
    convolution2dLayer(3,32)
    batchNormalizationLayer
    reluLayer
	
    maxPooling2dLayer(2,Stride=2)
	
    % training block 3
    convolution2dLayer(3,64)
    batchNormalizationLayer
    reluLayer
	
    maxPooling2dLayer(2,Stride=2)
	
    % training block 4
    convolution2dLayer(3,128)
    batchNormalizationLayer
    reluLayer
	
    % applies weight and bias vector
    fullyConnectedLayer(3)
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
%}

%% machine learning analysis using decision tree
