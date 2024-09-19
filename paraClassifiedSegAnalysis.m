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
writetable(SegmentClassPercentage,'evaluationExport/SegmentClassPercentage.csv','Delimiter',',','QuoteStrings','All');

%% apply CNN image analysis to see difference from paraInitialImageAnalysis.m
% silenced for speed
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
% split data for cross validation (0.7 training, 0.3 validation)
cv = cvpartition(length(paraClassSeg.ID),'HoldOut',0.3);
valIndex = cv.test;
dataTrain = paraClassSeg(~valIndex,:);
% split validation to half (0.7 training, 0.15 validation, 0.15 testing)
dataValidate = paraClassSeg(valIndex,:);
cv = cvpartition(length(dataValidate.ID),'HoldOut',0.5);
testIndex = cv.test;
dataValidate = paraClassSeg(~testIndex,:);
dataTest = paraClassSeg(testIndex,:);

% decision tree for mean speed
meanSpeedTree = fitctree(dataTrain.AvgSpeed,dataTrain.ClassNum);
% calculate accuracy by 1-inaccuracy
accMSTree = 1-loss(meanSpeedTree,dataValidate.AvgSpeed,dataValidate.ClassNum);

% decision tree for msd
MSDTree = fitctree(dataTrain.MSDPrev,dataTrain.ClassNum);
accMSDTree = 1-loss(MSDTree,dataValidate.MSDPrev,dataValidate.ClassNum);

% decision tree for tsd
TSDTree = fitctree(dataTrain.TSD,dataTrain.ClassNum);
accTSDTree = 1-loss(TSDTree,dataValidate.TSD,dataValidate.ClassNum);

% decision tree for mean angle in degrees
meanAngleTree = fitctree(dataTrain.AvgDegTheta,dataTrain.ClassNum);
accMATree = 1-loss(meanAngleTree,dataValidate.AvgDegTheta,dataValidate.ClassNum);

%% Bagged decision trees
% bagged decision tree for mean speed
meanSpeedTree = TreeBagger(50,dataTrain.AvgSpeed,dataTrain.ClassNum, ...
    'Method','classification','OOBPrediction','on');

% plot out-of-bag classification error
figure('Name','OOB Classification Error','Position',[300 300 1200 300]);
subplot(1,4,1);
plot(oobError(meanSpeedTree));
title('Mean Speed OOB Error');
xlabel('Number of Trees Grown');
ylabel('Out-of-bag Classification Error');
% calculate accuracy by 1-inaccuracy
accMSBagTree = mean(1-error(meanSpeedTree,dataValidate.AvgSpeed,dataValidate.ClassNum));

% bagged decision tree for msd
MSDTree = TreeBagger(50,dataTrain.MSDPrev,dataTrain.ClassNum, ...
    'Method','classification','OOBPrediction','on');

% plot out-of-bag classification error
subplot(1,4,2);
plot(oobError(MSDTree));
title('MSD OOB Error');
xlabel('Number of Trees Grown');
ylabel('Out-of-bag Classification Error');
% calculate accuracy by 1-inaccuracy
accMSDBagTree = mean(1-error(MSDTree,dataValidate.MSDPrev,dataValidate.ClassNum));

% bagged decision tree for tsd
TSDTree = TreeBagger(50,dataTrain.TSD,dataTrain.ClassNum, ...
    'Method','classification','OOBPrediction','on');

% plot out-of-bag classification error
subplot(1,4,3);
plot(oobError(TSDTree));
title('TSD OOB Error');
xlabel('Number of Trees Grown');
ylabel('Out-of-bag Classification Error');
% calculate accuracy by 1-inaccuracy
accTSDBagTree = mean(1-error(TSDTree,dataValidate.TSD,dataValidate.ClassNum));

% bagged decision tree for mean angle
meanAngleTree = TreeBagger(50,dataTrain.AvgDegTheta,dataTrain.ClassNum, ...
    'Method','classification','OOBPrediction','on');

% plot out-of-bag classification error
subplot(1,4,4);
plot(oobError(meanAngleTree));
title('Mean Angle OOB Error');
xlabel('Number of Trees Grown');
ylabel('Out-of-bag Classification Error');
% calculate accuracy by 1-inaccuracy
accMABagTree = mean(1-error(meanAngleTree,dataValidate.AvgDegTheta,dataValidate.ClassNum));

%% Classification Neural Network
% make table into array
XTrain = table2array(dataTrain(:,["AvgSpeed" "MSDPrev" "TSD" "AvgDegTheta"]));
TTrain = categorical(dataTrain.ClassNum);
XValid = table2array(dataValidate(:,["AvgSpeed" "MSDPrev" "TSD" "AvgDegTheta"]));
TValid = categorical(dataValidate.ClassNum);
XTest = table2array(dataTest(:,["AvgSpeed" "MSDPrev" "TSD" "AvgDegTheta"]));
TTest = categorical(dataTest.ClassNum);

% network layers
layers = [
    % 4 features (mean speed, MSD, TSD, mean angle)
    featureInputLayer(4)
    % fully connected layer, 3 output classes
    fullyConnectedLayer(3)
    softmaxLayer];

% training options
options = trainingOptions("sgdm", ...
    ExecutionEnvironment="cpu", ...
    InitialLearnRate=0.001, ...
    MaxEpochs=60, ...
    ValidationData={XValid,TValid}, ...
    Metrics='accuracy', ...
    Plots="training-progress");

net = trainnet(XTrain,TTrain,layers,'crossentropy',options);

%% test network with allocated testing data
YTest = minibatchpredict(net,XTest);
YTest = onehotdecode(YTest,categories(TTest),2);
nnetAcc = mean(YTest == TTest);

figure;
confusionchart(TTest,YTest);

%% compile and export accuracy data
Method = {'Decision Tree';'Bagged Tree'};
MeanSpeedAccuracy = [accMSTree;accMSBagTree];
MSDAccuracy = [accMSDTree;accMSDBagTree];
TSDAccuracy = [accTSDTree;accTSDBagTree];
MeanAngleAccuracy = [accMATree;accMABagTree];

AccuracyTable = table(Method,MeanSpeedAccuracy,MSDAccuracy,TSDAccuracy,MeanAngleAccuracy);
writetable(AccuracyTable,'evaluationExport/AccuracyTable.csv','Delimiter',',','QuoteStrings','All');