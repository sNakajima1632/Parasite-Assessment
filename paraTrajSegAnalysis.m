% Shido Nakajima
% Trajectory segmentation(splitting) to analyze trajectory dissociated with
% NINV/INV due to the calculated values clustered together in one dataset 
% shows insignificance. 
% To make machine learning meaningful, decided to segment trajectory to create
% more data to feed into machine learning. Also apply clustering and 
% manual classification if needed, expecting a more obvious difference between datasets.

clear;clc;close all;

%% import excel data, get index of movie and ID
% same as paraSort.m
paraData = readtable("data-SPZ-in-skin-to-analyze.xlsx");
paraData = sortrows(paraData,"movie");

% list index of where 'PARASITEID' value changes
parasiteidIndex = zeros(length(paraData.PARASITEID),1);
for i = 1:length(paraData.PARASITEID)
    parasiteidIndex(i) = str2double(extract(paraData.PARASITEID(i), digitsPattern(1,2)));
end
parasiteidIndex = ischange(parasiteidIndex);
parasiteidIndex = cat(1,1,find(parasiteidIndex));
parasiteidIndex = cat(1,parasiteidIndex,length(paraData.PARASITEID)+1);
% starting index of current parasite on plot
currParaIndex = 1;
nextParaIndex = parasiteidIndex(2)-1;

%% make index list for segmented dataset
seg30index = [];
for i=1:length(parasiteidIndex)-1
    % add first index of each parasite into seg30index
    seg30index = cat(1,seg30index,parasiteidIndex(i));

    % expected index of next segment (initial+30)
    nextSegInd = parasiteidIndex(i)+30;
    % gap between expected next segment index and index of next parasite ID
    segGap = parasiteidIndex(i+1)-nextSegInd;

    while segGap >= 15 
        % loops to get next segment if gap is larger than 30
        if segGap > 30
            seg30index = cat(1,seg30index,nextSegInd);
            nextSegInd = nextSegInd+30;
            segGap = parasiteidIndex(i+1)-nextSegInd;
        % adds next segment index and ends if gap is between and including 15 and 30
        elseif 15<=segGap && segGap<=30
            seg30index = cat(1,seg30index,nextSegInd);
            nextSegInd = nextSegInd+30;
            segGap = parasiteidIndex(i+1)-nextSegInd;
        end
        % ignores gap less than 15 including negative values
    end
end
seg30index = cat(1,seg30index,parasiteidIndex(end));

%% average/instantaneous speeds, mean sq displacement, trajectory step dispersion initiation
posX = paraData.x_micron_;
posY = paraData.y_micron_;
instSpeed = paraData.y_micron_;

SegmentNumber = {};
ID = {};
AvgSpeed = [];
MSDPrev = [];
TSD = [];
AvgCosTheta = [];
AvgDegTheta = [];
    
%% make new directory if needed
segImageDir = 'parasiteImages\segmented30\';
if ~exist(segImageDir,'dir')
    % make folders for NINV and INV
    mkdir parasiteImages\segmented30\;
end

% counter to indicate which segment of a parasite is shown, primarily for
% trajectory image naming convention
segmentCounter = 1;

%% loop to calculate all values of segmented trajectory for plotting (i=1:500)
% also plot and export trajectory image
for i = 1:length(seg30index)-1
    % index i integrated into parasiteidIndex for cleanliness
    inow = seg30index(i);
    ilast = seg30index(i+1)-1;
    currParaIndex = inow;

    % posXY and speedXY only feature coordinates of current parasiteID
    posXY = [posX(currParaIndex:ilast),posY(currParaIndex:ilast)];
    speedXY = gradient(posXY')'./gradient(paraData.t_sec_(inow:ilast));
    
    % reset segmentCounter with new parasiteID
    if segmentCounter > 1 && ismember(seg30index(i),parasiteidIndex)
        segmentCounter = 1;
    end
    % count number of images in the folder, and skip if images already exists to save time.
    segImageCount = dir(fullfile(segImageDir, '*.jpg'));
    if size(segImageCount,1) < length(seg30index)-1
        % plot and save segmented trajectory
        plotExport(posXY,paraData,seg30index,i,segmentCounter);
    end

    % list of instantaneous velocities for all points
    instSpeed(inow:ilast) = hypot(speedXY(:,1),speedXY(:,2));
    AvgSpeed = cat(1,AvgSpeed,mean(instSpeed(inow:ilast)));

    % save parasiteID for saving in new table
    ID = cat(1,ID,paraData.PARASITEID(inow));

    % save MSD by NINV or INV
    id = char(paraData.PARASITEID(inow));

    % MSD calculation
    msdp = mean(sum(diff(posXY).^2,2));
    MSDPrev = cat(1, MSDPrev, msdp);

    % TSD calculation
    stepLength = sqrt(sum(diff(posXY).^2,2));
    TSD = cat(1,TSD,sqrt(mean((stepLength(:)-mean(stepLength)).^2)));

    % angle between trajectories calculation
    cosTheta = [];
    degTheta = [];
    % loop to calculate angles
    for j=2:length(posXY)-1
        % vector of trajectory forming an angle (x2-x1,y2-y1)
        u=diff(posXY(j-1:j,:));
        v=diff(posXY(j:j+1,:)).*-1;

        % calculate mean angle change from vectors and convert to degrees
        cosTheta = cat(1,cosTheta,mean(max(min(dot(u,v)/(norm(u)*norm(v)),1),-1)));
        degTheta = cat(1,degTheta,mean(real(acosd(cosTheta))));
    end
    % take mean of all angles within current segment
    AvgCosTheta = cat(1,AvgCosTheta,mean(cosTheta));
    AvgDegTheta = cat(1,AvgDegTheta,mean(degTheta));

    % update segment number
    SegmentNumber = cat(1,SegmentNumber,segmentCounter);
    segmentCounter = segmentCounter+1;
end

% create table
MSDTable = table(ID,SegmentNumber,AvgSpeed,MSDPrev,TSD,AvgCosTheta,AvgDegTheta);

% create correlation matrix between mean speed, msdprev, msdorig, and tsd
rho = partialcorr([MSDTable.AvgSpeed MSDTable.MSDPrev MSDTable.TSD MSDTable.AvgDegTheta]);
% remove upper half
rho(logical(triu(ones(size(rho)),1))) = NaN;

%% plot data
% swarm plot for each data independently
f = figure('Name','Parasite Analysis','Position',[100 100 1200 800]);
axNum=tiledlayout('horizontal','Padding','compact');
% average speed plot
nexttile;
x=ones(length(MSDTable.AvgSpeed));
y=MSDTable.AvgSpeed;
swarmchart(x,y,10);
title('Mean Speed');

% MSD plot
nexttile;
y=MSDTable.MSDPrev;
swarmchart(x,y,10);
title('MSD');

% TSD from initial point plot
nexttile;
y=MSDTable.TSD;
swarmchart(x,y,10);
title('TSD');

% Mean angle between trajectories
nexttile;
y=MSDTable.AvgDegTheta;
swarmchart(x,y,10);
title('Mean Angle');

% plot heatmap of correlation coefficient between speed, MSD, and TSD
f2 = figure('Name','Pearson Correlation Coefficient between calculated values');
h = heatmap(rho,...
    'MissingDataColor','k', ...
    'XDisplayLabels',["Mean Speed", "MSD", "TSD", "Mean Angle"], ...
    'YDisplayLabels',["Mean Speed", "MSD", "TSD", "Mean Angle"]);

% 6 plots for: Mean Speed vs MSD, Mean Speed vs TSD, Mean Speed vs Mean Angle,
% MSD vs TSD, MSD vs Mean Angle, TSD vs Mean Angle
f3 = figure('Name','Analysis data combination plot','Position',[200 200 1400 800]);
% AvgSpeed vs MSDPrev plot
axComb=subplot(2,3,1);
x=MSDTable.AvgSpeed;
y=MSDTable.MSDPrev;
plot(axComb,x,y,'.');
xlabel(axComb,'AvgSpeed');
ylabel(axComb,'MSD');
title('AvgSpeed vs MSD');

% AvgSpeed vs TSD
axComb=subplot(2,3,2);
y=MSDTable.TSD;
plot(axComb,x,y,'.');
xlabel(axComb,'AvgSpeed');
ylabel(axComb,'TSD');
title('AvgSpeed vs TSD');

% AvgSpeed vs Mean Angle
axComb=subplot(2,3,3);
y=MSDTable.AvgDegTheta;
plot(axComb,x,y,'.');
xlabel(axComb,'AvgSpeed');
ylabel(axComb,'Mean Angle');
title('AvgSpeed vs Mean Angle');

% MSD vs TSD
axComb=subplot(2,3,4);
x=MSDTable.MSDPrev;
y=MSDTable.TSD;
plot(axComb,x,y,'.');
xlabel(axComb,'MSD');
ylabel(axComb,'TSD');
title('MSD vs TSD');

% MSD vs Mean Angle
axComb=subplot(2,3,5);
x=MSDTable.MSDPrev;
y=MSDTable.AvgDegTheta;
plot(axComb,x,y,'.');
xlabel(axComb,'MSD');
ylabel(axComb,'Mean Angle');
title('MSD vs Mean Angle');

% TSD vs Mean Angle
axComb=subplot(2,3,6);
x=MSDTable.TSD;
y=MSDTable.AvgDegTheta;
plot(axComb,x,y,'.');
xlabel(axComb,'TSD');
ylabel(axComb,'Mean Angle');
title('TSD vs Mean Angle');

%% Use Machine Learning (clustering) to evaluate data
% Hierarchical Clustering evaluation on single data to observe optimal number of clusters
% Dropped due to always resulting in maximum number of clusters
%{
eva = evalclusters(AvgSpeed,'linkage','CalinskiHarabasz','KList',1:125);
disp("Mean speed Hierarchical Clustering result");
disp(eva);
eva = evalclusters(MSDPrev,'linkage','CalinskiHarabasz','KList',1:125);
disp("MSD Hierarchical Clustering result");
disp(eva);
eva = evalclusters(MSDOrig,'linkage','CalinskiHarabasz','KList',1:125);
disp("MSD with reference at t=1 Hierarchical Clustering result");
disp(eva);
eva = evalclusters(TSD,'linkage','CalinskiHarabasz','KList',1:125);
disp("TSD Hierarchical Clustering result");
disp(eva);

% k-means clustering evaluation to observe optimal number of clusters
eva = evalclusters(AvgSpeed,'kmeans','CalinskiHarabasz','KList',1:125);
disp("Mean speed k-Means Clustering result");
disp(eva);
eva = evalclusters(MSDPrev,'kmeans','CalinskiHarabasz','KList',1:125);
disp("MSD k-Means Clustering result");
disp(eva);
eva = evalclusters(MSDOrig,'kmeans','CalinskiHarabasz','KList',1:125);
disp("MSD with reference at t=1 k-Means Clustering result");
disp(eva);
eva = evalclusters(TSD,'kmeans','CalinskiHarabasz','KList',1:125);
disp("TSD k-Means Clustering result");
disp(eva);
%}

% k-means clustering evaluation for combination using CalinskiHarabasz for speed
% all possible evaluation methods were used in the report
eva = evalclusters([AvgSpeed,MSDPrev],'kmeans','silhouette','KList',1:10);
disp("Mean speed vs MSD");
disp(eva);
eva = evalclusters([AvgSpeed,TSD],'kmeans','silhouette','KList',1:10);
disp("Mean speed vs TSD");
disp(eva);
eva = evalclusters([AvgSpeed,AvgDegTheta],'kmeans','silhouette','KList',1:10);
disp("Mean speed vs Mean Angle");
disp(eva);
eva = evalclusters([MSDPrev,TSD],'kmeans','silhouette','KList',1:10);
disp("MSD vs TSD");
disp(eva);
eva = evalclusters([MSDPrev,AvgDegTheta],'kmeans','silhouette','KList',1:10);
disp("MSD vs Mean Angle");
disp(eva);
eva = evalclusters([TSD,AvgDegTheta],'kmeans','silhouette','KList',1:10);
disp("TSD vs Mean Angle");
disp(eva);

%% plot kmeans colorcoded scatter plot of significant values from above
% MSD vs Mean angle with k=3
idc = kmeans([MSDPrev,AvgDegTheta],3);
figure('Name','k-means of MSD vs Mean Angle','Position',[500 250 1000 500]);
subplot(1,2,1);
gscatter(MSDPrev,AvgDegTheta,idc,'bgm');
legend('Cluster 1','Cluster 2','Cluster 3');
% plot evaluation
subplot(1,2,2);
plot(evalclusters([MSDPrev,AvgDegTheta],'kmeans','silhouette','KList',1:10));

% plot for k=7 and k=8 that were significant. Silenced due to non-stable
% evaluation method and result.
%{
% MSD vs Mean angle with k=7
idc = kmeans([MSDPrev,AvgDegTheta],7);
figure('Name','k-means of MSD vs Mean Angle','Position',[500 250 1000 500]);
subplot(1,2,1);
gscatter(MSDPrev,AvgDegTheta,idc,'rgbcmyk');
legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5','Cluster 6','Cluster 7');
% plot evaluation
subplot(1,2,2);
plot(evalclusters([MSDPrev,AvgDegTheta],'kmeans','CalinskiHarabasz','KList',1:10));

% Mean speed vs Mean anagle with k=7
idc = kmeans([AvgSpeed,AvgDegTheta],7);
figure('Name','k-means of Mean Speed vs Mean Angle','Position',[500 250 1000 500]);
subplot(1,2,1);
gscatter(AvgSpeed,AvgDegTheta,idc,'rgbcmyk');
legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5','Cluster 6','Cluster 7');
% plot evaluation
subplot(1,2,2);
plot(evalclusters([AvgSpeed,AvgDegTheta],'kmeans','DaviesBouldin','KList',1:10));

% MSD vs TSD with k=8
idc = kmeans([MSDPrev,TSD],8);
figure('Name','k-means of MSD vs TSD','Position',[500 250 1000 500]);
subplot(1,2,1);
gscatter(MSDPrev,TSD,idc,'rgbcmykw');
legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5','Cluster 6','Cluster 7','Cluster 8');
set(gca,'Color',[0.75 0.75 0.75]);
% plot evaluation
subplot(1,2,2);
plot(evalclusters([MSDPrev,TSD],'kmeans','DaviesBouldin','KList',1:10));
%}

%% classification id (1-3) based on k-means evaluation above added to MSDTable
ClassNum = idc;
Tcn = table(ClassNum);
MSDTable = [MSDTable, Tcn];
MSDTable = movevars(MSDTable,'ClassNum','Before','ID');

% make folder and export table for easier access. Comment out if wanted
mkdir evaluationExport\;
writetable(MSDTable,'evaluationExport/analysisData.csv','Delimiter',',','QuoteStrings','All');

%% function to plot segmented trajectory and save as jpg file
function plotExport(posXY,paraData,seg30index,i,segmentCounter)
    % plot parasite path and export
    plot(posXY(:,1),posXY(:,2),'-o','MarkerIndices',1:length(posXY)-1:length(posXY));
    axis off;
    hold on;
    plot(posXY(1,1),posXY(1,2),'o','Color','red');
    ax = gca;
    hold off;
    exportgraphics(ax,'parasiteImages/segmented30/'+string(paraData.PARASITEID(seg30index(i)))+'seg'+string(segmentCounter)+'.jpg');
end